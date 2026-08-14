/* read_filterbank.c
 *
 * Filterbank-input adapter for psrfits_subband: implements
 * filterbank_open()/filterbank_read_subint()/filterbank_read_part_DATA(),
 * drop-in replacements for psrfits_open()/psrfits_read_subint()/
 * psrfits_read_part_DATA() that populate the *same* struct psrfits fields
 * from a SIGPROC filterbank file instead of a PSRFITS file, so the rest of
 * psrfits_subband.c's subbanding/dedispersion/output-writing logic (which
 * only touches struct psrfits, never the input file format directly) runs
 * completely unmodified.
 *
 * The sigprocfb struct, and get_string()/strings_equal()/
 * read_filterbank_header(), are adapted from PRESTO's src/sigproc_fb.c
 * (Scott Ransom et al.), which itself notes "Much of this has been ripped
 * out of SIGPROC and then slightly modified. Thanks Dunc!". Reproduced
 * here (rather than linked from libpresto) because they are declared
 * `static` in PRESTO and not exported.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psrfits.h"

extern void unpack_2bit_to_8bit_unsigned(unsigned char *indata,
                                         unsigned char *outdata, int N);
extern void unpack_4bit_to_8bit_unsigned(unsigned char *indata,
                                         unsigned char *outdata, int N);

typedef struct SIGPROCFB {
    char inpfile[80];
    char source_name[80];
    double tstart;
    double tsamp;
    double src_raj;
    double src_dej;
    double az_start;
    double za_start;
    double fch1;
    double foff;
    double refdm;
    int machine_id;
    int telescope_id;
    int nchans;
    int nbits;
    int nifs;
    int nbeams;
    int ibeam;
    int sumifs;
    int signedints;
} sigprocfb;

/* Convert an MJD to a PSRFITS-style "YYYY-MM-DDTHH:MM:SS.SSS" DATE-OBS
 * string. Uses the standard Fliegel & Van Flandern Julian Day Number ->
 * Gregorian calendar algorithm (integer arithmetic, proleptic Gregorian).
 * Hand- and astropy-verified against MJD 54700.00417824074 ->
 * "2008-08-22T00:06:01.000". */
static void mjd_to_date_obs(double mjd, char *out, size_t outlen)
{
    long long jdn = (long long)floor(mjd) + 2400001LL;
    long long a = jdn + 32044;
    long long b = (4 * a + 3) / 146097;
    long long c = a - (146097 * b) / 4;
    long long d = (4 * c + 3) / 1461;
    long long e = c - (1461 * d) / 4;
    long long m = (5 * e + 2) / 153;
    int day = (int)(e - (153 * m + 2) / 5 + 1);
    int month = (int)(m + 3 - 12 * (m / 10));
    int year = (int)(100 * b + d - 4800 + m / 10);

    double frac = mjd - floor(mjd);
    double secs = frac * 86400.0;
    int hh = (int)(secs / 3600.0);
    int mm = (int)((secs - hh * 3600.0) / 60.0);
    double ss = secs - hh * 3600.0 - mm * 60.0;

    snprintf(out, outlen, "%04d-%02d-%02dT%02d:%02d:%06.3f",
             year, month, day, hh, mm, ss);
}

static size_t chkfread_local(void *data, size_t type, size_t number, FILE *stream)
{
    size_t num = fread(data, type, number, stream);
    if (num != number && ferror(stream)) {
        perror("\nError in read_filterbank.c's chkfread_local()");
        exit(-1);
    }
    return num;
}

static void get_string(FILE *inputfile, int *nbytes, char string[])
{
    int nchar;
    strcpy(string, "ERROR");
    chkfread_local(&nchar, sizeof(int), 1, inputfile);
    *nbytes = sizeof(int);
    if (feof(inputfile)) exit(0);
    if (nchar > 80 || nchar < 1)
        return;
    chkfread_local(string, nchar, 1, inputfile);
    string[nchar] = '\0';
    *nbytes += nchar;
}

static int strings_equal(char *string1, char *string2)
{
    return !strcmp(string1, string2);
}

/* Adapted from PRESTO's sigproc_fb.c:read_filterbank_header(). Returns the
 * header length in bytes (via *headerlen), or exits on malformed input. */
static int read_filterbank_header(sigprocfb *fb, FILE *inputfile, long long *headerlen)
{
    char string[80];
    int itmp, nbytes = 0, totalbytes;
    int expecting_rawdatafile = 0, expecting_source_name = 0;
    int barycentric, pulsarcentric;

    get_string(inputfile, &nbytes, string);
    if (!strings_equal(string, "HEADER_START")) {
        rewind(inputfile);
        return 0;
    }
    totalbytes = nbytes;

    fb->ibeam = 1;
    fb->signedints = 0;
    fb->sumifs = 1;
    fb->nifs = 1;
    fb->nbeams = 1;

    while (1) {
        get_string(inputfile, &nbytes, string);
        if (strings_equal(string, "HEADER_END"))
            break;
        totalbytes += nbytes;
        if (strings_equal(string, "rawdatafile")) {
            expecting_rawdatafile = 1;
        } else if (strings_equal(string, "source_name")) {
            expecting_source_name = 1;
        } else if (strings_equal(string, "az_start")) {
            chkfread_local(&(fb->az_start), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "za_start")) {
            chkfread_local(&(fb->za_start), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "src_raj")) {
            chkfread_local(&(fb->src_raj), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "src_dej")) {
            chkfread_local(&(fb->src_dej), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "tstart")) {
            chkfread_local(&(fb->tstart), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "tsamp")) {
            chkfread_local(&(fb->tsamp), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "fch1")) {
            chkfread_local(&(fb->fch1), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "foff")) {
            chkfread_local(&(fb->foff), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "refdm")) {
            chkfread_local(&(fb->refdm), sizeof(double), 1, inputfile);
            totalbytes += sizeof(double);
        } else if (strings_equal(string, "nchans")) {
            chkfread_local(&(fb->nchans), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "telescope_id")) {
            chkfread_local(&(fb->telescope_id), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "machine_id")) {
            chkfread_local(&(fb->machine_id), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "data_type")) {
            chkfread_local(&itmp, sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "nbits")) {
            chkfread_local(&(fb->nbits), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "barycentric")) {
            chkfread_local(&barycentric, sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "pulsarcentric")) {
            chkfread_local(&pulsarcentric, sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "nsamples")) {
            chkfread_local(&itmp, sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "nifs")) {
            chkfread_local(&(fb->nifs), sizeof(int), 1, inputfile);
            if (fb->nifs > 1) fb->sumifs = 0;
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "nbeams")) {
            chkfread_local(&(fb->nbeams), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "ibeam")) {
            chkfread_local(&(fb->ibeam), sizeof(int), 1, inputfile);
            totalbytes += sizeof(int);
        } else if (strings_equal(string, "signed")) {
            char tmp;
            chkfread_local(&tmp, sizeof(char), 1, inputfile);
            fb->signedints = tmp;
            totalbytes += sizeof(char);
        } else if (expecting_rawdatafile) {
            strcpy(fb->inpfile, string);
            expecting_rawdatafile = 0;
        } else if (expecting_source_name) {
            strcpy(fb->source_name, string);
            expecting_source_name = 0;
        } else {
            fprintf(stderr,
                    "ERROR: read_filterbank_header - unknown parameter: %s\n",
                    string);
            exit(1);
        }
    }
    totalbytes += nbytes;  // HEADER_END itself
    *headerlen = totalbytes;
    return 1;
}

/* --- Minimal telescope/machine name lookup, matching presto.sigproc's
 * telescope_ids/machine_ids -- purely informational (written into the
 * output PSRFITS header), not used by the subbanding/dedispersion math. */
static const char *telescope_name(int id)
{
    switch (id) {
        case 0: return "Fake";
        case 1: return "Arecibo";
        case 2: return "Ooty";
        case 3: return "Nancay";
        case 4: return "Parkes";
        case 5: return "Jodrell";
        case 6: return "GBT";
        case 7: return "GMRT";
        case 8: return "Effelsberg";
        case 9: return "ATA";
        case 10: return "SRT";
        case 11: return "LOFAR";
        case 12: return "VLA";
        case 20: return "CHIME";
        case 21: return "FAST";
        case 30: return "MWA";
        case 64: return "MeerKAT";
        case 65: return "KAT-7";
        default: return "Unknown";
    }
}

static const char *machine_name(int id)
{
    switch (id) {
        case 0: return "FAKE";
        case 1: return "PSPM";
        case 2: return "WAPP";
        case 3: return "AOFTM";
        case 4: return "BCPM1";
        case 5: return "OOTY";
        case 6: return "SCAMP";
        case 7: return "SPIGOT";
        case 11: return "BG/P";
        case 12: return "PDEV";
        default: return "FILTERBANK";
    }
}

/* --- Module state for the (possibly multi-file) filterbank sequence.
 * The whole sequence is treated as one continuous virtual sample stream:
 * unlike PSRFITS (where NSBLK guarantees every row divides evenly within
 * a file), a filterbank file's sample count need not be a multiple of our
 * chosen block size, so a single block read can legitimately span two
 * files. fb_read_global() below is the one place that crosses that
 * boundary; everything else just asks it for samples starting at a given
 * position in the concatenated stream. --- */
static FILE *fb_file = NULL;
static int fb_open_idx = -1;                /* index into fb_filenames of fb_file, -1 = none open */
static char **fb_filenames = NULL;          /* alias of pf->filenames, not owned */
static int fb_numfiles = 0;
static long long *fb_file_headerlen = NULL; /* per-file header length, bytes */
static long long *fb_file_nsamples = NULL;  /* per-file sample count */
static long long *fb_file_start = NULL;     /* prefix sums; size fb_numfiles+1 */
static sigprocfb fb_hdr;                    /* file 0's header -- canonical format */
static long long fb_bytes_per_sample_row = 0;  /* nchans*nifs*nbits/8 */
static long long fb_global_pos = 0;         /* "committed" position, in spectra, across all files */
static long long fb_total_samples = 0;      /* sum across all files */


static void fb_ensure_open(int idx)
{
    if (fb_open_idx == idx && fb_file != NULL)
        return;
    if (fb_file) {
        fclose(fb_file);
        fb_file = NULL;
    }
    fb_file = fopen(fb_filenames[idx], "rb");
    if (fb_file == NULL) {
        fprintf(stderr, "Error: could not open filterbank file '%s'\n",
                fb_filenames[idx]);
        exit(1);
    }
    fb_open_idx = idx;
}


static int fb_file_for_global(long long global_pos)
{
    int i;
    for (i = 0; i < fb_numfiles; i++)
        if (global_pos < fb_file_start[i + 1])
            return i;
    return fb_numfiles;  /* past the end of the last file */
}


/* Read `nspec` sample-rows starting at global sample index `global_pos`
 * (0-based, across the whole concatenated file sequence) into `buf`,
 * transparently crossing file boundaries as needed. Returns the number of
 * rows actually read (< nspec at true end-of-sequence or a short read),
 * matching fread()'s partial-read convention. Purely a function of
 * `global_pos` -- no side effects on any "current position" state -- so
 * it works equally for a committed read (advance afterward) or a
 * non-committing peek (don't). */
static long long fb_read_global(long long global_pos, long long nspec, unsigned char *buf)
{
    long long remaining = nspec, pos = global_pos;
    unsigned char *ptr = buf;

    while (remaining > 0) {
        int idx = fb_file_for_global(pos);
        long long within_file, avail_in_file, this_read, nread;
        if (idx >= fb_numfiles)
            break;
        fb_ensure_open(idx);
        within_file = pos - fb_file_start[idx];
        avail_in_file = fb_file_nsamples[idx] - within_file;
        this_read = (remaining < avail_in_file) ? remaining : avail_in_file;
        fseeko(fb_file, fb_file_headerlen[idx] + within_file * fb_bytes_per_sample_row,
               SEEK_SET);
        nread = fread(ptr, fb_bytes_per_sample_row, this_read, fb_file);
        ptr += nread * fb_bytes_per_sample_row;
        pos += nread;
        remaining -= nread;
        if (nread < this_read)
            break;  /* short read / unexpected EOF within a file */
    }
    return nspec - remaining;
}


int filterbank_open(struct psrfits *pf)
{
    struct hdrinfo *hdr = &(pf->hdr);
    int i;

    if (pf->numfiles == 0) {
        fprintf(stderr, "Error: -filterbank requires at least one input file.\n");
        pf->status = 1;
        return pf->status;
    }

    fb_numfiles = pf->numfiles;
    fb_filenames = pf->filenames;
    fb_file_headerlen = (long long *)malloc(sizeof(long long) * fb_numfiles);
    fb_file_nsamples = (long long *)malloc(sizeof(long long) * fb_numfiles);
    fb_file_start = (long long *)malloc(sizeof(long long) * (fb_numfiles + 1));

    for (i = 0; i < fb_numfiles; i++) {
        sigprocfb this_hdr;
        long long this_headerlen, filelen, this_bytes_per_row;
        FILE *f = fopen(fb_filenames[i], "rb");
        if (f == NULL) {
            fprintf(stderr, "Error: could not open filterbank file '%s'\n",
                    fb_filenames[i]);
            pf->status = 1;
            return pf->status;
        }
        if (!read_filterbank_header(&this_hdr, f, &this_headerlen)) {
            fprintf(stderr,
                    "Error: '%s' does not look like a SIGPROC filterbank file "
                    "(no HEADER_START).\n", fb_filenames[i]);
            fclose(f);
            pf->status = 1;
            return pf->status;
        }
        if (i == 0) {
            fb_hdr = this_hdr;  /* canonical format reference for the whole sequence */
        } else {
            /* Structural mismatches would misalign every read after this
             * file (wrong byte stride/channel count) -- hard error. */
            if (this_hdr.nchans != fb_hdr.nchans || this_hdr.nbits != fb_hdr.nbits ||
                this_hdr.nifs != fb_hdr.nifs) {
                fprintf(stderr,
                        "Error: '%s' (nchans=%d nbits=%d nifs=%d) doesn't match "
                        "'%s' (nchans=%d nbits=%d nifs=%d) -- can't treat these "
                        "as one contiguous observation.\n",
                        fb_filenames[i], this_hdr.nchans, this_hdr.nbits, this_hdr.nifs,
                        fb_filenames[0], fb_hdr.nchans, fb_hdr.nbits, fb_hdr.nifs);
                fclose(f);
                pf->status = 1;
                return pf->status;
            }
            /* Frequency-labeling mismatches don't affect the byte layout,
             * so just warn (mirroring PSRFITS's own soft warnings for
             * mismatched TELESCOP etc. across a multi-file set). */
            if (fabs(this_hdr.fch1 - fb_hdr.fch1) > 1e-6 ||
                fabs(this_hdr.foff - fb_hdr.foff) > 1e-9) {
                fprintf(stderr,
                        "Warning: '%s' (fch1=%.6f foff=%.6f) doesn't match '%s' "
                        "(fch1=%.6f foff=%.6f); using file 0's frequencies for "
                        "the whole sequence.\n",
                        fb_filenames[i], this_hdr.fch1, this_hdr.foff,
                        fb_filenames[0], fb_hdr.fch1, fb_hdr.foff);
            }
        }

        this_bytes_per_row = ((long long)this_hdr.nchans * this_hdr.nifs * this_hdr.nbits) / 8;
        fseeko(f, 0, SEEK_END);
        filelen = ftello(f);
        fclose(f);

        fb_file_headerlen[i] = this_headerlen;
        fb_file_nsamples[i] = (filelen - this_headerlen) / this_bytes_per_row;
        printf("Found filterbank file '%s': %lld spectra\n",
               fb_filenames[i], (long long)fb_file_nsamples[i]);
    }

    fb_file_start[0] = 0;
    for (i = 0; i < fb_numfiles; i++)
        fb_file_start[i + 1] = fb_file_start[i] + fb_file_nsamples[i];
    fb_total_samples = fb_file_start[fb_numfiles];
    strncpy(pf->filename, fb_filenames[0], 200);
    printf("Opened %d filterbank file(s), %lld total spectra.\n",
           fb_numfiles, (long long)fb_total_samples);

    if (fb_hdr.nbits != 8) {
        fprintf(stderr,
                "Warning: filterbank nbits=%d; this reader has only been "
                "validated against 8-bit data. Proceed with caution.\n",
                fb_hdr.nbits);
    }

    /* --- Populate struct hdrinfo, mirroring psrfits_open()'s field set --- */
    strcpy(hdr->obs_mode, "SEARCH");
    strncpy(hdr->telescope, telescope_name(fb_hdr.telescope_id), 24);
    strcpy(hdr->observer, "unknown");
    strncpy(hdr->source, fb_hdr.source_name[0] ? fb_hdr.source_name : "unknown", 24);
    hdr->frontend[0] = '\0';
    strncpy(hdr->backend, machine_name(fb_hdr.machine_id), 24);
    hdr->project_id[0] = '\0';
    mjd_to_date_obs(fb_hdr.tstart, hdr->date_obs, sizeof(hdr->date_obs));

    /* src_raj/src_dej are SIGPROC's HHMMSS.SSSS / DDMMSS.SSSS encoding */
    {
        double raj = fb_hdr.src_raj, dej = fb_hdr.src_dej;
        int rh, rm, dd, dm;
        double rs, ds;
        int dsign = (dej < 0) ? -1 : 1;
        dej = fabs(dej);
        rh = (int)(raj / 10000.0);
        rm = (int)((raj - rh * 10000.0) / 100.0);
        rs = raj - rh * 10000.0 - rm * 100.0;
        dd = (int)(dej / 10000.0);
        dm = (int)((dej - dd * 10000.0) / 100.0);
        ds = dej - dd * 10000.0 - dm * 100.0;
        snprintf(hdr->ra_str, 16, "%02d:%02d:%07.4f", rh, rm, rs);
        snprintf(hdr->dec_str, 16, "%s%02d:%02d:%06.3f",
                 dsign < 0 ? "-" : "", dd, dm, ds);
    }

    strcpy(hdr->poln_type, "LIN");
    strcpy(hdr->poln_order, fb_hdr.nifs == 1 ? "AA+BB" : "");
    hdr->summed_polns = (fb_hdr.nifs == 1) ? 1 : 0;
    strcpy(hdr->track_mode, "TRACK");
    strcpy(hdr->cal_mode, "OFF");
    hdr->feed_mode[0] = '\0';

    hdr->MJD_epoch = (long double)fb_hdr.tstart;
    hdr->start_day = (int)floorl(hdr->MJD_epoch);
    hdr->start_sec = (double)((hdr->MJD_epoch - hdr->start_day) * 86400.0L);
    hdr->start_lst = 0.0;
    hdr->dt = fb_hdr.tsamp;

    /* Keep foff's sign: sub.dat_freqs[] will be filled fch1 + i*foff, same
     * raw physical channel order as the file -- see the comment in
     * filterbank_read_subint() below for why this matters. */
    hdr->df = fb_hdr.foff;
    hdr->orig_df = fb_hdr.foff;
    hdr->BW = fb_hdr.nchans * fb_hdr.foff;
    hdr->fctr = fb_hdr.fch1 + (fb_hdr.nchans - 1) / 2.0 * fb_hdr.foff;
    hdr->orig_nchan = fb_hdr.nchans;
    hdr->nchan = fb_hdr.nchans;
    hdr->chan_dm = 0.0;

    hdr->ra2000 = 0.0; hdr->dec2000 = 0.0;
    hdr->azimuth = fb_hdr.az_start; hdr->zenith_ang = fb_hdr.za_start;
    hdr->beam_FWHM = 0.0;
    hdr->cal_freq = 0.0; hdr->cal_dcyc = 0.0; hdr->cal_phs = 0.0;
    hdr->feed_angle = 0.0; hdr->scanlen = 0.0;
    hdr->fd_sang = 0.0; hdr->fd_xyph = 0.0;
    hdr->fd_hand = 1; hdr->be_phase = 1;
    hdr->scan_number = 1;

    hdr->nbits = fb_hdr.nbits;
    hdr->orig_nbits = fb_hdr.nbits;
    hdr->nbin = 0;
    hdr->npol = fb_hdr.nifs;
    hdr->rcvr_polns = fb_hdr.nifs;
    hdr->offset_subint = 0;
    hdr->onlyI = 0;
    hdr->ds_time_fact = 1;
    hdr->ds_freq_fact = 1;

    /* Choose an internal "subint" block size. This is a free parameter for
     * filterbank input (unlike PSRFITS, where NSBLK is baked into the raw
     * file); 2048 spectra/block is a reasonable, arbitrary default. */
    hdr->nsblk = 2048;

    fb_bytes_per_sample_row = ((long long)hdr->nchan * hdr->npol * hdr->nbits) / 8;
    /* Mirrors psrfits_open()'s SEARCH_MODE computation of sub->bytes_per_subint;
     * init_subbanding() allocates pfi->sub.rawdata with exactly this size, so
     * leaving it unset (garbage stack value) causes a heap buffer overflow
     * the moment filterbank_read_subint() fread()s a full row into it. */
    pf->sub.bytes_per_subint = (int)(fb_bytes_per_sample_row * hdr->nsblk);

    pf->rownum = 1;
    pf->tot_rows = 0;
    pf->rows_per_file = (int)(fb_total_samples / hdr->nsblk);
    pf->N = 0;
    pf->T = 0.0;
    pf->status = 0;
    fb_global_pos = 0;
    fb_open_idx = -1;

    return 0;
}


/* Convert `n` spectra of raw packed filterbank samples (already read into
 * `raw`, one row of n*nchan*npol samples, nbits-packed) into floats,
 * mirroring apply_scales_and_offsets() with scale=1, offset=0 (filterbank
 * has no per-channel scale/offset concept -- the raw value IS the
 * calibrated value), while honoring the header's signed/unsigned flag. */
static void fb_samples_to_float(int nchan, int npol, int n, int nbits,
                                int is_signed, unsigned char *raw, float *out)
{
    long long ii;
    long long N = (long long)nchan * npol * n;

    if (nbits == 8) {
        if (is_signed) {
            signed char *p = (signed char *)raw;
            for (ii = 0; ii < N; ii++) out[ii] = (float)p[ii];
        } else {
            for (ii = 0; ii < N; ii++) out[ii] = (float)raw[ii];
        }
    } else if (nbits == 16) {
        if (is_signed) {
            short *p = (short *)raw;
            for (ii = 0; ii < N; ii++) out[ii] = (float)p[ii];
        } else {
            unsigned short *p = (unsigned short *)raw;
            for (ii = 0; ii < N; ii++) out[ii] = (float)p[ii];
        }
    } else if (nbits == 32) {
        float *p = (float *)raw;
        for (ii = 0; ii < N; ii++) out[ii] = p[ii];
    } else if (nbits == 4) {
        unsigned char *tmp = (unsigned char *)malloc(N);
        unpack_4bit_to_8bit_unsigned(raw, tmp, N);
        for (ii = 0; ii < N; ii++) out[ii] = (float)tmp[ii];
        free(tmp);
    } else if (nbits == 2) {
        unsigned char *tmp = (unsigned char *)malloc(N);
        unpack_2bit_to_8bit_unsigned(raw, tmp, N);
        for (ii = 0; ii < N; ii++) out[ii] = (float)tmp[ii];
        free(tmp);
    } else {
        fprintf(stderr, "Error: unsupported filterbank nbits=%d\n", nbits);
        exit(1);
    }
}


int filterbank_read_subint(struct psrfits *pf)
{
    struct hdrinfo *hdr = &(pf->hdr);
    struct subint *sub = &(pf->sub);
    int ii;
    long long nread;

    if (fb_global_pos + hdr->nsblk > fb_total_samples) {
        printf("Finished with filterbank input file(s).\n");
        pf->status = 1;
        return pf->status;
    }

    for (ii = 0; ii < hdr->nchan; ii++)
        sub->dat_freqs[ii] = (float)(fb_hdr.fch1 + ii * fb_hdr.foff);
    for (ii = 0; ii < hdr->nchan; ii++)
        sub->dat_weights[ii] = 1.0;
    for (ii = 0; ii < hdr->nchan * hdr->npol; ii++) {
        sub->dat_offsets[ii] = 0.0;
        sub->dat_scales[ii] = 1.0;
    }

    /* May transparently cross into the next file if hdr->nsblk doesn't
     * evenly divide the current file's remaining samples -- see
     * fb_read_global()'s comment. */
    nread = fb_read_global(fb_global_pos, hdr->nsblk, sub->rawdata);
    if (nread != hdr->nsblk) {
        printf("Finished with filterbank input file(s) (short read).\n");
        pf->status = 1;
        return pf->status;
    }

    sub->tsubint = hdr->nsblk * hdr->dt;
    sub->offs = (pf->rownum - 1 + 0.5) * sub->tsubint;
    sub->lst = 0.0; sub->ra = 0.0; sub->dec = 0.0;
    sub->glon = 0.0; sub->glat = 0.0;
    sub->feed_ang = 0.0; sub->pos_ang = 0.0; sub->par_ang = 0.0;
    sub->tel_az = 0.0; sub->tel_zen = 0.0;

    fb_global_pos += hdr->nsblk;
    pf->rownum++;
    pf->tot_rows++;
    pf->status = 0;
    return 0;
}


int filterbank_read_part_DATA(struct psrfits *pf, int N, int numunsigned,
                              float *fbuffer)
{
    struct hdrinfo *hdr = &(pf->hdr);
    long long bytes_to_read;
    unsigned char *buffer;
    long long nread;
    (void)numunsigned;  /* filterbank data is already single-Stokes; no
                          * poln-dependent (un)signedness split needed */

    bytes_to_read = ((long long)hdr->nchan * hdr->npol * N * hdr->nbits) / 8;
    buffer = (unsigned char *)malloc(bytes_to_read);

    if (fb_global_pos + N > fb_total_samples) {
        free(buffer);
        pf->status = 1;
        return pf->status;
    }

    /* Peek at the next N spectra *without* advancing the committed file
     * position -- filterbank_read_subint() (called right after this, in
     * get_current_row()) is what actually commits the advance, exactly
     * matching psrfits_read_part_DATA()/psrfits_read_subint()'s division
     * of labor for PSRFITS input. fb_read_global() takes an explicit
     * position and has no side effects, so "peeking" is simply passing
     * fb_global_pos without updating it afterward -- no save/restore of
     * file-open state needed even if the peek spans a file boundary. */
    nread = fb_read_global(fb_global_pos, N, buffer);
    if (nread != N) {
        free(buffer);
        pf->status = 1;
        return pf->status;
    }

    fb_samples_to_float(hdr->nchan, hdr->npol, N, hdr->nbits,
                        fb_hdr.signedints, buffer, fbuffer);
    free(buffer);
    return 0;
}


int filterbank_close(void)
{
    if (fb_file) {
        fclose(fb_file);
        fb_file = NULL;
    }
    fb_open_idx = -1;
    free(fb_file_headerlen); fb_file_headerlen = NULL;
    free(fb_file_nsamples); fb_file_nsamples = NULL;
    free(fb_file_start); fb_file_start = NULL;
    return 0;
}
