// This code is to partially de-disperse and subband
// PSRFITS search-mode data.  Currently it is specifically
// for GUPPI data, however, I intend to make it more general
// eventually.   S. Ransom  Oct 2008
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include "psrfits.h"
#include "psrfits_subband_cmd.h"

// This tests to see if 2 times are within 100ns of each other
#define TEST_CLOSE(a, b) (fabs((a)-(b)) <= 1e-7 ? 1 : 0)

extern double delay_from_dm(double dm, double freq_emitted);
extern int split_root_suffix(char *input, char **root, char **suffix);
extern void avg_std(char *x, int n, double *mean, double *std, int stride);

struct subband_info {
    int nsub;
    int nchan;
    int chan_per_sub;
    int npol;
    int max_early;
    int max_late;
    int max_overlap;
    int buflen;  // Number of samples (in time) in a block
    int bufwid;  // Number of channels times number of polns
    double dm;
    double sub_df;
    double *sub_delays;
    double *chan_delays;
    int *idelays;
    float *sub_freqs;
    float *weights;
    float *offsets;
    float *scales;
    float *chan_avgs;
    float *chan_stds;
    unsigned char *buffer;
    unsigned char *outbuffer;
};


static void print_percent_complete(int current, int number, int reset)
{
   static int newper = 0, oldper = -1;

   if (reset) {
      oldper = -1;
      newper = 0;
   } else {
      newper = (int) (current / (float) (number) * 100.0);
      if (newper < 0)
         newper = 0;
      if (newper > 100)
         newper = 100;
      if (newper > oldper) {
         printf("\r%3d%%", newper);
         fflush(stdout);
         oldper = newper;
      }
   }
}

void get_chan_stats(struct psrfits *pfi, struct subband_info *si){   
    int ii, offset=0;
    double avg, std;

    for (ii = 0 ; ii < si->bufwid ; ii++, offset++) {
        // Only use 1/8 of the total length in order to speed things up
        avg_std((char *)pfi->sub.data+offset, si->buflen/8, &avg, &std, si->bufwid);
        si->chan_avgs[ii] = avg;
        si->chan_stds[ii] = std;
    }
}


void get_sub_stats(struct psrfits *pfo, struct subband_info *si) {
    int ii, stride=si->bufwid/si->chan_per_sub;
    double avg, std;

    for (ii = 0 ; ii < stride ; ii++) {
        avg_std((char *)pfo->sub.data+ii, si->buflen, &avg, &std, stride);
    }
}


int get_current_row(struct psrfits *pfi, struct subband_info *si) {
    static int firsttime = 1, num_pad_blocks = 0;
    static double last_offs, row_duration;
    double diff_offs, dnum_blocks;
    int ii, jj;
    
    if (firsttime) {
        row_duration = pfi->sub.tsubint;
        last_offs = pfi->sub.offs-row_duration;
        firsttime = 0;
    }

    print_percent_complete(pfi->rownum-1, pfi->rows_per_file, 
                           pfi->rownum==1 ? 1 : 0);

    if (num_pad_blocks==0) {  // Try to read the PSRFITS file

        // Read the current row of data
        psrfits_read_subint(pfi);
        diff_offs = pfi->sub.offs - last_offs;

        if (!TEST_CLOSE(diff_offs, row_duration) || pfi->status) {
            if (pfi->status) { // End of the files
                num_pad_blocks = 1;
            } else { // Missing row(s)
                dnum_blocks = diff_offs/row_duration - 1.0;
                num_pad_blocks = (int)(dnum_blocks + 1e-12);
#if 1
                printf("At row %d, found %d dropped rows.\n", 
                       pfi->rownum, num_pad_blocks);
#endif        
                pfi->rownum--;   // Will re-read when no more padding
                pfi->tot_rows--; // Only count "real" rows towards tot_rows
                pfi->N -= pfi->hdr.nsblk;  // Will be re-added below for padding
            }
            // Now fill the main part of si->buffer with the chan_avgs so that
            // it acts like a correctly read block (or row)
            for (ii = 0 ; ii < si->bufwid ; ii++) {
                for (jj = 0 ; jj < si->buflen ; jj++) {
                    pfi->sub.data[jj*si->bufwid+ii] = (char)rint(si->chan_avgs[ii]);
                }
            }

        } else { // Return the row from the file
            // Determine channel statistics
            get_chan_stats(pfi, si);
            last_offs = pfi->sub.offs;
            return 0;
        }
    }

    // Return the same padding as before
    last_offs += row_duration;
    pfi->N += pfi->hdr.nsblk;
    pfi->T = pfi->N * pfi->hdr.dt;
    num_pad_blocks--;
    return num_pad_blocks;
}


/* Average adjacent frequency channels, including dispersion, together
 * to make de-dispersed subbands.  Note: this only works properly for
 * 8-bit data currently.
 */
void make_subbands(struct psrfits *pfi, struct subband_info *si) {
    int ii, jj, kk, itmp;
    char *indata = (char *)pfi->sub.data;
    char *outdata = (char *)si->outbuffer;
    const int dsfact = si->chan_per_sub;
    const int in_bufwid = si->bufwid;
    const int out_bufwid = si->bufwid / si->chan_per_sub;
    const int offset = dsfact >> 1;
    const int shift = ffs(dsfact) - 1; // dsfact's power of 2

    // Iterate over the times 
    for (ii = 0 ; ii < si->buflen ; ii++) {
        int *idelays = si->idelays;
        // Iterate over the output chans/pols 
        for (jj = 0 ; jj < out_bufwid ; jj++) {
            itmp = offset;  // starting value for "rounding"
            // Iterate over the input chans/pols 
            for (kk = 0 ; kk < dsfact ; kk++, idelays++, indata++) {
                itmp += *(indata + *idelays * in_bufwid);
            }
            // The following adds 1/2 of dsfact to the total (which allows
            // for rounding-type behavior) and then divides by dsfact (via
            // a bit shift of the appropriate number of bits, since dsfact
            // should be a power-of-two).
            *outdata++ = itmp >> shift;
        }
    }
}


void init_subbanding(int nsub, double dm, 
                     struct psrfits *pfi, 
                     struct subband_info *si) {
    int ii, jj, kk, cindex;
    double lofreq, dtmp;
    
    si->nsub = nsub;
    si->nchan = pfi->hdr.nchan;
    si->npol = pfi->hdr.npol;
    si->chan_per_sub = si->nchan / si->nsub;
    si->bufwid = si->nchan * si->npol; // Freq * polns
    si->buflen = pfi->hdr.nsblk;  // Time
    if (si->nchan % si->nsub) {
        fprintf(stderr, 
                "Error!  %d channels is not evenly divisible by %d subbands!\n", 
                si->nchan, si->nsub);
        exit(1);
    }
    si->dm = dm;
    si->sub_df = pfi->hdr.orig_df * si->chan_per_sub;
    si->sub_freqs = (float *)malloc(sizeof(float) * si->nsub);
    si->chan_delays = (double *)malloc(sizeof(double) * si->nchan);
    si->sub_delays = (double *)malloc(sizeof(double) * si->nsub);
    si->idelays = (int *)malloc(sizeof(int) * si->nchan * si->npol);
    si->weights = (float *)malloc(sizeof(float) * si->nsub);
    si->offsets = (float *)malloc(sizeof(float) * si->nsub * si->npol);
    si->scales = (float *)malloc(sizeof(float) * si->nsub * si->npol);
    si->chan_avgs = (float *)malloc(sizeof(float) * si->bufwid);
    si->chan_stds = (float *)malloc(sizeof(float) * si->bufwid);

    /* Alloc data buffers for the input PSRFITS file */
    pfi->sub.dat_freqs = (float *)malloc(sizeof(float) * pfi->hdr.nchan);
    pfi->sub.dat_weights = (float *)malloc(sizeof(float) * pfi->hdr.nchan);
    pfi->sub.dat_offsets = (float *)malloc(sizeof(float)
                                          * pfi->hdr.nchan * pfi->hdr.npol);
    pfi->sub.dat_scales  = (float *)malloc(sizeof(float)
                                          * pfi->hdr.nchan * pfi->hdr.npol);
    // This is temporary...
    pfi->sub.data = (unsigned char *)malloc(pfi->sub.bytes_per_subint);
        
    // Read the first row of data
    psrfits_read_subint(pfi);

    // Reset the read counters since we'll re-read
    pfi->rownum--;
    pfi->tot_rows--;
    pfi->N -= pfi->hdr.nsblk;

    // Compute the subband properties, DM delays and offsets
    lofreq = pfi->sub.dat_freqs[0] - pfi->hdr.orig_df * 0.5;
    for (ii = 0, cindex  = 0 ; ii < si->nsub ; ii++) {
        dtmp = lofreq + ((double)ii + 0.5) * si->sub_df;
        si->sub_freqs[ii] = dtmp;
        si->sub_delays[ii] = delay_from_dm(si->dm, dtmp);
// TODO:  set these (weights, offsets, scales) based on input data
// TODO:  also, adjust for the U and V 0.5 offsets?
        si->weights[ii] = 1.0;
        for (jj = 0 ; jj < si->npol ; jj++) {
            si->offsets[jj*si->nsub+ii] = 0.0;
            si->scales[jj*si->nsub+ii] = 1.0;
        }
        // Determine the dispersion delays and convert them
        // to offsets in units of sample times
        for (jj = 0 ; jj < si->chan_per_sub ; jj++, cindex++) {
            si->chan_delays[cindex] = delay_from_dm(si->dm, 
                                                    pfi->sub.dat_freqs[cindex]);
            si->chan_delays[cindex] -= si->sub_delays[ii];
            si->idelays[cindex] = (int)rint(si->chan_delays[cindex] / pfi->hdr.dt);
            // Copy the delays if we have more than 1 poln
            for (kk = 1 ; kk < si->npol ; kk++) {
                si->idelays[cindex+kk*si->nchan] = si->idelays[cindex];
            }
        }
    }

    // Now determine the earliest and latest delays
    si->max_early = si->max_late = 0;
    for (ii = 0 ; ii < si->nchan ; ii++) {
        if (si->idelays[ii] < si->max_early)
            si->max_early = si->idelays[ii];
        if (si->idelays[ii] > si->max_late)
            si->max_late = si->idelays[ii];
    }
    si->max_overlap = abs(si->max_early) + si->max_late;

    // This buffer will hold the input data, plus the bits of data from
    // the previous and next blocks
    si->buffer = (unsigned char *)calloc((si->buflen + 2 * si->max_overlap) * 
                                         si->bufwid, sizeof(unsigned char));
    // The input data will be stored directly in the buffer space
    // So the following is really just an offset into the bigger buffer
    free(pfi->sub.data);  // Free the temporary allocation from above
    pfi->sub.data = si->buffer + si->max_overlap * si->bufwid * sizeof(unsigned char);
    // We need the following since we do out-of-place subbanding
    si->outbuffer = (unsigned char *)calloc(si->nsub * si->npol * si->buflen, 
                                            sizeof(unsigned char));
    
    // re-read the first row (i.e. for "real" this time)
    get_current_row(pfi, si);
    
    // Now fill the first part of si->buffer with the chan_avgs so that
    // it acts like a previously read block (or row)
    for (ii = 0 ; ii < si->bufwid ; ii++) {
        for (jj = 0 ; jj < si->max_overlap ; jj++) {
            si->buffer[jj*si->bufwid+ii] = (char)rint(si->chan_avgs[ii]);
        }
    }
}


void set_output_vals(struct psrfits *pfi, 
                     struct psrfits *pfo, 
                     struct subband_info *si) {
    // Copy everything
    *pfo = *pfi;
    pfo->filenum = 0; // This causes the output files to be created
    pfo->filename[0] = '\0';
    pfo->rownum = 1;
    pfo->tot_rows = 0;
    pfo->N = 0;
    sprintf(pfo->basefilename, "%s_subs", pfi->basefilename);
    // Reset different params
    pfo->sub.dat_freqs = si->sub_freqs;
    pfo->sub.dat_weights = si->weights;
    pfo->sub.dat_offsets = si->offsets;
    pfo->sub.dat_scales  = si->scales;
    pfo->hdr.ds_freq_fact = si->chan_per_sub;
    pfo->hdr.chan_dm = si->dm;
    //pfo->hdr.nchan = si->nsub;
    //pfo->hdr.npol = si->npol;
    //pfo->sub.bytes_per_subint = (pfo->hdr.nbits * pfo->hdr.nchan * 
    //                             pfo->hdr.npol * pfo->hdr.nsblk) / 8;
    pfo->sub.data = si->outbuffer;
}


int main(int argc, char *argv[]) {
    Cmdline *cmd;
    struct psrfits pfi, pfo;
    struct subband_info si;
    int stat=0, padding=0;

    // Call usage() if we have no command line arguments
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(0);
    }
    
    // Parse the command line using the excellent program Clig
    cmd = parseCmdline(argc, argv);

    // Open the input PSRFITs files
    pfi.tot_rows = pfi.N = pfi.T = pfi.status = 0;
    pfi.filenum = 1;
    pfi.filename[0] = '\0';
    sprintf(pfi.basefilename, cmd->argv[0]);
    int rv = psrfits_open(&pfi);
    if (rv) { fits_report_error(stderr, rv); exit(1); }

    // Initialize the subbanding
    // (including reading the first row of data and
    //  putting it in si->buffer)
    init_subbanding(cmd->nsub, cmd->dm, &pfi, &si);

    // Update the output PSRFITS structure
    set_output_vals(&pfi, &pfo, &si);

    // Loop through the data
    do {
        // Put the overlapping parts from the next block into si->buffer
        char *ptr = (char *)(pfi.sub.data + si.buflen * si.bufwid);
        if (padding > 0)
            stat = psrfits_read_part_DATA(&pfi, si.max_overlap, ptr);
        if (stat || padding) { // Need to use padding since we ran out of data
            printf("Adding next-block padding...\n");
            int ii, jj;
            // Now fill the last part of si->buffer with the chan_avgs so that
            // it acts like a correctly read block (or row)
            for (ii = 0 ; ii < si.bufwid ; ii++) {
                for (jj = 0 ; jj < si.max_overlap ; jj++) {
                    ptr[jj*si.bufwid+ii] = (char)rint(si.chan_avgs[ii]);
                }
            }
        }

        // Now create the subbanded row in the output buffer
        make_subbands(&pfi, &si);
        //get_sub_stats(&pfo, &si);

        // Write the new row to the output file
// TODO:  update the other row-based params as well
        pfo.sub.offs = ((pfo.rownum-1)+0.5) * pfo.sub.tsubint;
        psrfits_write_subint(&pfo);

        // Break out of the loop here if stat is set
        if (stat) break;

        // shift the last part of the current row into the "last-row" 
        // part of the data buffer
        memcpy(si.buffer, 
               pfi.sub.data + (si.buflen - si.max_overlap) * si.bufwid,
               si.max_overlap * si.bufwid);

        // Read the next row (or padding)
        padding = get_current_row(&pfi, &si);

    } while (pfi.status == 0);

    psrfits_close(&pfi);
    if (rv) { fits_report_error(stderr, rv); }
    psrfits_close(&pfo);
    if (rv) { fits_report_error(stderr, rv); }
    exit(0);
}
