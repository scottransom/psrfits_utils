#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <pthread.h>
#include "psrfits.h"
#include "merge_psrfits_cmd.h"
#include "merge_psrfits.h"


void reorder_data(unsigned char* outbuf, unsigned char *inbuf, int nband, int nspec, int npol, int nchan, int nbits) {
    int band, spec, pol, inoff = 0, outoff = 0;
    int spband = nspec * npol * nchan;
    int spspec = npol * nchan;

    for (spec = 0 ; spec < nspec ; spec++) {
	for (pol = 0 ; pol < npol ; pol++) {
	    for (band = 0 ; band < nband ; band++) {
		inoff = (band * spband + pol * nchan + spec * spspec) * nbits/8.;
		memcpy(outbuf + outoff, inbuf + inoff, nchan * nbits/8.);
		outoff += nchan * nbits/8.;
	    }
	}
    }
}

void copy_subint_params(struct psrfits *pfo, struct psrfits *pfi) {


    struct subint *sub_in, *sub_out;
    sub_in = &(pfi->sub);
    sub_out = &(pfo->sub);

    sub_out->tsubint = sub_in->tsubint;
    sub_out->offs = sub_in->offs;
    sub_out->lst = sub_in->lst;
    sub_out->ra = sub_in->ra;
    sub_out->dec = sub_in->dec;
    sub_out->glon = sub_in->glon;
    sub_out->glat = sub_in->glat;
    sub_out->feed_ang = sub_in->feed_ang;
    sub_out->pos_ang = sub_in->pos_ang;
    sub_out->par_ang = sub_in->par_ang;
    sub_out->tel_az = sub_in->tel_az;
    sub_out->tel_zen = sub_in->tel_zen;
/*
    sub_out.dat_freqs = (float *) &pf.sub.dat_freqs[fargs[i].chan_id * nchan];
    sub_out.dat_weights = (float *) &pf.sub.dat_weights[fargs[i].chan_id * nchan];
    sub_out.dat_offsets = (float *) &pf.sub.dat_offsets[fargs[i].chan_id * nchan * npol];
    sub_out.dat_scales  = (float *) &pf.sub.dat_scales[fargs[i].chan_id * nchan * npol];
    sub_out.data = (unsigned char *) &tmpbuf[fargs[i].chan_id * bytes_per_subint];
    */
}

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
         printf("\r%3d%% ", newper);
	 if (newper==100) printf("\n");
         fflush(stdout);
         oldper = newper;
      }
   }
}



int main(int argc, char *argv[]) {

    int i, j, status=0;
    Cmdline *cmd;

    // Call usage() if we have no command line arguments
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(0);
    }
    
    // Parse the command line using the excellent program Clig
    cmd = parseCmdline(argc, argv);

    int numfiles = cmd->argc;
    pthread_t threads[numfiles];  // Thread ids
    thread_args fargs[numfiles];  // Arguments passed to the threads

    // -- Read headers of the PSRFITS files --
    // TODO :  Need to tune basefilename depending on the original name file
    // ie, for Nancay, there is an extra _freq_ term in the fits filename : 'nuppi_55529_0355+54_{scan}_{freq}_0001.fits'
    char *adr, filename[128];
    struct psrfits pf;
    for (i=0; i<numfiles; i++) {
	sprintf(filename, "%s", cmd->argv[i]);
	adr = (char *) strrchr(filename, '_');
	sprintf(adr, "\0");
        sprintf(fargs[i].pf.basefilename, "%s", filename);
	fargs[i].pf.filenum = 1;
        fargs[i].pf.filename[0] = '\0';  
	fargs[i].pf.status = 0;
        status = psrfits_open(&fargs[i].pf);
        if (status) fits_report_error(stderr, status);
    }

    // -- Determine lowest and highest frequency, then set up the middle freq --
    // -- Also determine frequency order -- 
    double flow, fhigh, fctr;
    flow = fhigh = fargs[0].pf.hdr.fctr;
    for (i=1; i<numfiles; i++) {
      if (flow > fargs[i].pf.hdr.fctr) flow = fargs[i].pf.hdr.fctr;
      if (fhigh < fargs[i].pf.hdr.fctr) fhigh = fargs[i].pf.hdr.fctr;
    }
    fctr = flow + (fhigh - flow)/2.;
    if (cmd->verboseP) printf("f_ctr_low = %lf   f_ctr_high = %lf   fctr = %lf\n", flow, fhigh, fctr);

    // -- Identify the channels --
    for (i=0; i<numfiles; i++) {
        fargs[i].chan_id = (int)((fargs[i].pf.hdr.fctr - flow)/fargs[i].pf.hdr.BW);
        if (cmd->verboseP) printf("File #%d  : chan_id = %d\n", i, fargs[i].chan_id);
    }

    // -- Check if we are missing subbands --
    int nbands;
    nbands = (int)((fhigh - flow)/fargs[0].pf.hdr.BW) + 1;
    if (nbands != numfiles) printf("Are we missing one or more subband file ? (%d subbands and %d numfiles)\n", nbands, numfiles);
    

    // -- Check if same source --
    char source[16];
    sprintf(source, "%s", fargs[0].pf.hdr.source);
    if (cmd->verboseP) printf("Source = %s\n", source);
    for (i=1; i<numfiles; i++) {
        if (strcmp(source, fargs[0].pf.hdr.source)!=0) {
	    fprintf(stderr, " File %s has not the same source as file %s (%s and %s)\n", fargs[i].pf.basefilename, fargs[0].pf.basefilename, fargs[i].pf.hdr.source, fargs[0].pf.hdr.source);
	    exit(-1);
	}    
    }

    // -- Check if same number of channels --
    int nchan;
    nchan = fargs[0].pf.hdr.nchan;
    if(cmd->verboseP) printf("Nchan = %d\n", nchan);
    for (i=1; i<numfiles; i++) {
        if (nchan != fargs[i].pf.hdr.nchan) {
	    fprintf(stderr, " File %s has not the same number of channels as file %s (%d and %d)\n", fargs[i].pf.basefilename, fargs[0].pf.basefilename, fargs[i].pf.hdr.nchan, fargs[0].pf.hdr.nchan);
	    exit(-1);
	}    
    }

    // -- Vars --
    int bytes_per_subint, npol;
    bytes_per_subint = fargs[0].pf.sub.bytes_per_subint;
    npol = fargs[0].pf.hdr.npol;

    // -- Copy header for the output file --
    memcpy(&pf, &fargs[0].pf, sizeof(fargs[0].pf));

    // -- Modify some parameters of the merged PSRFITS file from the non merged --
    adr = (char *) strrchr(pf.basefilename, '_');
    sprintf(adr, "\0");
    
    if(cmd->outdir) {
	char tmpbasefilename[200];
        sprintf(tmpbasefilename, "%s/%s", cmd->outdir, pf.basefilename);
	strncpy(pf.basefilename, tmpbasefilename, 200);
    }	


    sprintf(pf.hdr.obs_mode, "SEARCH");
    pf.tot_rows = pf.N = pf.T = pf.status = 0;
    pf.filenum = 0;
    pf.filename[0] = '\0';
    pf.hdr.orig_nchan = nchan * nbands;
    pf.hdr.nchan = nchan * nbands;
    pf.hdr.fctr = fctr;
    pf.hdr.BW = fargs[0].pf.hdr.BW * nbands;
    pf.sub.bytes_per_subint = bytes_per_subint * nbands;
    unsigned long long filelen = (PSRFITS_MAXFILELEN_SEARCH<<30);  
    if (cmd->verboseP) printf("bytes per subint = %d filelen = %lld freq_fact = %d nbits=%d\n", pf.sub.bytes_per_subint, filelen, pf.hdr.ds_freq_fact, pf.hdr.nbits);
    pf.rows_per_file = filelen / pf.sub.bytes_per_subint;

    // -- Create the new merged psrfits --
    status = psrfits_create(&pf);

    // -- For in-memory transpose of data --
    unsigned char *tmpbuf = NULL;
    tmpbuf = (unsigned char *)malloc(pf.sub.bytes_per_subint);

    // -- Alloc data buffers for the output PSRFITS files --
    pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
    pf.sub.dat_scales  = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
    pf.sub.data = (unsigned char *)malloc(pf.sub.bytes_per_subint);

    // -- Given the frequency, set pointers to the correct position --
    for (i=0; i<numfiles; i++) {
        fargs[i].pf.sub.dat_freqs = (float *) &pf.sub.dat_freqs[fargs[i].chan_id * nchan];
        fargs[i].pf.sub.dat_weights = (float *) &pf.sub.dat_weights[fargs[i].chan_id * nchan];
        fargs[i].pf.sub.dat_offsets = (float *) &pf.sub.dat_offsets[fargs[i].chan_id * nchan * npol];
        fargs[i].pf.sub.dat_scales  = (float *) &pf.sub.dat_scales[fargs[i].chan_id * nchan * npol];
        fargs[i].pf.sub.data = (unsigned char *) &tmpbuf[fargs[i].chan_id * bytes_per_subint];
    }

    // -- Loop through the data --
    int statsum = 0;
    while (!statsum) {

        // -- Create reading threads -- 
	for (i=0; i<numfiles; i++) {
	    status = pthread_create(&threads[i], NULL, &merge_psrfits_thread, &fargs[i]);
	}    

	// -- Waiting for threads to finish --
	for (i=0; i<numfiles; i++) {
	    pthread_join(threads[i], NULL);
		if (cmd->verboseP) {
		    for(j=0;j<fargs[i].pf.hdr.nchan; j++) printf("File #%d  Freq[%03d]=%f\n", i, j, fargs[i].pf.sub.dat_freqs[j]);
		}
	}
 	for (i=0; i<numfiles; i++)  statsum += fargs[i].status;
	if (statsum) break;

	copy_subint_params(&pf, &fargs[0].pf);
	reorder_data(pf.sub.data, tmpbuf, nbands, pf.hdr.nsblk, pf.hdr.npol, nchan, pf.hdr.nbits);

	status = psrfits_write_subint(&pf);

	print_percent_complete(fargs[0].pf.tot_rows, fargs[0].pf.rows_per_file, fargs[0].pf.tot_rows == 1 ? 1:0);
    } 

    // -- Closing all files --
    status = psrfits_close(&pf);
    for (i=0; i<numfiles; i++) {
        status = psrfits_close(&fargs[i].pf);
    }	

    exit(0);
}
