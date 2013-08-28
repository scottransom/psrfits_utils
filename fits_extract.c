#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "fitsio.h"

#define NFIELDS_I      6

void usage() {
    printf("Usage:  fits_extract inputfile outputfile [options]\n"
    	   "Options:\n"
	   "   -h, --help			Print this\n"
	   "   -s, --start			Start of the ouput file\n"
	   "   -e, --end 	 		End of the output file\n"
	   );

}    


int main(int argc, char *argv[])
{

    static struct option long_opts[] = {
        {"start",  1, NULL, 's'},
        {"end",    1, NULL, 'e'},
	{0,0,0,0}
    };	

    fitsfile *infptr, *outfptr;   /* FITS file pointers defined in fitsio.h */
    int opt, opti;
    int i, status = 0, ii = 1, jj=1;       /* status must always be initialized = 0  */
    int lo_row, hi_row;
    unsigned char *buffer = 0;

    long nrows;
    int ncols;
    int tfields;
    long  pcount;
    char  extname[FLEN_VALUE];
    char  *ttype[NFIELDS_I];
    char  *tform[NFIELDS_I];
    char  *tunit[NFIELDS_I];
    char  *tcomm[NFIELDS_I];

    float start=0.0, end=1.0;
    while ((opt=getopt_long(argc,argv,"s:e:h",long_opts,&opti))!=-1) {
        switch (opt) {
	    case 's':
	        start = atof(optarg);
		break;
	    case 'e':
	        end = atof(optarg);
		break;
	    case 'h':
	    default:
	        usage();
		exit(0);
		break;
	}	

    }

    if (argc < 3)
    {
      usage();
      return(0);
    }

    //printf("> %s %s\n", argv[1], argv[2]);
    //printf("> %s %s\n", argv[optind], argv[2]);


   /* Allocate space for the table parameters and initialize */
   for (i=0; i<NFIELDS_I; i++)
   {
      ttype[i] = (char *) malloc(FLEN_VALUE*sizeof(char));
      tform[i] = (char *) malloc(FLEN_VALUE*sizeof(char));
      tunit[i] = (char *) malloc(FLEN_VALUE*sizeof(char));
      tcomm[i] = (char *) malloc(FLEN_CARD*sizeof(char));
      strcpy(ttype[i], " ");
   }


    /* Open the input file */
    if ( !fits_open_file(&infptr, argv[optind], READONLY, &status) )
    {
      /* Create the output file */
      if ( !fits_create_file(&outfptr, argv[optind+1], &status) )
      {
        /* Copy every HDU until we get an error */
        //while( !fits_movabs_hdu(infptr, ii++, NULL, &status) )
	
	/* Copy first HDU */
        fits_copy_hdu(infptr, outfptr, 0, &status);

	fits_movnam_hdu(infptr, BINARY_TBL, "SUBINT", 0, &status);

	fits_get_num_rows(infptr, &nrows, &status);
	printf("Number of rows = %ld\n", nrows);
	
	fits_get_num_cols(infptr, &ncols, &status);
	printf("Number of cols = %d\n", ncols);
		
	long width;
	fits_read_key(infptr, TLONG, "NAXIS1", &width, NULL, &status);
	printf("Width = %ld\n", width);
	buffer = (unsigned char *) malloc(width);
	
	/* Copy the fits header */
	fits_copy_header(infptr, outfptr, &status);

//#if 0
	/* Determine low and high row number ot copy */
	lo_row = 1 + start * nrows;
	hi_row = end * nrows;
	printf("lo_row=%d hi_row=%d\n", lo_row, hi_row);

	for (ii=lo_row, jj=1; ii<=hi_row; ii++, jj++) {
	    printf("\rRead %d/%ld", ii, nrows);
	    fits_read_tblbytes( infptr,  ii, 1, width, buffer, &status);

	    fits_write_tblbytes(outfptr, jj, 1, width, buffer, &status);
	}
	nrows = (hi_row-lo_row)+1;
	fits_update_key(outfptr, TLONG, "NAXIS2", &nrows, 0, &status);
	printf("\nDone\n");
//#endif
        /* Reset status after normal error */
        //if (status == END_OF_FILE) status = 0;

        fits_close_file(outfptr,  &status);
      }
      fits_close_file(infptr, &status);
    }

    /* if error occured, print out error message */
    if (status) fits_report_error(stderr, status);
    return(status);
}
