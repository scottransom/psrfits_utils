#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <pthread.h>
#include "psrfits.h"
#include "merge_psrfits.h"


void *merge_psrfits_thread(void *_args) {

    int status = 0;

    // -- Get arguments --
    thread_args *args = (thread_args *) _args;

    pthread_mutex_lock(&lock_read);
    args->status = psrfits_read_subint(&args->pf);
    pthread_mutex_unlock(&lock_read);
    //printf("offs = %lf\n", args->pf.sub.offs);

    pthread_exit(NULL);

}

