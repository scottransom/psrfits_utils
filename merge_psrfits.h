#include "psrfits.h"
typedef struct {
    int status;
    int chan_id;
    char filename[128];
    struct psrfits pf;
    unsigned char *data;


} thread_args;

void *merge_psrfits_thread(void *args);

pthread_mutex_t lock_read;
