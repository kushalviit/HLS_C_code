#ifndef HMMDECODER_H_
#define HMMDECODER_H_
#define LEN 128
#define NOS 3
#define NOE 8

typedef int seq_t;
typedef float mat_t;

void hmm_decoder(seq_t sequence[LEN],
                mat_t tran_mat[NOS][NOS],
                mat_t emi_mat[NOS][NOE],mat_t output[NOS][LEN]);

#endif

