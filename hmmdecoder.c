#include "hmmdecoder.h"

void hmm_decoder(seq_t sequence[LEN],mat_t tran_mat[NOS][NOS],mat_t  emi_mat[NOS][NOE],mat_t  output[NOS][LEN])
{

mat_t forward_seq[NOS][LEN+1];
mat_t backward_seq[NOS][LEN+1];
mat_t scaler[LEN+1];
/*mat_t temp_out[NOS][LEN];*/
int i;
int j;
int k;
int index;
float temp;

forward_seq[0][0]=1;
forward_seq[1][0]=0;
forward_seq[2][0]=0;
backward_seq[0][LEN]=1;
backward_seq[1][LEN]=1;
backward_seq[2][LEN]=1;
scaler[0]=1;

Forward_outer_Loop: for(i=1;i<=LEN;i++)
{
  index=sequence[i-1];
  index-=1;
  
 Forward_Inner_LOOP: for(j=0;j<NOS;j++)
                    {
                       temp=0;
                       Forward_MAC_LOOP: for(k=0;k<NOS;k++)
                                       {temp+=forward_seq[k][i-1]*tran_mat[k][j];
                                       }
                     forward_seq[j][i]=temp*emi_mat[j][index];
                     scaler[i]+=forward_seq[j][i];
                     }
   
Forward_Scaling_Loop: for(j=0;j<NOS;j++)
                       {
                         forward_seq[j][i]/=scaler[i];
                       }
}

Backward_outer_LOOP:for(i=LEN;i>0;i--)
{
             index=sequence[i-1];
             index-=1;
Backward_Inner_LOOP:for(j=0;j<NOS;j++)
             {
               temp=0;
               Backward_MAC_LOOP:for(k=0;k<NOS;k++)
                  temp+=tran_mat[j][k]*backward_seq[k][i]*emi_mat[k][index];
               backward_seq[j][i-1]=temp/scaler[i];
              }
}

Output_outer_LOOP:for(i=1;i<=LEN;i++)
{
	Output_inner_LOOP:for(k=0;k<NOS;k++)
         output[k][i-1]=forward_seq[k][i]*backward_seq[k][i];
}

}
