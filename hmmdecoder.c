#include "hmmdecoder.h"
#include <stdio.h>
void hmm_decoder(seq_t sequence[LEN],mat_t tran_mat[NOS][NOS],mat_t  emi_mat[NOS][NOE],mat_t  output[NOS][LEN])
{

mat_t forward_seq[NOS][LEN+1];
mat_t backward_seq[NOS][LEN+1];
mat_t scaler[LEN+1];
mat_t prod[NOS][NOS];
mat_t prod1[NOS][NOS];
mat_t sum[NOS],sum1[NOS];
/*mat_t temp_out[NOS][LEN];*/
int i,j,k;
int index,index1;
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
  
for(j=0;j<NOS;j++)
for(k=0;k<NOS;k++)
prod[k][j]=forward_seq[k][i-1]*tran_mat[k][j];

temp=0;
for(j=0;j<NOS;j++)
{
 sum[j]=prod[0][j]+prod[1][j]+prod[2][j];
}

for(j=0;j<NOS;j++)
{
 sum[j]*=emi_mat[j][index];
}

for(j=0;j<NOS;j++)
temp+=sum[j];

scaler[i]=temp;

for(j=0;j<NOS;j++)
forward_seq[j][i]=sum[j]/temp;

index1=sequence[LEN-i];
index1-=1;

for(j=0;j<NOS;j++)
for(k=0;k<NOS;k++)
prod1[k][j]=tran_mat[j][k]*backward_seq[k][LEN-i+1]*emi_mat[k][index];



for(j=0;j<NOS;j++)
{
 sum1[j]=prod1[0][j]+prod1[1][j]+prod1[2][j];
}

for(j=0;j<NOS;j++)
backward_seq[j][LEN-i]=sum1[j];

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
