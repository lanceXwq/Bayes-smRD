#include "math.h"
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   double *lh, *icp, *tp, *r;
   double *A, *c;

   lh = mxGetDoubles(prhs[0]);  /* pointer to first input matrix */
   icp = mxGetDoubles(prhs[1]); /* pointer to second input matrix */
   tp = mxGetDoubles(prhs[2]);  /* pointer to third input matrix */
   r = mxGetDoubles(prhs[3]);

   mwSize K, M;
   K = mxGetM(prhs[0]);
   M = mxGetN(prhs[0]);

   plhs[0] = mxCreateDoubleMatrix(K, 1, mxREAL);
   c = mxGetDoubles(plhs[0]);
   plhs[1] = mxCreateDoubleMatrix(K, M, mxREAL);
   A = mxGetDoubles(plhs[1]);

   int mK;
   mwIndex k, m, mm;
   double p[M], ss, s = 0;

   for (m = 0; m < M; m++)
   {
      mK = m * K;
      A[mK] = lh[mK] * icp[m];
      s += A[mK];
   }
   for (m = 0; m < M; m++)
   {
      A[m * K] /= s;
   }

   for (k = 1; k < K; k++)
   {
      s = 0;
      for (m = 0; m < M; m++)
      {
         ss = 0;
         mK = m * K;
         for (mm = 0; mm < M; mm++)
         {
            ss += A[k - 1 + mm * K] * tp[mm + m * M];
         }
         A[k + m * K] = ss * lh[k + mK];
         s += A[k + mK];
      }
      for (m = 0; m < M; m++)
      {
         A[k + m * K] /= s;
      }
   }

   s = A[K - 1];
   m = 0;
   while (r[K - 1] > s)
   {
      m++;
      s += A[K - 1 + m * K];
   }
   c[K - 1] = m + 1;

   for (k = K - 1; k-- > 0;)
   {
      s = 0;
      for (m = 0; m < M; m++)
      {
         p[m] = A[k + m * K] * tp[m + (int)(c[k + 1] - 1) * M];
         s += p[m];
      }
      m = 0;
      r[k] *= s;
      ss = p[0];
      while (r[k] > ss)
      {
         m++;
         ss += p[m];
      }
      c[k] = m + 1;
   }
}