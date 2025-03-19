#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int K;
   double *phi, *xi, *f;
   double *X;

   K = mxGetScalar(prhs[0]);    /* pointer to first input matrix */
   phi = mxGetDoubles(prhs[1]); /* pointer to second input matrix */
   xi = mxGetDoubles(prhs[2]);  /* pointer to third input matrix */
   f = mxGetDoubles(prhs[3]);   /* pointer to fourth input matrix */

   plhs[0] = mxCreateDoubleMatrix(K, 3, mxREAL);
   X = mxGetDoubles(plhs[0]);

   int k, s, ks;
   double w, xip[K];
   for (int d = 0; d < 3; d++)
   {
      s = K * d;
      w = phi[s];
      X[s] = f[s] / w;

      for (k = 1; k <= K - 1; k++)
      {
         ks = k + s;
         xip[k - 1] = xi[ks - 1] / w;
         w = phi[ks] - xi[ks - 1] * xip[k - 1];
         X[ks] = (f[ks] - xi[ks - 1] * X[ks - 1]) / w;
      }

      for (k = K - 2; k >= 0; k--)
      {
         ks = k + s;
         X[ks] = X[ks] - xip[k] * X[ks + 1];
      }
   }
}