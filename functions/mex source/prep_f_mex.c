#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int K;
   double *X, *P, h, *phi_m, *xi, *X0, *inv_var;
   double *f;

   K = mxGetScalar(prhs[0]);
   X = mxGetDoubles(prhs[1]);
   P = mxGetDoubles(prhs[2]);
   h = mxGetScalar(prhs[3]);
   phi_m = mxGetDoubles(prhs[4]);
   xi = mxGetDoubles(prhs[5]);
   X0 = mxGetDoubles(prhs[6]);
   inv_var = mxGetDoubles(prhs[7]);

   plhs[0] = mxCreateDoubleMatrix(K, 3, mxREAL);
   f = mxGetDoubles(plhs[0]);

   int k = 0;
   for (int d = 0; d < 3; d++)
   {
      f[k] = phi_m[k] * X[k] + 2 * P[k] - h * X0[k] * inv_var[k] - xi[k] * X[k + 1];
      k++;
      for (k; k < d * K + K - 1; k++)
      {
         f[k] = phi_m[k] * X[k] + 2 * P[k] - xi[k - 1] * X[k - 1] - xi[k] * X[k + 1];
      }
      f[k] = phi_m[k] * X[k] + 2 * P[k] - xi[k - 1] * X[k - 1];
      k++;
   }
}