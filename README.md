# Bayes-smRD

This software package includes several Matlab scripts and auxiliary functions, which implement the computational algorithms for the framework Bayes-smRD.

## Quick Start for Linux Users

1) To use this software please place all the files in a directory contained in MATLAB's search path and run `run_main.m`.
2) Running `run_main.m` will automatically generate a synthetic dataset and run analysis. It typically takes about a day to get enough samples.
3) We are working on Bayes-smRD's documentation and highly encourage you to reach out to us at: <weiqing1@asu.edu>, or <spresse@asu.edu>
4) This code is developed and tested in Matlab R2022a.
5) This code includes Matlab MEX functions compiled with GCC and Linux, for other platforms, please reach out to us, and we will be happy to provide source files of these functions!

## Other Operating Systems

For other systems, the MEX functions need to be recompiled. The C source codes for these functions are in [this directory](functions/mex%20source/). To compile and get executables you can simply run `mex filename.c` in the Matlab command window. For more information please refer to the [official documentation](https://www.mathworks.com/help/matlab/ref/mex.html).

## Cite Our Work
