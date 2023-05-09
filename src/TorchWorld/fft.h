//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//
// These functions and variables are defined to use FFT as well as FFTW
// Please see fft.cpp to show the detailed information
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_FFT_H_
#define TORCHWORLD_FFT_H_

#include "TorchWorld/macrodefinitions.h"

// Commands for FFT (This is the same as FFTW)
#define TW_FFT_FORWARD 1
#define TW_FFT_BACKWARD 2
#define TW_FFT_ESTIMATE 3

namespace tw::FFT {
// Complex number for FFT
    typedef double fft_complex[2];
// Struct used for FFT
    typedef struct {
        int n;
        int sign;
        unsigned int flags;
        fft_complex *c_in;
        double *in;
        fft_complex *c_out;
        double *out;
        double *input;
        int *ip;
        double *w;
    } fft_plan;

    fft_plan fft_plan_dft_1d(int n, fft_complex *in, fft_complex *out, int sign,
                             unsigned int flags);

    fft_plan fft_plan_dft_c2r_1d(int n, fft_complex *in, double *out,
                                 unsigned int flags);

    fft_plan fft_plan_dft_r2c_1d(int n, double *in, fft_complex *out,
                                 unsigned int flags);

    void fft_execute(fft_plan p);

    void fft_destroy_plan(fft_plan p);
}

#endif  // TORCHWORLD_FFT_H_
