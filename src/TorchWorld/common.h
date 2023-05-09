//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_COMMON_H_
#define TORCHWORLD_COMMON_H_

#include <vector>
#include <torch/torch.h>

#include "TorchWorld/fft.h"
#include "TorchWorld/macrodefinitions.h"

namespace tw::Common {

//-----------------------------------------------------------------------------
// Structs on FFT
//-----------------------------------------------------------------------------
// Forward FFT in the real sequence
    struct ForwardRealFFT {
        int fft_size;
        double *waveform;
        FFT::fft_complex *spectrum;
        FFT::fft_plan forward_fft;
    } ;

// Inverse FFT in the real sequence
    struct InverseRealFFT {
        int fft_size;
        double *waveform;
        FFT::fft_complex *spectrum;
        FFT::fft_plan inverse_fft;
    } ;

// Inverse FFT in the complex sequence
    struct InverseComplexFFT {
        int fft_size;
        FFT::fft_complex *input;
        FFT::fft_complex *output;
        FFT::fft_plan inverse_fft;
    };

// Minimum phase analysis from logarithmic power spectrum
    struct MinimumPhaseAnalysis {
        int fft_size;
        double *log_spectrum;
        FFT::fft_complex *minimum_phase_spectrum;
        FFT::fft_complex *cepstrum;
        FFT::fft_plan inverse_fft;
        FFT::fft_plan forward_fft;
    };

//-----------------------------------------------------------------------------
// GetSuitableFFTSize() calculates the suitable FFT size.
// The size is defined as the minimum length whose length is longer than
// the input sample.
//
// Input:
//   sample : Length of the input signal
//
// Output:
//   Suitable FFT size
//-----------------------------------------------------------------------------
    int GetSuitableFFTSize(int sample);

//-----------------------------------------------------------------------------
// These four functions are simple max() and min() function
// for "int" and "double" type.
//-----------------------------------------------------------------------------
    inline int MyMaxInt(int x, int y) {
        return x > y ? x : y;
    }

    inline double MyMaxDouble(double x, double y) {
        return x > y ? x : y;
    }

    inline int MyMinInt(int x, int y) {
        return x < y ? x : y;
    }

    inline double MyMinDouble(double x, double y) {
        return x < y ? x : y;
    }

//-----------------------------------------------------------------------------
// These functions are used in at least two different .cpp files

//-----------------------------------------------------------------------------
// DCCorrection interpolates the power under f0 Hz
// and is used in CheapTrick() and D4C().
//-----------------------------------------------------------------------------
    void DCCorrection(const double *input, double current_f0, int fs, int fft_size,
                      double *output);

//-----------------------------------------------------------------------------
// LinearSmoothing() carries out the spectral smoothing by rectangular window
// whose length is width Hz and is used in CheapTrick() and D4C().
//-----------------------------------------------------------------------------
    void LinearSmoothing(const double *input, double width, int fs, int fft_size,
                         double *output);

//-----------------------------------------------------------------------------
// NuttallWindow() calculates the coefficients of Nuttall window whose length
// is y_length and is used in Dio(), Harvest() and D4C().
//-----------------------------------------------------------------------------
    void NuttallWindow(int y_length, double *y);

//-----------------------------------------------------------------------------
// GetSafeAperiodicity() limit the range of aperiodicity from 0.001 to
// 0.999999999999 (1 - TorchWorld::kMySafeGuardMinimum).
//-----------------------------------------------------------------------------
    inline double GetSafeAperiodicity(double x) {
        return MyMaxDouble(0.001, MyMinDouble(0.999999999999, x));
    }

//-----------------------------------------------------------------------------
// These functions are used to speed up the processing.
// Forward FFT
    void InitializeForwardRealFFT(int fft_size, ForwardRealFFT *forward_real_fft);

    void DestroyForwardRealFFT(ForwardRealFFT *forward_real_fft);

// Inverse FFT
    void InitializeInverseRealFFT(int fft_size, InverseRealFFT *inverse_real_fft);

    void DestroyInverseRealFFT(InverseRealFFT *inverse_real_fft);

// Inverse FFT (Complex)
    void InitializeInverseComplexFFT(int fft_size,
                                     InverseComplexFFT *inverse_complex_fft);

    void DestroyInverseComplexFFT(InverseComplexFFT *inverse_complex_fft);

// Minimum phase analysis (This analysis uses FFT)
    void InitializeMinimumPhaseAnalysis(int fft_size,
                                        MinimumPhaseAnalysis *minimum_phase);

    void GetMinimumPhaseSpectrum(const MinimumPhaseAnalysis *minimum_phase);

    void DestroyMinimumPhaseAnalysis(MinimumPhaseAnalysis *minimum_phase);

}

#endif  // WORLD_COMMON_H_
