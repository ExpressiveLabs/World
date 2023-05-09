//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_STONEMASK_H_
#define TORCHWORLD_STONEMASK_H_

#include <TorchWorld/macrodefinitions.h>
#include <TorchWorld/common.h>
#include <TorchWorld/fft.h>

namespace tw {
    class StoneMask {
//-----------------------------------------------------------------------------
// StoneMask() refines the estimated F0 by Dio()
//
// Input:
//   x                      : Input signal
//   x_length               : Length of the input signal
//   fs                     : Sampling frequency
//   time_axis              : Temporal information
//   f0                     : f0 contour
//   f0_length              : Length of f0
//
// Output:
//   refined_f0             : Refined F0
//-----------------------------------------------------------------------------
        StoneMask(const double *x, int x_length, int fs, const double *temporal_positions, const double *f0, int f0_length, double *refined_f0);
        ~StoneMask() = default;

        static void GetBaseIndex(double current_position, const double *base_time, int base_time_length, int fs, int *index_raw);
        static void GetMainWindow(double current_position, const int *index_raw, int base_time_length, int fs, double window_length_in_time, double *main_window);
        static void GetDiffWindow(const double *main_window, int base_time_length, double *diff_window);
        static void GetSpectra(const double *x, int x_length, int fft_size, const int *index_raw, const double *main_window, const double *diff_window, int base_time_length, const Common::ForwardRealFFT *forward_real_fft, FFT::fft_complex *main_spectrum, FFT::fft_complex *diff_spectrum);
        static double FixF0(const double *power_spectrum, const double *numerator_i, int fft_size, int fs, double initial_f0, int number_of_harmonics);
        static double GetTentativeF0(const double *power_spectrum, const double *numerator_i, int fft_size, int fs, double initial_f0);
        static double GetMeanF0(const double *x, int x_length, int fs, double current_position, double initial_f0, int fft_size, double window_length_in_time, const double *base_time, int base_time_length);
        static double GetRefinedF0(const double *x, int x_length, int fs, double current_potision, double initial_f0);
    };
}

#endif  // WORLD_STONEMASK_H_
