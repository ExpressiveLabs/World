//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef WORLD_SYNTHESIS_H_
#define WORLD_SYNTHESIS_H_

#include <torch/torch.h>
#include "world/macrodefinitions.h"

WORLD_BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// Synthesis() synthesize the voice based on f0, spectrogram and
// aperiodicity (not excitation signal).
//
// Input:
//   f0                   : f0 contour
//   f0_length            : Length of f0
//   spectrogram          : Spectrogram estimated by CheapTrick
//   fft_size             : FFT size
//   aperiodicity         : Aperiodicity spectrogram based on D4C
//   frame_period         : Temporal period used for the analysis
//   fs                   : Sampling frequency
//   y_length             : Length of the output signal (Memory of y has been
//                          allocated in advance)
// Output:
//   y                    : Calculated speech
//-----------------------------------------------------------------------------
void Synthesis(const torch::Tensor& f0, int f0_length,
    const torch::Tensor& spectrogram, const torch::Tensor& aperiodicity,
    int fft_size, double frame_period, int fs, int y_length, std::vector<double>& y);

WORLD_END_C_DECLS

#endif  // WORLD_SYNTHESIS_H_
