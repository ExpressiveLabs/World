//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_CODEC_H_
#define TORCHWORLD_CODEC_H_

#include <torch/torch.h>
#include "TorchWorld/macrodefinitions.h"

namespace tw {
    class Codec {
        public:
            static void InitializeAperiodicity(int f0_length, int fft_size, torch::Tensor &aperiodicity);
            static int GetNumberOfAperiodicities(int fs);
            static void CodeAperiodicity(const double *const *aperiodicity, int f0_length, int fs, int fft_size, double **coded_aperiodicity);
            static void DecodeAperiodicity(const torch::Tensor &coded_aperiodicity, int f0_length, int fs, int fft_size, torch::Tensor &aperiodicity);
            static void GetAperiodicity(const double *coarse_frequency_axis, const double *coarse_aperiodicity, int number_of_aperiodicities, const double *frequency_axis, int fft_size, const torch::Tensor &aperiodicity);

            static void CodeSpectralEnvelope(const double *const *spectrogram, int f0_length, int fs, int fft_size, int number_of_dimensions, double **coded_spectral_envelope);
            static void DecodeSpectralEnvelope(const double *const *coded_spectral_envelope, int f0_length, int fs, int fft_size, int number_of_dimensions, double **spectrogram);

            static int CheckVUV(const torch::Tensor &coarse_aperiodicity, int number_of_aperiodicities, double *tmp_aperiodicity);
    };



//-----------------------------------------------------------------------------
// DecodeSpectralEnvelope decodes the coded spectral envelope.
//
// Input:
//   coded_aperiodicity   : Coded aperiodicity
//   f0_length            : Length of F0 contour
//   fs                   : Sampling frequency
//   fft_size             : FFT size
//   number_of_dimensions : Parameter for compression
//
// Output:
//   spectrogram
//-----------------------------------------------------------------------------

}

#endif  // WORLD_CODEC_H_
