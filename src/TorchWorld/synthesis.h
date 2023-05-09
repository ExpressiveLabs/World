//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_SYNTHESIS_H_
#define TORCHWORLD_SYNTHESIS_H_

#include <TorchWorld/macrodefinitions.h>
#include <TorchWorld/common.h>

#include <torch/torch.h>

#include <utility>

namespace tw {
    class Synthesis {
        public:
            struct Options {
                int fs;
                int fft_size;
                double frame_period;

                int y_length;
                int f0_length;

                Options(int _fs, int _fft_size, double _frame_period, int _y_length, int _f0_length) : fs(_fs), fft_size(_fft_size), frame_period(_frame_period), y_length(_y_length), f0_length(_f0_length) {}
            };

            struct TimeBaseCollection {

            };

            Synthesis(const torch::Tensor& _f0, const torch::Tensor& _spectrogram, const torch::Tensor& _aperiodicity, std::shared_ptr<std::vector<double>> _y, Options& _options) : f0(_f0), spectrogram(_spectrogram), aperiodicity(_aperiodicity), y(std::move(_y)), options(_options) {}
            ~Synthesis() = default;

            void initializeFFTs();

            void run();
            void step();

            int getTimeBase();
            void getTemporalParameters();
            int getPulseLocations();

        private:
            const torch::Tensor &f0, &spectrogram, &aperiodicity;
            std::vector<double> impulse_response, interpolated_vuv, interpolated_f0, dc_remover;

            std::vector<double> pulse_locations, pulse_locations_time_shift;
            std::vector<int> pulse_locations_index;
            std::vector<double> time_axis, coarse_time_axis, coarse_f0, coarse_vuv;

            std::shared_ptr<std::vector<double>> y;

            torch::Device device = torch::kCPU;

            Options& options;

            Common::MinimumPhaseAnalysis* minimum_phase = nullptr;
            Common::InverseRealFFT* inverse_real_fft = nullptr;
            Common::ForwardRealFFT* forward_real_fft = nullptr;

            double lowest_f0 = 0.0;
            double current_vuv = 0.0, current_time = 0.0, fractional_time_shift = 0.0;
            int noise_size = 0;
    };
}

#endif  // TORCHWORLD_SYNTHESIS_H_
