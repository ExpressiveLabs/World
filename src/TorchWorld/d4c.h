//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
//
// C++/Torch version of WORLD
// Port by: Daniel Kamp (daniel [at] expressivelabs.net)
// Last update: 2023/05/03
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_D4C_H_
#define TORCHWORLD_D4C_H_

#include <TorchWorld/macrodefinitions.h>
#include <TorchWorld/constantnumbers.h>
#include <TorchWorld/common.h>

#include <torch/torch.h>

namespace tw {
    class D4C {
        public:
            struct Options {
                double threshold;
                int fs;
                int fft_size;

                Options(int _fs, int _fft_size) : threshold(world::kThreshold), fs(_fs), fft_size(_fft_size) {}
            };

            D4C(std::shared_ptr<std::vector<double>> _x, int _x_length, std::shared_ptr<std::vector<double>> _temporal_positions, torch::Tensor& _f0, int _f0_length, torch::Tensor& _aperiodicity, const Options &options) :
                x(std::move(_x)),
                x_length(_x_length),
                temporal_positions(std::move(_temporal_positions)),
                f0(_f0),
                f0_length(_f0_length),
                aperiodicity(_aperiodicity),
                options(options) {}

            ~D4C() = default;

            void run();
            void step();

            void applyWindow(int window_type, double window_length_ratio);
            void staticCentroid();
            void getCentroid(std::vector<double>& centroid);
            void smoothedPowerSpectrum();
            void staticGroupDelay();
            void coarseAperiodicity();
            void getAperiodicity(const torch::Tensor& ap);

            void LoveTrain(std::vector<double>& aperiodicity0);
            double LoveTrainSub(int b0, int b1, int b2);

            void initializeFFTs();
            void initializeTensor();

            void generateWindow(int half_window_length, int window_type, double window_length_ratio, std::vector<int>& base_index, std::vector<int>& safe_index);

        private:
            double current_position = 0.0, current_f0 = 0.0;
            int x_length, f0_length, window_length = 0, num_aperiodicities = 0;

            std::shared_ptr<std::vector<double>> x;
            std::shared_ptr<std::vector<double>> temporal_positions;

            std::vector<double> coarse_aperiodicity;
            std::vector<double> coarse_frequency_axis;
            std::vector<double> frequency_axis;
            std::vector<double> static_centroid;
            std::vector<double> smoothed_power_spectrum;
            std::vector<double> static_group_delay;

            std::vector<double> window;

            torch::Tensor& f0;
            torch::Tensor& aperiodicity;

            Common::ForwardRealFFT* forward_real_fft = nullptr;

            const Options& options;
    };
}

#endif  // TORCHWORLD_D4C_H_
