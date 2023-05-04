//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
//
// C++/Torch version of WORLD
// Port by: Daniel Kamp (daniel [at] expressivelabs.net)
// Last update: 2023/05/03
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_CHEAPTRICK_H_
#define TORCHWORLD_CHEAPTRICK_H_

#include <TorchWorld/macrodefinitions.h>
#include <TorchWorld/common.h>

#include <torch/torch.h>
#include <utility>

namespace tw {
    class CheapTrick {
        struct Options {
            double q1;
            double f0_floor;
            int fft_size;
            int fs;

            explicit Options(int _fs);
            void calcFFTSize();
            [[nodiscard]] double getF0Floor() const;
        };

        public:
            // TODO: Add more Torchifications (i.e. torch::Tensor instead of std::vector)
            CheapTrick(std::shared_ptr<std::vector<double>> _x, int _x_length, std::shared_ptr<std::vector<double>> _temporal_positions,  const torch::Tensor &_f0, int _f0_length, const torch::Tensor &_spectrogram, const Options& _option) :
                x(std::move(_x)),
                x_length(_x_length),
                temporal_positions(std::move(_temporal_positions)),
                f0(_f0),
                f0_length(_f0_length),
                spectrogram(_spectrogram),
                options(_option) {}

            ~CheapTrick() = default;

            void run();
            void initializeFFTs();

            void step();

            void applyWindowing();
            void calculatePowerSpectrum();
            void addNoise();
            void smoothingWithRecovery();

            void generateWindow(int half_window_length, std::vector<double>& window, std::vector<int>& base_index, std::vector<int>& safe_index) const;

        private:
            double current_f0 = 0.0;
            double current_position = 0.0;

            int x_length, f0_length;

            std::shared_ptr<std::vector<double>> x = nullptr;
            std::shared_ptr<std::vector<double>> temporal_positions = nullptr;
            std::vector<double> spectral_envelope;

            const torch::Tensor &f0;
            const torch::Tensor &spectrogram;

            ForwardRealFFT* forward_real_fft = nullptr;
            InverseRealFFT* inverse_real_fft = nullptr;

            const Options& options;
    };
}

#endif  // TORCHWORLD_CHEAPTRICK_H_
