//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
//
// C++/Torch version of WORLD
// Port by: Daniel Kamp (daniel [at] expressivelabs.net)
// Last update: 2023/05/03
//
// Band-aperiodicity estimation on the basis of the idea of D4C.
//-----------------------------------------------------------------------------

#include <TorchWorld/d4c.h>
#include <TorchWorld/common.h>
#include <TorchWorld/constantnumbers.h>
#include <TorchWorld/matlabfunctions.h>

#include <cmath>
#include <algorithm>  // for std::sort()



namespace tw {
//-----------------------------------------------------------------------------
// SetParametersForGetWindowedWaveform()
//-----------------------------------------------------------------------------
    void D4C::generateWindow(int half_window_length, int window_type, double window_length_ratio, std::vector<int>& base_index, std::vector<int>& safe_index) {
        for (int i = -half_window_length; i <= half_window_length; ++i) {
            base_index[i + half_window_length] = i;
        }

        int origin = matlab_round(current_position * options.fs + 0.001);
        for (int i = 0; i <= half_window_length * 2; ++i) {
            safe_index[i] = MyMinInt(x_length - 1, MyMaxInt(0, origin + base_index[i]));
        }

        // Designing of the window function
        double position;
        if (window_type == world::kHanning) {  // Hanning window
            for (int i = 0; i <= half_window_length * 2; ++i) {
                position = (2.0 * base_index[i] / window_length_ratio) / options.fs;
                window[i] = 0.5 * cos(world::kPi * position * current_f0) + 0.5;
            }
        } else {  // Blackman window
            for (int i = 0; i <= half_window_length * 2; ++i) {
                position = (2.0 * base_index[i] / window_length_ratio) / options.fs;
                window[i] = 0.42 + 0.5 * cos(world::kPi * position * current_f0) +
                            0.08 * cos(world::kPi * position * current_f0 * 2);
            }
        }
    }

//-----------------------------------------------------------------------------
// GetWindowedWaveform() windows the waveform by F0-adaptive window
// In the variable window_type, 1: hanning, 2: blackman
//-----------------------------------------------------------------------------
    void D4C::applyWindow(int window_type, double window_length_ratio) {
        int half_window_length = matlab_round(window_length_ratio * options.fs / current_f0 / 2.0);

        std::vector<int> base_index(half_window_length * 2 + 1);
        std::vector<int> safe_index(half_window_length * 2 + 1);
        window = std::vector<double>(half_window_length * 2 + 1);

        generateWindow(half_window_length, window_type, window_length_ratio, base_index, safe_index);

        // F0-adaptive windowing
        for (int i = 0; i <= half_window_length * 2; ++i) {
            forward_real_fft->waveform[i] = x->at(safe_index[i]) * window[i] + randn() * world::kMySafeGuardMinimum;
        }

        double tmp_weight1 = 0;
        double tmp_weight2 = 0;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            tmp_weight1 += forward_real_fft->waveform[i];
            tmp_weight2 += window[i];
        }

        double weighting_coefficient = tmp_weight1 / tmp_weight2;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            forward_real_fft->waveform[i] -= window[i] * weighting_coefficient;
        }
    }

//-----------------------------------------------------------------------------
// GetCentroid() calculates the energy centroid (see the book, time-frequency
// analysis written by L. Cohen).
//-----------------------------------------------------------------------------
    void D4C::getCentroid(std::vector<double>& centroid) {
        for (int i = 0; i < options.fft_size; ++i) {
            forward_real_fft->waveform[i] = 0.0;
        }

        applyWindow(world::kBlackman, 4.0);

        double power = 0.0;
        for (int i = 0; i <= matlab_round(2.0 * options.fs / current_f0) * 2; ++i) {
            power += forward_real_fft->waveform[i] * forward_real_fft->waveform[i];
        }

        for (int i = 0; i <= matlab_round(2.0 * options.fs / current_f0) * 2; ++i) {
            forward_real_fft->waveform[i] /= sqrt(power);
        }

        fft_execute(forward_real_fft->forward_fft);

        std::vector<double> tmp_real(options.fft_size / 2 + 1);
        std::vector<double> tmp_imag(options.fft_size / 2 + 1);

        for (int i = 0; i <= options.fft_size / 2; ++i) {
            tmp_real[i] = forward_real_fft->spectrum[i][0];
            tmp_imag[i] = forward_real_fft->spectrum[i][1];
        }

        for (int i = 0; i < options.fft_size; ++i) {
            forward_real_fft->waveform[i] *= i + 1.0;
        }

        fft_execute(forward_real_fft->forward_fft);

        for (int i = 0; i <= options.fft_size / 2; ++i) {
            centroid[i] = forward_real_fft->spectrum[i][0] * tmp_real[i] + tmp_imag[i] * forward_real_fft->spectrum[i][1];
        }
    }

//-----------------------------------------------------------------------------
// GetStaticCentroid() calculates the temporally static energy centroid.
// Basic idea was proposed by H. Kawahara.
//-----------------------------------------------------------------------------
    void D4C::staticCentroid() {
        std::vector<double> centroid1(options.fft_size / 2 + 1);
        std::vector<double> centroid2(options.fft_size / 2 + 1);

        getCentroid(centroid1);
        getCentroid(centroid2);

        for (int i = 0; i <= options.fft_size / 2; ++i)
            static_centroid[i] = centroid1[i] + centroid2[i];

        // TODO: C++-ify this
        DCCorrection(static_centroid.data(), current_f0, options.fs, options.fft_size, static_centroid.data());
    }

//-----------------------------------------------------------------------------
// GetSmoothedPowerSpectrum() calculates the smoothed power spectrum.
// The parameters used for smoothing are optimized in advance.
//-----------------------------------------------------------------------------
    void D4C::smoothedPowerSpectrum() {
        for (int i = 0; i < options.fft_size; ++i) {
            forward_real_fft->waveform[i] = 0.0;
        }

        applyWindow(world::kHanning, 4.0);

        fft_execute(forward_real_fft->forward_fft);

        for (int i = 0; i <= options.fft_size / 2; ++i) {
            smoothed_power_spectrum[i] = forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] + forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];
        }

        DCCorrection(smoothed_power_spectrum.data(), current_f0, options.fs, options.fft_size,
                     smoothed_power_spectrum.data());
        LinearSmoothing(smoothed_power_spectrum.data(), current_f0, options.fs, options.fft_size,
                        smoothed_power_spectrum.data());
    }

//-----------------------------------------------------------------------------
// GetStaticGroupDelay() calculates the temporally static group delay.
// This is the fundamental parameter in D4C.
//-----------------------------------------------------------------------------
    void D4C::staticGroupDelay() {
        for (int i = 0; i <= options.fft_size / 2; ++i)
            static_group_delay[i] = static_centroid[i] / smoothed_power_spectrum[i];
        LinearSmoothing(static_group_delay.data(), current_f0 / 2.0, options.fs, options.fft_size, static_group_delay.data());

        std::vector<double> smoothed_group_delay(options.fft_size / 2 + 1);
        LinearSmoothing(static_group_delay.data(), current_f0, options.fs, options.fft_size, smoothed_group_delay.data());

        for (int i = 0; i <= options.fft_size / 2; ++i) {
            static_group_delay[i] -= smoothed_group_delay[i];
        }
    }

//-----------------------------------------------------------------------------
// GetCoarseAperiodicity() calculates the aperiodicity in multiples of 3 kHz.
// The upper limit is given based on the sampling frequency.
//-----------------------------------------------------------------------------
    void D4C::coarseAperiodicity() {
        int boundary = matlab_round(options.fft_size * 8.0 / window_length);
        int half_window_length = window_length / 2;

        for (int i = 0; i < options.fft_size; ++i) forward_real_fft->waveform[i] = 0.0;

        std::vector<double> power_spectrum(options.fft_size / 2 + 1);
        int center;

        for (int i = 0; i < num_aperiodicities; ++i) {
            center = static_cast<int>(world::kFrequencyInterval * (i + 1) * options.fft_size / options.fs);

            for (int j = 0; j <= half_window_length * 2; ++j) {
                forward_real_fft->waveform[j] = static_group_delay[center - half_window_length + j] * window[j];
            }

            fft_execute(forward_real_fft->forward_fft);

            for (int j = 0; j <= options.fft_size / 2; ++j) {
                power_spectrum[j] = forward_real_fft->spectrum[j][0] * forward_real_fft->spectrum[j][0] + forward_real_fft->spectrum[j][1] * forward_real_fft->spectrum[j][1];
            }

            std::sort(power_spectrum.begin(), power_spectrum.begin() + options.fft_size / 2 + 1);

            for (int j = 1; j <= options.fft_size / 2; ++j) {
                power_spectrum[j] += power_spectrum[j - 1];
            }

            coarse_aperiodicity[i + 1] = 10 * log10(power_spectrum[options.fft_size / 2 - boundary - 1] / power_spectrum[options.fft_size / 2]);
        }
    }

    double D4C::LoveTrainSub(int boundary0, int boundary1, int boundary2) {
        std::vector<double> power_spectrum(options.fft_size, 0.0);

        window_length = matlab_round(1.5 * options.fs / current_f0) * 2 + 1;

        applyWindow(world::kBlackman, 3.0);

        for (int i = window_length; i < options.fft_size; ++i)
            forward_real_fft->waveform[i] = 0.0;
        fft_execute(forward_real_fft->forward_fft);

        for (int i = boundary0 + 1; i < options.fft_size / 2 + 1; ++i) {
            power_spectrum[i] = forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] + forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];
        }
        for (int i = boundary0; i <= boundary2; ++i) {
            power_spectrum[i] += +power_spectrum[i - 1];
        }

        return power_spectrum[boundary1] / power_spectrum[boundary2];
    }

//-----------------------------------------------------------------------------
// D4CLoveTrain() determines the aperiodicity with VUV detection.
// If a frame was determined as the unvoiced section, aperiodicity is set to
// very high value as the safeguard.
// If it was voiced section, the aperiodicity of 0 Hz is set to -60 dB.
//-----------------------------------------------------------------------------
    void D4C::LoveTrain(std::vector<double>& aperiodicity0) {
        double lowest_f0 = 40.0;
        int fft_size = static_cast<int>(pow(2.0, 1.0 + static_cast<int>(log(3.0 * options.fs / lowest_f0 + 1) / world::kLog2)));

        // Cumulative powers at 100, 4000, 7900 Hz are used for VUV identification.
        int boundary0 = static_cast<int>(ceil(100.0 * fft_size / options.fs));
        int boundary1 = static_cast<int>(ceil(4000.0 * fft_size / options.fs));
        int boundary2 = static_cast<int>(ceil(7900.0 * fft_size / options.fs));

        for (int i = 0; i < f0_length; ++i) {
            auto val = *f0[i].data_ptr<double>();

            if (val == 0.0) {
                aperiodicity0[i] = 0.0;
                continue;
            }

            aperiodicity0[i] = LoveTrainSub(boundary0, boundary1, boundary2);
        }
    }

//-----------------------------------------------------------------------------
// D4CGeneralBody() calculates a spectral envelope at a temporal
// position. This function is only used in D4C().
// Caution:
//   forward_fft is allocated in advance to speed up the processing.
//-----------------------------------------------------------------------------
    void D4C::step() {
        static_centroid = std::vector<double>(options.fft_size / 2 + 1);
        smoothed_power_spectrum = std::vector<double>(options.fft_size / 2 + 1);
        static_group_delay = std::vector<double>(options.fft_size / 2 + 1);

        staticCentroid();
        smoothedPowerSpectrum();
        staticGroupDelay();
        coarseAperiodicity();

        // Revision of the result based on the F0
        for (int i = 0; i < num_aperiodicities; ++i) {
            coarse_aperiodicity[i + 1] = MyMinDouble(0.0, coarse_aperiodicity[i + 1] + (current_f0 - 100) / 50.0);
        }
    }

    void D4C::getAperiodicity(const torch::Tensor& ap) {
        interp1(coarse_frequency_axis.data(), coarse_aperiodicity.data(), num_aperiodicities + 2, frequency_axis.data(), options.fft_size / 2 + 1, ap.to(c10::kDouble).data_ptr<double>());

        for (int i = 0; i <= options.fft_size / 2; ++i) {
            ap[i] = pow(10.0, ap[i] / 20.0);
        }
    }

    void D4C::initializeFFTs() {
        forward_real_fft = {nullptr};
        InitializeForwardRealFFT(options.fft_size, forward_real_fft);
    }

    void D4C::initializeTensor() {
        for (int i = 0; i < f0_length; ++i) {
            for (int j = 0; j < options.fft_size / 2 + 1; ++j) {
                aperiodicity[i][j] = 1.0 - world::kMySafeGuardMinimum;
            }
        }
    }

    void D4C::run() {
        randn_reseed();

        initializeFFTs();
        initializeTensor();

        int fft_size_d4c = static_cast<int>(pow(2.0, 1.0 + static_cast<int>(log(4.0 * options.fs / world::kFloorF0D4C + 1) / world::kLog2)));
        num_aperiodicities = static_cast<int>(MyMinDouble(world::kUpperLimit, options.fs / 2.0 - world::kFrequencyInterval) / world::kFrequencyInterval);

        // Since the window function is common in D4CGeneralBody(), it is designed here to speed up.
        window_length = (int) (world::kFrequencyInterval * fft_size_d4c / options.fs) * 2 + 1;
        window = std::vector<double>(window_length);
        NuttallWindow(window_length, window.data());

        // D4C Love Train (Aperiodicity of 0 Hz is given by the different algorithm)
        std::vector<double> aperiodicity0(f0_length);
        LoveTrain(aperiodicity0);

        coarse_aperiodicity = std::vector<double>(num_aperiodicities + 2);
        coarse_aperiodicity[0] = -60.0;
        coarse_aperiodicity[num_aperiodicities + 1] = -world::kMySafeGuardMinimum;

        coarse_frequency_axis = std::vector<double>(num_aperiodicities + 2);
        for (int i = 0; i <= num_aperiodicities; ++i) {
            coarse_frequency_axis[i] = i * world::kFrequencyInterval;
        }
        coarse_frequency_axis[num_aperiodicities + 1] = options.fs / 2.0;

        frequency_axis = std::vector<double>(options.fft_size / 2 + 1);
        for (int i = 0; i <= options.fft_size / 2; ++i) {
            frequency_axis[i] = static_cast<double>(i) * options.fs / options.fft_size;
        }

        for (int i = 0; i < f0_length; ++i) {
            if (*f0[i].data_ptr<double>() == 0 || aperiodicity0[i] <= options.threshold) continue;

            current_f0 = MyMaxDouble(world::kFloorF0D4C, *f0[i].data_ptr<double>());
            current_position = temporal_positions->at(i);

            // NOTE: if this doesn't work, check what coarse_aperiodicity[1] actually does
            step();

            // Linear interpolation to convert the coarse aperiodicity into its
            // spectral representation.
            getAperiodicity(aperiodicity[i]);
        }
    }
}  // namespace