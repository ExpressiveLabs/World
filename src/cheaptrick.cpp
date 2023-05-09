//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//
// Spectral envelope estimation on the basis of the idea of CheapTrick.
//-----------------------------------------------------------------------------
#include <TorchWorld/cheaptrick.h>
#include <cmath>

#include <TorchWorld/common.h>
#include <TorchWorld/constantnumbers.h>
#include <TorchWorld/matlabfunctions.h>

namespace tw {
    // Function implementations for CheapTrick::Options struct
    CheapTrick::Options::Options(int _fs) {
        q1 = -0.15;
        f0_floor = world::kFloorF0;
        fft_size = 0;
        fs = _fs;

        calcFFTSize();
    }

    void CheapTrick::Options::calcFFTSize() {
        fft_size = static_cast<int>(pow(2.0, 1.0 + static_cast<int>(log(3.0 * fs / f0_floor + 1) / world::kLog2)));
    }

    double CheapTrick::Options::getF0Floor() const {
        return 3.0 * fs / (fft_size - 3.0);
    }

    // Function implementations for CheapTrick class
    void CheapTrick::run() {
        int fft_size = options.fft_size;

        MatlabFunctions::randn_reseed();

        double f0_floor = options.getF0Floor();
        spectral_envelope = std::vector<double>(fft_size);

        initializeFFTs();

        for (int i = 0; i < f0_length; ++i) {
            auto val = f0.index({i}).item<double>();

            current_f0 = val <= f0_floor ? world::kDefaultF0 : val;
            current_position = temporal_positions->at(i);

            step();

            for (int j = 0; j <= fft_size / 2; ++j)
                spectrogram.index({i, j}) = spectral_envelope[j];
        }

        DestroyForwardRealFFT(forward_real_fft);
        DestroyInverseRealFFT(inverse_real_fft);
    }

    //-----------------------------------------------------------------------------
    // SmoothingWithRecovery() carries out the spectral smoothing and spectral
    // recovery on the Cepstrum domain.
    //-----------------------------------------------------------------------------
    void CheapTrick::smoothingWithRecovery() {
        std::vector<double> smoothing_lifter(options.fft_size);
        std::vector<double> compensation_lifter(options.fft_size);
        double q1 = options.q1;

        smoothing_lifter[0] = 1.0;
        compensation_lifter[0] = (1.0 - 2.0 * q1) + 2.0 * q1;
        double quefrency;
        for (int i = 1; i <= forward_real_fft->fft_size / 2; ++i) {
            quefrency = static_cast<double>(i) / options.fs;
            smoothing_lifter[i] = sin(world::kPi * current_f0 * quefrency) /
                                  (world::kPi * current_f0 * quefrency);
            compensation_lifter[i] = (1.0 - 2.0 * q1) + 2.0 * q1 * cos(2.0 * world::kPi * quefrency * current_f0);
        }

        for (int i = 0; i <= options.fft_size / 2; ++i)
            forward_real_fft->waveform[i] = log(forward_real_fft->waveform[i]);
        for (int i = 1; i < options.fft_size / 2; ++i)
            forward_real_fft->waveform[options.fft_size - i] = forward_real_fft->waveform[i];
        fft_execute(forward_real_fft->forward_fft);

        for (int i = 0; i <= options.fft_size / 2; ++i) {
            inverse_real_fft->spectrum[i][0] = forward_real_fft->spectrum[i][0] *
                                               smoothing_lifter[i] * compensation_lifter[i] / options.fft_size;
            inverse_real_fft->spectrum[i][1] = 0.0;
        }
        fft_execute(inverse_real_fft->inverse_fft);

        for (int i = 0; i <= options.fft_size / 2; ++i)
            spectral_envelope[i] = exp(inverse_real_fft->waveform[i]);
    }

//-----------------------------------------------------------------------------
// GetPowerSpectrum() calculates the power_spectrum with DC correction.
// DC stands for Direct Current. In this case, the component from 0 to F0 Hz
// is corrected.
//-----------------------------------------------------------------------------
    void CheapTrick::calculatePowerSpectrum() {
        int half_window_length = MatlabFunctions::matlab_round(1.5 * options.fs / current_f0);

        // FFT
        for (int i = half_window_length * 2 + 1; i < options.fft_size; ++i)
            forward_real_fft->waveform[i] = 0.0;
        fft_execute(forward_real_fft->forward_fft);

        // Calculation of the power spectrum.
        double *power_spectrum = forward_real_fft->waveform;
        for (int i = 0; i <= options.fft_size / 2; ++i)
            power_spectrum[i] =
                    forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] +
                    forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];

        // DC correction
        Common::DCCorrection(power_spectrum, current_f0, options.fs, options.fft_size, power_spectrum);
    }

//-----------------------------------------------------------------------------
// SetParametersForGetWindowedWaveform()
//-----------------------------------------------------------------------------
    void CheapTrick::generateWindow(int half_window_length, std::vector<double>& window, std::vector<int>& safe_index, std::vector<int>& base_index) const {
        for (int i = -half_window_length; i <= half_window_length; ++i) {
            base_index.at(i + half_window_length) = i;
        }

        int origin = MatlabFunctions::matlab_round(current_position * options.fs + 0.001);
        for (int i = 0; i <= half_window_length * 2; ++i)
            safe_index.at(i) = Common::MyMinInt(x_length - 1, Common::MyMaxInt(0, origin + base_index.at(i)));

        // Designing of the window function
        double average = 0.0;
        double position;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            position = base_index.at(i) / 1.5 / options.fs;
            window.at(i) = 0.5 * cos(world::kPi * position * current_f0) + 0.5;
            average += window.at(i) * window.at(i);
        }
        average = sqrt(average);
        for (int i = 0; i <= half_window_length * 2; ++i) {
            window.at(i) /= average;
        }
    }

//-----------------------------------------------------------------------------
// GetWindowedWaveform() windows the waveform by F0-adaptive window
//-----------------------------------------------------------------------------
    void CheapTrick::applyWindowing() {
        int half_window_length = MatlabFunctions::matlab_round(1.5 * options.fs / current_f0);

        auto base_index = std::vector<int>(half_window_length * 2 + 1);
        auto safe_index = std::vector<int>(half_window_length * 2 + 1);
        auto window = std::vector<double>(half_window_length * 2 + 1);

        generateWindow(half_window_length, window, base_index, safe_index);

        // F0-adaptive windowing
        double *waveform = forward_real_fft->waveform;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            waveform[i] = x->at(safe_index.at(i)) * window.at(i) + MatlabFunctions::randn() * world::kMySafeGuardMinimum;
        }

        double tmp_weight1 = 0;
        double tmp_weight2 = 0;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            tmp_weight1 += waveform[i];
            tmp_weight2 += window.at(i);
        }

        double weighting_coefficient = tmp_weight1 / tmp_weight2;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            waveform[i] -= window.at(i) * weighting_coefficient;
        }
    }

    // Add noise to the spectrum.
    void CheapTrick::addNoise() {
        for (int i = 0; i <= options.fft_size / 2; ++i)
            forward_real_fft->waveform[i] = forward_real_fft->waveform[i] + fabs(MatlabFunctions::randn()) * world::kEps;
    }


    // Calculate a spectral envelope at the current temporal position.
    void CheapTrick::step() {
        // F0-adaptive windowing
        applyWindowing();

        // Calculate power spectrum with DC correction
        // Note: The calculated power spectrum is stored in an array for waveform.
        // In this implementation, power spectrum is transformed by FFT (NOT IFFT).
        // However, the same result is obtained.
        // This is tricky but important for simple implementation.
        calculatePowerSpectrum();

        // Smoothing of the power (linear axis)
        // forward_real_fft.waveform is the power spectrum.
        // TODO: Bring this up to code with the rest of the updated library.
        Common::LinearSmoothing(forward_real_fft->waveform, current_f0 * 2.0 / 3.0,
                        options.fs, options.fft_size, forward_real_fft->waveform);

        // Add infinitesimal noise
        // This is a safeguard to avoid including zero in the spectrum.
        addNoise();

        // Smoothing (log axis) and spectral recovery on the cepstrum domain.
        smoothingWithRecovery();
    }

    void CheapTrick::initializeFFTs() {
        forward_real_fft = {nullptr};
        InitializeForwardRealFFT(options.fft_size, forward_real_fft);

        inverse_real_fft = {nullptr};
        InitializeInverseRealFFT(options.fft_size, inverse_real_fft);
    }


}  // namespace