//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//
// Voice synthesis based on f0, spectrogram and aperiodicity.
// forward_real_fft, inverse_real_fft and minimum_phase are used to speed up.
//-----------------------------------------------------------------------------
#include "TorchWorld/synthesis.h"

#include <math.h>

#include <ext/tqdm.h>
#include "TorchWorld/common.h"
#include "TorchWorld/constantnumbers.h"
#include "TorchWorld/matlabfunctions.h"

namespace tw {
    static void GetNoiseSpectrum(int noise_size, int fft_size, const Common::ForwardRealFFT *forward_real_fft) {
        double average = 0.0;
        for (int i = 0; i < noise_size; ++i) {
            forward_real_fft->waveform[i] = MatlabFunctions::randn();
            average += forward_real_fft->waveform[i];
        }

        average /= noise_size;
        for (int i = 0; i < noise_size; ++i)
            forward_real_fft->waveform[i] -= average;
        for (int i = noise_size; i < fft_size; ++i)
            forward_real_fft->waveform[i] = 0.0;
        fft_execute(forward_real_fft->forward_fft);
    }

//-----------------------------------------------------------------------------
// GetAperiodicResponse() calculates an aperiodic response.
//-----------------------------------------------------------------------------
    static void GetAperiodicResponse(int noise_size, int fft_size,
                                     const double *spectrum, const double *aperiodic_ratio, double current_vuv,
                                     const Common::ForwardRealFFT *forward_real_fft,
                                     const Common::InverseRealFFT *inverse_real_fft,
                                     const Common::MinimumPhaseAnalysis *minimum_phase, double *aperiodic_response) {
        GetNoiseSpectrum(noise_size, fft_size, forward_real_fft);

        if (current_vuv != 0.0)
            for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
                minimum_phase->log_spectrum[i] =
                        log(spectrum[i] * aperiodic_ratio[i]) / 2.0;
        else
            for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
                minimum_phase->log_spectrum[i] = log(spectrum[i]) / 2.0;
        GetMinimumPhaseSpectrum(minimum_phase);

        for (int i = 0; i <= fft_size / 2; ++i) {
            inverse_real_fft->spectrum[i][0] =
                    minimum_phase->minimum_phase_spectrum[i][0] *
                    forward_real_fft->spectrum[i][0] -
                    minimum_phase->minimum_phase_spectrum[i][1] *
                    forward_real_fft->spectrum[i][1];
            inverse_real_fft->spectrum[i][1] =
                    minimum_phase->minimum_phase_spectrum[i][0] *
                    forward_real_fft->spectrum[i][1] +
                    minimum_phase->minimum_phase_spectrum[i][1] *
                    forward_real_fft->spectrum[i][0];
        }
        fft_execute(inverse_real_fft->inverse_fft);
        MatlabFunctions::fftshift(inverse_real_fft->waveform, fft_size, aperiodic_response);
    }

//-----------------------------------------------------------------------------
// RemoveDCComponent()
//-----------------------------------------------------------------------------
    static void RemoveDCComponent(const double *periodic_response, int fft_size,
                                  const double *dc_remover, double *new_periodic_response) {
        double dc_component = 0.0;
        for (int i = fft_size / 2; i < fft_size; ++i)
            dc_component += periodic_response[i];
        for (int i = 0; i < fft_size / 2; ++i)
            new_periodic_response[i] = -dc_component * dc_remover[i];
        for (int i = fft_size / 2; i < fft_size; ++i)
            new_periodic_response[i] -= dc_component * dc_remover[i];
    }

//-----------------------------------------------------------------------------
// GetSpectrumWithFractionalTimeShift() calculates a periodic spectrum with
// the fractional time shift under 1/fs.
//-----------------------------------------------------------------------------
    static void GetSpectrumWithFractionalTimeShift(int fft_size, double coefficient, const Common::InverseRealFFT *inverse_real_fft) {
        double re, im, re2, im2;
        for (int i = 0; i <= fft_size / 2; ++i) {
            re = inverse_real_fft->spectrum[i][0];
            im = inverse_real_fft->spectrum[i][1];
            re2 = cos(coefficient * i);
            im2 = sqrt(1.0 - re2 * re2);  // sin(pshift)

            inverse_real_fft->spectrum[i][0] = re * re2 + im * im2;
            inverse_real_fft->spectrum[i][1] = im * re2 - re * im2;
        }
    }

//-----------------------------------------------------------------------------
// GetPeriodicResponse() calculates a periodic response.
//-----------------------------------------------------------------------------
    static void GetPeriodicResponse(int fft_size, const double *spectrum, const double *aperiodic_ratio, double current_vuv, const Common::InverseRealFFT *inverse_real_fft, const Common::MinimumPhaseAnalysis *minimum_phase, const double *dc_remover, double fractional_time_shift, int fs, double *periodic_response) {
        if (current_vuv <= 0.5 || aperiodic_ratio[0] > 0.999) {
            for (int i = 0; i < fft_size; ++i) periodic_response[i] = 0.0;
            return;
        }

        for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
            minimum_phase->log_spectrum[i] =
                    log(spectrum[i] * (1.0 - aperiodic_ratio[i]) +
                        world::kMySafeGuardMinimum) / 2.0;
        GetMinimumPhaseSpectrum(minimum_phase);

        for (int i = 0; i <= fft_size / 2; ++i) {
            inverse_real_fft->spectrum[i][0] =
                    minimum_phase->minimum_phase_spectrum[i][0];
            inverse_real_fft->spectrum[i][1] =
                    minimum_phase->minimum_phase_spectrum[i][1];
        }

        // apply fractional time delay of fractional_time_shift seconds
        // using linear phase shift
        double coefficient =
                2.0 * world::kPi * fractional_time_shift * fs / fft_size;
        GetSpectrumWithFractionalTimeShift(fft_size, coefficient, inverse_real_fft);

        fft_execute(inverse_real_fft->inverse_fft);
        MatlabFunctions::fftshift(inverse_real_fft->waveform, fft_size, periodic_response);
        RemoveDCComponent(periodic_response, fft_size, dc_remover,
                          periodic_response);
    }

    static void GetSpectralEnvelope(double current_time, double frame_period,
                                    int f0_length, const torch::Tensor &spectrogram, int fft_size,
                                    double *spectral_envelope) {
        int current_frame_floor = Common::MyMinInt(f0_length - 1,
                                           static_cast<int>(floor(current_time / frame_period)));
        int current_frame_ceil = Common::MyMinInt(f0_length - 1,
                                          static_cast<int>(ceil(current_time / frame_period)));
        double interpolation = current_time / frame_period - current_frame_floor;

        if (current_frame_floor == current_frame_ceil)
            for (int i = 0; i <= fft_size / 2; ++i)
                spectral_envelope[i] = fabs(spectrogram[current_frame_floor][i].item<float>());
        else
            for (int i = 0; i <= fft_size / 2; ++i)
                spectral_envelope[i] =
                        (1.0 - interpolation) * fabs(spectrogram[current_frame_floor][i].item<float>()) +
                        interpolation * fabs(spectrogram[current_frame_ceil][i].item<float>());
    }

    static void GetAperiodicRatio(double current_time, double frame_period,
                                  int f0_length, const torch::Tensor &aperiodicity, int fft_size,
                                  double *aperiodic_spectrum) {
        int current_frame_floor = Common::MyMinInt(f0_length - 1, static_cast<int>(floor(current_time / frame_period)));
        int current_frame_ceil = Common::MyMinInt(f0_length - 1,
                                          static_cast<int>(ceil(current_time / frame_period)));
        double interpolation = current_time / frame_period - current_frame_floor;

        if (current_frame_floor == current_frame_ceil)
            for (int i = 0; i <= fft_size / 2; ++i)
                aperiodic_spectrum[i] =
                        pow(Common::GetSafeAperiodicity(aperiodicity[current_frame_floor][i].item<float>()), 2.0);
        else
            for (int i = 0; i <= fft_size / 2; ++i)
                aperiodic_spectrum[i] = pow((1.0 - interpolation) * Common::GetSafeAperiodicity(aperiodicity[current_frame_floor][i].item<float>()) + interpolation * Common::GetSafeAperiodicity(aperiodicity[current_frame_ceil][i].item<float>()), 2.0);
    }

//-----------------------------------------------------------------------------
// GetOneFrameSegment() calculates a periodic and aperiodic response at a time.
//-----------------------------------------------------------------------------
    void Synthesis::step() {
        std::vector<double> aperiodic_response(options.fft_size);
        std::vector<double> periodic_response(options.fft_size);

        std::vector<double> spectral_envelope(options.fft_size);
        std::vector<double> aperiodic_ratio(options.fft_size);

        GetSpectralEnvelope(current_time, options.frame_period, options.f0_length, spectrogram,
                            options.fft_size, spectral_envelope.data());
        GetAperiodicRatio(current_time, options.frame_period, options.f0_length, aperiodicity,
                          options.fft_size, aperiodic_ratio.data());

        // Synthesis of the periodic response
        GetPeriodicResponse(options.fft_size, spectral_envelope.data(), aperiodic_ratio.data(),
                            current_vuv, inverse_real_fft, minimum_phase, dc_remover.data(),
                            fractional_time_shift, options.fs, periodic_response.data());

        // Synthesis of the aperiodic response
        GetAperiodicResponse(noise_size, options.fft_size, spectral_envelope.data(),
                             aperiodic_ratio.data(), current_vuv, forward_real_fft,
                             inverse_real_fft, minimum_phase, aperiodic_response.data());

        double sqrt_noise_size = sqrt(static_cast<double>(noise_size));
        for (int i = 0; i < options.fft_size; ++i)
            impulse_response[i] = (periodic_response[i] * sqrt_noise_size + aperiodic_response[i]) /
                    options.fft_size;
    }

    void Synthesis::getTemporalParameters() {
        for (int i = 0; i < options.y_length; ++i) {
            time_axis[i] = i / static_cast<double>(options.fs);
        }

        // the array 'coarse_time_axis' is supposed to have 'f0_length + 1' positions
        for (int i = 0; i < options.f0_length; ++i) {
            coarse_time_axis[i] = i * options.frame_period;
            coarse_f0[i] = f0[i].item<double>() < lowest_f0 ? 0.0 : f0[i].item<double>();
            coarse_vuv[i] = coarse_f0[i] == 0.0 ? 0.0 : 1.0;
        }

        coarse_time_axis[options.f0_length] = options.f0_length * options.frame_period;

        coarse_f0[options.f0_length] = coarse_f0[options.f0_length - 1] * 2 - coarse_f0[options.f0_length - 2];

        coarse_vuv[options.f0_length] = coarse_vuv[options.f0_length - 1] * 2 - coarse_vuv[options.f0_length - 2];
    }

    int Synthesis::getPulseLocations() {
        std::vector<double> total_phase(options.y_length);
        std::vector<double> wrap_phase(options.y_length);
        std::vector<double> wrap_phase_abs(options.y_length - 1);

        total_phase[0] = 2.0 * world::kPi * interpolated_f0[0] / options.fs;
        wrap_phase[0] = fmod(total_phase[0], 2.0 * world::kPi);

        for (int i = 1; i < options.y_length; ++i) {
            total_phase[i] = total_phase[i - 1] + 2.0 * world::kPi * interpolated_f0[i] / options.fs;
            wrap_phase[i] = fmod(total_phase[i], 2.0 * world::kPi);
            wrap_phase_abs[i - 1] = fabs(wrap_phase[i] - wrap_phase[i - 1]);
        }

        int number_of_pulses = 0;
        for (int i = 0; i < options.y_length - 1; ++i) {
            if (wrap_phase_abs[i] > world::kPi) {
                pulse_locations[number_of_pulses] = time_axis[i];
                pulse_locations_index[number_of_pulses] = i;

                // calculate the time shift in seconds between exact fractional pulse
                // position and the integer pulse position (sample i)
                // as we don't have access to the exact pulse position, we infer it
                // from the point between sample i and sample i + 1 where the
                // accumulated phase cross a multiple of 2pi
                // this point is found by solving y1 + x * (y2 - y1) = 0 for x, where y1
                // and y2 are the phases corresponding to sample i and i + 1, offset so
                // they cross zero; x >= 0
                double y1 = wrap_phase[i] - 2.0 * world::kPi;
                double y2 = wrap_phase[i + 1];
                double x = -y1 / (y2 - y1);
                pulse_locations_time_shift[number_of_pulses] = x / options.fs;

                ++number_of_pulses;
            }
        }

        return number_of_pulses;
    }

    int Synthesis::getTimeBase() {
        time_axis = std::vector<double>(options.y_length);
        coarse_time_axis = std::vector<double>(options.f0_length + 1);
        coarse_f0 = std::vector<double>(options.f0_length + 1);
        coarse_vuv = std::vector<double>(options.f0_length + 1);
        interpolated_f0 = std::vector<double>(options.y_length);

        getTemporalParameters();

        MatlabFunctions::interp1(coarse_time_axis.data(), coarse_f0.data(), options.f0_length + 1,
                time_axis.data(), options.y_length, interpolated_f0.data());
        MatlabFunctions::interp1(coarse_time_axis.data(), coarse_vuv.data(), options.f0_length + 1,
                time_axis.data(), options.y_length, interpolated_vuv.data());

        for (int i = 0; i < options.y_length; ++i) {
            interpolated_vuv[i] = interpolated_vuv[i] > 0.5 ? 1.0 : 0.0;
            interpolated_f0[i] = interpolated_vuv[i] == 0.0 ? world::kDefaultF0 : interpolated_f0[i];
        }

        return getPulseLocations();
    }

    static void GetDCRemover(int fft_size, std::vector<double>& dc_remover) {
        double dc_component = 0.0;

        for (int i = 0; i < fft_size / 2; ++i) {
            dc_remover[i] = 0.5 - 0.5 * cos(2.0 * world::kPi * (i + 1.0) / (1.0 + fft_size));
            dc_remover[fft_size - i - 1] = dc_remover[i];
            dc_component += dc_remover[i] * 2.0;
        }
        for (int i = 0; i < fft_size / 2; ++i) {
            dc_remover[i] /= dc_component;
            dc_remover[fft_size - i - 1] = dc_remover[i];
        }
    }

    void Synthesis::initializeFFTs() {
        minimum_phase = {0};
        InitializeMinimumPhaseAnalysis(options.fft_size, minimum_phase);

        inverse_real_fft = {0};
        InitializeInverseRealFFT(options.fft_size, inverse_real_fft);

        forward_real_fft = {0};
        InitializeForwardRealFFT(options.fft_size, forward_real_fft);
    }

    void Synthesis::run() {
        MatlabFunctions::randn_reseed();
        initializeFFTs();

        lowest_f0 = options.fs / options.fft_size + 1.0;

        impulse_response = std::vector<double>(options.y_length, 0.0);
        interpolated_vuv = std::vector<double>(options.y_length, 0.0);

        pulse_locations = std::vector<double>(options.y_length, 0.0);
        pulse_locations_index = std::vector<int>(options.y_length, 0);
        pulse_locations_time_shift = std::vector<double>(options.y_length, 0.0);

        int number_of_pulses = getTimeBase();

        dc_remover = std::vector<double>(options.fft_size, 0.0);
        GetDCRemover(options.fft_size, dc_remover);

        options.frame_period /= 1000.0;
        int index, offset, lower_limit, upper_limit;

        for (int i: tqdm::range(number_of_pulses)) {
            noise_size = pulse_locations_index[Common::MyMinInt(number_of_pulses - 1, i + 1)] - pulse_locations_index[i];

            step();

            offset = pulse_locations_index[i] - options.fft_size / 2 + 1;
            lower_limit = Common::MyMaxInt(0, -offset);
            upper_limit = Common::MyMinInt(options.fft_size, options.y_length - offset);

            for (int j = lower_limit; j < upper_limit; ++j) {
                index = j + offset;
                y->at(index) += impulse_response[j];
            }
        }

        DestroyMinimumPhaseAnalysis(minimum_phase);
        DestroyInverseRealFFT(inverse_real_fft);
        DestroyForwardRealFFT(forward_real_fft);
    }

}  // namespace