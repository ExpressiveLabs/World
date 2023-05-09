//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_DIO_H_
#define TORCHWORLD_DIO_H_

#include "TorchWorld/macrodefinitions.h"
#include "fft.h"

namespace tw {
    class Dio {
        public:
    //-----------------------------------------------------------------------------
    // Struct for DIO
    //-----------------------------------------------------------------------------
            struct Options {
                double f0_floor;
                double f0_ceil;
                double channels_in_octave;
                double frame_period;  // msec
                int speed;  // (1, 2, ..., 12)
                double allowed_range;  // Threshold used for fixing the F0 contour.

                Options();
            };

            //-----------------------------------------------------------------------------
            // struct for GetFourZeroCrossingIntervals()
            // "negative" means "zero-crossing point going from positive to negative"
            // "positive" means "zero-crossing point going from negative to positive"
            //-----------------------------------------------------------------------------
            struct ZeroCrossings {
                double *negative_interval_locations;
                double *negative_intervals;
                int number_of_negatives;
                double *positive_interval_locations;
                double *positive_intervals;
                int number_of_positives;
                double *peak_interval_locations;
                double *peak_intervals;
                int number_of_peaks;
                double *dip_interval_locations;
                double *dip_intervals;
                int number_of_dips;
            };

    //-----------------------------------------------------------------------------
    // DIO
    //
    // Input:
    //   x                    : Input signal
    //   x_length             : Length of x
    //   fs                   : Sampling frequency
    //   option               : Struct to order the parameter for DIO
    //
    // Output:
    //   temporal_positions   : Temporal positions.
    //   f0                   : F0 contour.
    //-----------------------------------------------------------------------------
            Dio(const double *x, int x_length, int fs, const Options *option,
                     double *temporal_positions, double *f0);
            ~Dio() = default;

            static void GeneralBody(const double *x, int x_length, int fs, double frame_period, double f0_floor, double f0_ceil, double channels_in_octave, int speed, double allowed_range, double *temporal_positions, double *f0);
    //-----------------------------------------------------------------------------
    // InitializeDioOption allocates the memory to the struct and sets the
    // default parameters.
    //
    // Output:
    //   option   : Struct for the optional parameter.
    //-----------------------------------------------------------------------------

    //-----------------------------------------------------------------------------
    // GetSamplesForDIO() calculates the number of samples required for Dio().
    //
    // Input:
    //   fs             : Sampling frequency [Hz]
    //   x_length       : Length of the input signal [Sample].
    //   frame_period   : Frame shift [msec]
    //
    // Output:
    //   The number of samples required to store the results of Dio()
    //-----------------------------------------------------------------------------
            static int GetSamplesForDIO(int fs, int x_length, double frame_period);
            static void DesignLowCutFilter(int N, int fft_size, double *low_cut_filter);
            static void GetFilteredSignal(int half_average_length, int fft_size, const FFT::fft_complex *y_spectrum, int y_length, double *filtered_signal);
            static void GetSpectrumForEstimation(const double *x, int x_length, int y_length, double actual_fs, int fft_size, int decimation_ratio, FFT::fft_complex *y_spectrum);

            static void GetF0CandidateFromRawEvent(double boundary_f0, double fs, const FFT::fft_complex *y_spectrum, int y_length, int fft_size, double f0_floor, double f0_ceil, const double *temporal_positions, int f0_length, double *f0_score, double *f0_candidate);
            static void GetF0CandidatesAndScores(const double *boundary_f0_list, int number_of_bands, double actual_fs, int y_length, const double *temporal_positions, int f0_length, const FFT::fft_complex *y_spectrum, int fft_size, double f0_floor, double f0_ceil, double **raw_f0_candidates, double **raw_f0_scores);

            static int ZeroCrossingEngine(const double *filtered_signal, int y_length, double fs, double *interval_locations, double *intervals);
            static void GetFourZeroCrossingIntervals(double *filtered_signal, int y_length, double actual_fs, ZeroCrossings *zero_crossings);
            static void DestroyZeroCrossings(ZeroCrossings *zero_crossings);

            static void GetF0CandidateContour(const ZeroCrossings *zero_crossings, double boundary_f0, double f0_floor, double f0_ceil, const double *temporal_positions, int f0_length, double *f0_candidate, double *f0_score);
            static void GetF0CandidateContourSub(const double *const *interpolated_f0_set, int f0_length, double f0_floor, double f0_ceil, double boundary_f0, double *f0_candidate, double *f0_score);
            static void FixF0Contour(double frame_period, int number_of_candidates, int fs, const double *const *f0_candidates, const double *best_f0_contour, int f0_length, double f0_floor, double allowed_range, double *fixed_f0_contour);
            static void GetBestF0Contour(int f0_length, const double *const *f0_candidates, const double *const *f0_scores, int number_of_bands, double *best_f0_contour);

            static void FixStep4(const double *f0_step3, int f0_length, const double *const *f0_candidates, int number_of_candidates, double allowed_range, const int *positive_index, int positive_count, double *f0_step4);
            static void FixStep3(const double *f0_step2, int f0_length, const double *const *f0_candidates, int number_of_candidates, double allowed_range, const int *negative_index, int negative_count, double *f0_step3);
            static void FixStep2(const double *f0_step1, int f0_length, int voice_range_minimum, double *f0_step2);
            static void FixStep1(const double *best_f0_contour, int f0_length, int voice_range_minimum, double allowed_range, double *f0_step1);
            static double SelectBestF0(double current_f0, double past_f0, const double *const *f0_candidates, int number_of_candidates, int target_index, double allowed_range);

            static void GetNumberOfVoicedSections(const double *f0, int f0_length, int *positive_index, int *negative_index, int *positive_count, int *negative_count);
    };
}

#endif  // WORLD_DIO_H_
