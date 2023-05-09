//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_SYNTHESISREALTIME_H_
#define TORCHWORLD_SYNTHESISREALTIME_H_

#include "TorchWorld/common.h"
#include "TorchWorld/macrodefinitions.h"

namespace tw {
    class SynthesisRealtime {
        public:

    //-----------------------------------------------------------------------------
    // A struct for real-time synthesis.
    // I use no class for compatibility in C language.
    // Please make a class for encapsulating it as needed.
    // This synthesizer uses a ring buffer.
    //-----------------------------------------------------------------------------
            struct WorldSynthesizer {
                // Basic parameters
                int fs;
                double frame_period;
                int buffer_size;
                int number_of_pointers;
                int fft_size;

                // Sound buffer for output. The length is buffer_size [sample].
                double *buffer;
                int current_pointer;
                int i;

                // For DC removal
                double *dc_remover;

                //---------------------------------------------------------------------------
                // Followings are internal parameters.
                // You should not modify them if you are not expert.

                // Speech parameters in each pointer.
                int *f0_length;
                int *f0_origin;
                double ***spectrogram;
                double ***aperiodicity;


                // Note:
                // This is an extremely rough implementation.
                // I should optimize this implementation.
                int current_pointer2;
                int head_pointer;
                int synthesized_sample;

                // Internal parameters.
                int handoff;
                double handoff_phase;
                double handoff_f0;
                int last_location;

                int cumulative_frame;
                int current_frame;

                double **interpolated_vuv;
                double **pulse_locations;
                int **pulse_locations_index;
                int *number_of_pulses;

                double *impulse_response;

                // FFT
                Common::MinimumPhaseAnalysis minimum_phase;
                Common::InverseRealFFT inverse_real_fft;
                Common::ForwardRealFFT forward_real_fft;

                WorldSynthesizer(int fs, double frame_period, int fft_size, int buffer_size, int number_of_pointers);
                ~WorldSynthesizer();

                void refresh();
                void clearRingBuffer();
            };

    //-----------------------------------------------------------------------------
    // AddParameters() attempts to add speech parameters.
    // You can add several frames at the same time.
    //
    // Input:
    //   f0                   : F0 contour with length of f0_length
    //   f0_length            : This is associated with the number of frames
    //   spectrogram          : Spectrogram
    //   aperiodicity         : Aperiodicity
    //
    // Output:
    //   synth                : Synthesizer
    //
    // Return value:
    //   1: True, 0: False.
    //-----------------------------------------------------------------------------
            static int AddParameters(double *f0, int f0_length, double **spectrogram, double **aperiodicity, WorldSynthesizer *synth);

    //-----------------------------------------------------------------------------
    // RefreshSynthesizer() sets the parameters to default.
    //-----------------------------------------------------------------------------
            static void RefreshSynthesizer(WorldSynthesizer *synth);

    //-----------------------------------------------------------------------------
    // DestroySynthesizer() release the memory.
    //-----------------------------------------------------------------------------
            static void DestroySynthesizer(WorldSynthesizer *synth);

    //-----------------------------------------------------------------------------
    // IsLocked() checks whether the synthesizer is locked or not.
    // "Lock" is defined as the situation that the ring buffer cannot add
    // parameters and cannot synthesize the waveform.
    // It will be caused when the duration calculated by the number of added frames
    // is below 1 / F0 + buffer_size / fs.
    // If this function returns True, please refresh the synthesizer.
    //
    // Input:
    //   Synth            : Synthesizer (pointer)
    //
    // Output:
    //   1: True, 0: False.
    //-----------------------------------------------------------------------------
            static int IsLocked(WorldSynthesizer *synth);
            static int Synthesis2(WorldSynthesizer *synth);

            static void ClearRingBuffer(int start, int end, WorldSynthesizer *synth);
            static int SeekSynthesizer(double current_location, WorldSynthesizer *synth);
            static void SearchPointer(int frame, WorldSynthesizer *synth, int flag, double **front, double **next);

            static void GetDCRemover(int fft_size, double *dc_remover);
            static void RemoveDCComponent(const double *periodic_response, int fft_size, const double *dc_remover, double *new_periodic_response);

            static void GetPeriodicResponse(int fft_size, const double *spectrum, const double *aperiodic_ratio, double current_vuv, const Common::InverseRealFFT *inverse_real_fft, const Common::MinimumPhaseAnalysis *minimum_phase, const double *dc_remover, double *periodic_response);
            static void GetSpectralEnvelope(double current_location, WorldSynthesizer *synth, double *spectral_envelope);
            static void GetAperiodicRatio(double current_location, WorldSynthesizer *synth, double *aperiodic_spectrum);
            static double GetCurrentVUV(int current_location, WorldSynthesizer *synth);

            static void GetOneFrameSegment(int noise_size, int current_location, WorldSynthesizer *synth);
            static void GetTemporalParametersForTimeBase(const double *f0, int f0_length, WorldSynthesizer *synth, double *coarse_time_axis, double *coarse_f0, double *coarse_vuv);
            static void GetPulseLocationsForTimeBase(const double *interpolated_f0, const double *time_axis, int number_of_samples, double origin, WorldSynthesizer *synth);
            static void GetTimeBase(const double *f0, int f0_length, int start_sample, int number_of_samples, WorldSynthesizer *synth);
            static int GetNextPulseLocationIndex(WorldSynthesizer *synth);

            static int UpdateSynthesizer(int current_location, WorldSynthesizer *synth);
            static int CheckSynthesizer(WorldSynthesizer *synth);
    };
}

#endif  // WORLD_SYNTHESISREALTIME_H_
