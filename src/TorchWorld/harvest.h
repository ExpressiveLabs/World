//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_HARVEST_H_
#define TORCHWORLD_HARVEST_H_

#include <TorchWorld/macrodefinitions.h>
#include <TorchWorld/constantnumbers.h>

namespace tw {
    class Harvest {
    public:

//-----------------------------------------------------------------------------
// Struct for Harvest
//-----------------------------------------------------------------------------
        struct Options {
            double f0_floor;
            double f0_ceil;
            double frame_period;

            Options() {
                f0_ceil = world::kCeilF0;
                f0_floor = world::kFloorF0;
                frame_period = 5;
            }
        };

//-----------------------------------------------------------------------------
// Harvest
//
// Input:
//   x                    : Input signal
//   x_length             : Length of x
//   fs                   : Sampling frequency
//   option               : Struct to order the parameter for Harvest
//
// Output:
//   temporal_positions   : Temporal positions.
//   f0                   : F0 contour.
//-----------------------------------------------------------------------------
        Harvest(const double *x, int x_length, int fs, const Options *option, double *temporal_positions, double *f0);

        static void GeneralBody(const double *x, int x_length, int fs, int frame_period, double f0_floor, double f0_ceil, double channels_in_octave, int speed, double *temporal_positions, double *f0);

        static int GetSamplesForHarvest(int fs, int x_length, double frame_period);

    };
}

#endif  // TORCHWORLD_HARVEST_H_
