//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef TORCHWORLD_STONEMASK_H_
#define TORCHWORLD_STONEMASK_H_

#include "TorchWorld/macrodefinitions.h"

namespace tw {
    TW_WORLD_BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// StoneMask() refines the estimated F0 by Dio()
//
// Input:
//   x                      : Input signal
//   x_length               : Length of the input signal
//   fs                     : Sampling frequency
//   time_axis              : Temporal information
//   f0                     : f0 contour
//   f0_length              : Length of f0
//
// Output:
//   refined_f0             : Refined F0
//-----------------------------------------------------------------------------
    void StoneMask(const double *x, int x_length, int fs,
                   const double *temporal_positions, const double *f0, int f0_length,
                   double *refined_f0);

    TW_WORLD_END_C_DECLS
}

#endif  // WORLD_STONEMASK_H_
