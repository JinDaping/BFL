// $Id$
// Copyright (C) 2003 Klaas Gadeyne <first dot last at gmail dot com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
//

#include "conditionalgaussian_additivenoise.h"
#include <cmath>
#include "../wrappers/rng/rng.h"

namespace BFL
{

  using namespace MatrixWrapper;

#define CondGausAddNoise ConditionalGaussianAdditiveNoise

  /// Constructor
  /** @note Because of the virtual public inheritance, we create the Pdf
      baseclass ourselves
  */
  CondGausAddNoise::CondGausAddNoise(const Gaussian& additiveNoise,
				     int num_conditional_arguments)
    : ConditionalGaussian(additiveNoise.DimensionGet(),
			  num_conditional_arguments)
    , _additiveNoise_Mu   (additiveNoise.ExpectedValueGet())
    , _additiveNoise_Sigma(additiveNoise.CovarianceGet())
  {}

  CondGausAddNoise::CondGausAddNoise(int dim,
				     int num_conditional_arguments)
    : ConditionalGaussian(dim, num_conditional_arguments)
  {
    _additiveNoise_Mu.resize(dim);
    _additiveNoise_Sigma.resize(dim);
  }

  /// Destructor
  CondGausAddNoise::~CondGausAddNoise()
  {}

  SymmetricMatrix
  CondGausAddNoise::CovarianceGet() const
  {
    return AdditiveNoiseSigmaGet();
  }

  const ColumnVector&
  CondGausAddNoise::AdditiveNoiseMuGet() const
  {
    return _additiveNoise_Mu;
  }

  const SymmetricMatrix&
  CondGausAddNoise::AdditiveNoiseSigmaGet() const
  {
    return _additiveNoise_Sigma;
  }

  void
  CondGausAddNoise::AdditiveNoiseMuSet(const ColumnVector& mu)
  {
    _additiveNoise_Mu = mu;
  }

  void
  CondGausAddNoise::AdditiveNoiseSigmaSet(const SymmetricMatrix& sigma)
  {
    _additiveNoise_Sigma = sigma;
  }

} // End namespace
