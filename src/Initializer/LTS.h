/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/
 
#ifndef INITIALIZER_LTS_H_
#define INITIALIZER_LTS_H_

#include <Initializer/typedefs.hpp>
#include <generated_code/tensor.h>
#include "tree/Common.hpp"

#if CONVERGENCE_ORDER < 2 || CONVERGENCE_ORDER > 8
#error Preprocessor flag CONVERGENCE_ORDER is not in {2, 3, 4, 5, 6, 7, 8}.
#endif

/* ORIGINAL CODE
#   define MEMKIND_GLOBAL   seissol::memory::HighBandwidth
#if CONVERGENCE_ORDER <= 7
#   define MEMKIND_TIMEDOFS seissol::memory::HighBandwidth
#else
#   define MEMKIND_TIMEDOFS seissol::memory::Standard
#endif
#if CONVERGENCE_ORDER <= 4
#   define MEMKIND_CONSTANT seissol::memory::HighBandwidth
#else
#   define MEMKIND_CONSTANT seissol::memory::Standard
#endif
#if CONVERGENCE_DOFS <= 3
#   define MEMKIND_DOFS     seissol::memory::HighBandwidth
#else
#   define MEMKIND_DOFS     seissol::memory::Standard
#endif
*/

//----------------------------------------------------------------
// DEBUGING: porting to gpu
#   define MEMKIND_GLOBAL   seissol::memory::Standard  //debugging: must be DeviceGlobalMemory at the end
#if CONVERGENCE_ORDER <= 7
#   define MEMKIND_TIMEDOFS seissol::memory::Standard  //debugging: must be both Standard and DeviceGlobalMemory at the end
#else
#   define MEMKIND_TIMEDOFS seissol::memory::Standard
#endif
#if CONVERGENCE_ORDER <= 4
#   define MEMKIND_CONSTANT seissol::memory::Standard
#else
#   define MEMKIND_CONSTANT seissol::memory::Standard
#endif
#if CONVERGENCE_DOFS <= 3
#   define MEMKIND_DOFS     seissol::memory::Standard
#else
#   define MEMKIND_DOFS     seissol::memory::Standard
#endif
//----------------------------------------------------------------


namespace seissol {
  namespace initializers {
    struct LTS;
    class  LTSTree;
  }
}

struct seissol::initializers::LTS {

  // define all variables needed for computataions within one element
  seissol::initializers::Variable<real[tensor::Q::size()]>       dofs;
#if NUMBER_OF_RELAXATION_MECHANISMS > 0
  seissol::initializers::Variable<real[tensor::Qane::size()]>    dofsAne;
#endif
  seissol::initializers::Variable<real*>                         buffers;
  seissol::initializers::Variable<real*>                         derivatives;
  seissol::initializers::Variable<CellLocalInformation>          cellInformation;
  seissol::initializers::Variable<real*[4]>                      faceNeighbors;
  seissol::initializers::Variable<LocalIntegrationData>          localIntegration;
  seissol::initializers::Variable<NeighboringIntegrationData>    neighboringIntegration;
  seissol::initializers::Variable<CellMaterialData>              material;
  seissol::initializers::Variable<PlasticityData>                plasticity;
  seissol::initializers::Variable<CellDRMapping[4]>              drMapping;
  seissol::initializers::Variable<real[3]>                       energy;
  seissol::initializers::Variable<real[7]>                       pstrain;
  seissol::initializers::Variable<real*>                         displacements;
  seissol::initializers::Bucket                                  buffersDerivatives;
  seissol::initializers::Bucket                                  displacementsBuffer;

  void addTo(LTSTree& tree);
};
#endif
