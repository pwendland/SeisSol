/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Local kernel of SeisSol.
 **/

#include "Kernels/Local.h"

#ifndef NDEBUG
#pragma message "compiling local kernel with assertions"
#endif

#include <yateto.h>


#include <cassert>
#include <stdint.h>
#include <cstring>

//DEBUGGING
#include <functional>
#include <iostream>
#include <stdio.h>
#include <vector>


void seissol::kernels::Local::setGlobalData(GlobalData const* global) {
#ifndef NDEBUG
  for (unsigned stiffness = 0; stiffness < 3; ++stiffness) {
    assert( ((uintptr_t)global->stiffnessMatrices(stiffness)) % ALIGNMENT == 0 );
  }
  for (unsigned flux = 0; flux < 4; ++flux) {
    assert( ((uintptr_t)global->localChangeOfBasisMatricesTransposed(flux)) % ALIGNMENT == 0 );
    assert( ((uintptr_t)global->changeOfBasisMatrices(flux)) % ALIGNMENT == 0 );
  }
#endif

  m_volumeKernelPrototype.kDivM = global->stiffnessMatrices;
  m_localFluxKernelPrototype.rDivM = global->changeOfBasisMatrices;
  m_localFluxKernelPrototype.fMrT = global->localChangeOfBasisMatricesTransposed;
}


void seissol::kernels::Local::computeIntegralModified(real *i_timeIntegratedDegreesOfFreedom[tensor::I::size()],
                                                      LocalData::Loader& loader,
                                                      seissol::initializers::DeviceVarInfo& manager,
                                                      const size_t num_cells,
                                                      LocalTmp& )
{

  auto data = loader.entry(0);
  // assert alignments
#ifndef NDEBUG
  //assert( ((uintptr_t)i_timeIntegratedDegreesOfFreedom) % ALIGNMENT == 0 );
  //assert( ((uintptr_t)data.dofs)              % ALIGNMENT == 0 );
#endif

  // compute the volume kernel on GPU
  kernel::volume volKrnl = m_volumeKernelPrototype;
  volKrnl.star(0) = manager.getDevicePointer(STARS_ID);
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    volKrnl.star(i) = volKrnl.star(i - 1) + tensor::star::Size[i];
  }

  volKrnl.Q = manager.getDevicePointer(DOFS_ID);
  volKrnl.I = manager.getDevicePointer(INTEGRATED_DOFS_ID);

  tensor::num_elements_in_cluster = num_cells;
  volKrnl.execute();



  //-------------------------------------------------------------------------------------------------------------------
  size_t AplusT_size = 4 * tensor::AplusT::Size;
  real *d_AplusT = manager.getDevicePointer(APLUST_ID);

  kernel::localFlux lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Q = manager.getDevicePointer(DOFS_ID);
  lfKrnl.I = manager.getDevicePointer(INTEGRATED_DOFS_ID);
  const unsigned num_faces = 4;
  unsigned *QI_strides[4];
  unsigned *AplusT_strides[4];
  for (unsigned face = 0; face < num_faces; ++face) {
    const std::vector<unsigned> bin = manager.getFaceElementIndices(LOCAL_FLUX, face);
    if (!bin.empty()) {
      // allocate indices in a device
      QI_strides[face] = (unsigned *) device_malloc(
          bin.size() * sizeof(unsigned));    // TODO: remove allocation of memory
      device_copy_to(QI_strides[face], bin.data(), bin.size() * sizeof(unsigned));
      device_scale_array(tensor::Q::Size, bin.size(), QI_strides[face]);

      AplusT_strides[face] = (unsigned *) device_malloc(
          bin.size() * sizeof(unsigned));    // TODO: remove allocation of memory
      device_copy_to(AplusT_strides[face], bin.data(), bin.size() * sizeof(unsigned));
      device_scale_array(AplusT_size, bin.size(), AplusT_strides[face]);

      // attach strides to the kernel
      tensor::Q::stride_ptr = QI_strides[face];
      tensor::I::stride_ptr = QI_strides[face];
      tensor::AplusT::stride_ptr = AplusT_strides[face];
      tensor::fMrT::stride_ptr[face] = seissol::tensor::device_zeros;
      tensor::rDivM::stride_ptr[face] = seissol::tensor::device_zeros;

      // attach a flux solver to the kernel
      lfKrnl.AplusT = d_AplusT + (face * tensor::AplusT::Size);

      // specify the number of faces to compute
      tensor::num_elements_in_cluster = bin.size();

      // execute the kernel for a particular face
      lfKrnl.execute(face);
    }
  }

  // release resources
  for (unsigned face = 0; face < num_faces; ++face) {
    device_free((void*)QI_strides[face]);  // TODO: remove de-allocation of memory
    device_free((void*)AplusT_strides[face]);  // TODO: remove de-allocation of memory
  }
}


void seissol::kernels::Local::computeIntegral(  real       i_timeIntegratedDegreesOfFreedom[tensor::I::size()],
                                                LocalData& data,
                                                LocalTmp& )
{
  // assert alignments
#ifndef NDEBUG
  assert( ((uintptr_t)i_timeIntegratedDegreesOfFreedom) % ALIGNMENT == 0 );
  assert( ((uintptr_t)data.dofs)              % ALIGNMENT == 0 );
#endif

  kernel::volume volKrnl = m_volumeKernelPrototype;
  volKrnl.Q = data.dofs;
  volKrnl.I = i_timeIntegratedDegreesOfFreedom;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    volKrnl.star(i) = data.localIntegration.starMatrices[i];
  }

  volKrnl.execute();

  kernel::localFlux lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Q = data.dofs;
  lfKrnl.I = i_timeIntegratedDegreesOfFreedom;
  lfKrnl._prefetch.I = i_timeIntegratedDegreesOfFreedom + tensor::I::size();
  lfKrnl._prefetch.Q = data.dofs + tensor::Q::size();
  

  
  for( unsigned int face = 0; face < 4; face++ ) {
    // no element local contribution in the case of dynamic rupture boundary conditions
    if(data.cellInformation.faceTypes[face] != dynamicRupture) {
      lfKrnl.AplusT = data.localIntegration.nApNm1[face];
      lfKrnl.execute(face);
    }
  }
}

void seissol::kernels::Local::flopsIntegral(  enum faceType const i_faceTypes[4],
                                              unsigned int        &o_nonZeroFlops,
                                              unsigned int        &o_hardwareFlops )
{
  o_nonZeroFlops = seissol::kernel::volume::NonZeroFlops;
  o_hardwareFlops = seissol::kernel::volume::HardwareFlops;

  for( unsigned int face = 0; face < 4; ++face ) {
    if( i_faceTypes[face] != dynamicRupture ) {
      o_nonZeroFlops  += seissol::kernel::localFlux::nonZeroFlops(face);
      o_hardwareFlops += seissol::kernel::localFlux::hardwareFlops(face);
    }
  }
}

unsigned seissol::kernels::Local::bytesIntegral()
{
  unsigned reals = 0;

  // star matrices load
  reals += yateto::computeFamilySize<tensor::star>();
  // flux solvers
  reals += 4 * tensor::AplusT::size();

  // DOFs write
  reals += tensor::Q::size();
  
  return reals * sizeof(real);
}
