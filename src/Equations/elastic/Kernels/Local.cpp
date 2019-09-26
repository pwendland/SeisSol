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
  real *Q = (real*)device_malloc(num_cells * tensor::Q::Size * sizeof(real));  // TODO: remove allocation of memory
  real *I = (real*)device_malloc(num_cells * tensor::I::Size * sizeof(real));  // TODO: remove allocation of memory

  device_copy_to((void*)I,
                 (void*)i_timeIntegratedDegreesOfFreedom[0],
                 num_cells * tensor::I::Size * sizeof(real));




  // allocate memory for star matrices for all elements
  size_t stars_size =  yateto::numFamilyMembers<tensor::star>() * tensor::star::Size[0];
  real *d_stars = (real*)device_malloc(stars_size * num_cells * sizeof(real));  // TODO: remove allocation of memory


  // copy star matrices to a device (extracting matrices from and array of structures)
  device_copy_2D_to((void*)d_stars,
                    stars_size * sizeof(real),
                    (void*)data.localIntegration.starMatrices,
                    sizeof(LocalIntegrationData),
                    stars_size * sizeof(real),
                    num_cells);


  volKrnl.Q = Q;
  volKrnl.I = I;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    //seissol::tensor::star::jump_to_next[i] = stars_size;  // TODO: this line must be enable i.e a programmer must control this
    volKrnl.star(i) = d_stars + i * tensor::star::Size[0];
  }

  tensor::num_elements_in_cluster = num_cells;
  volKrnl.execute();

  //-------------------------------------------------------------------------------------------------------------------

  size_t AplusT_size = 4 * tensor::AplusT::Size;
  real *d_AplusT = (real*)device_malloc(AplusT_size * num_cells * sizeof(real));  // TODO: remove allocation of memory


  // copy star matrices to a device (extracting matrices from and array of structures)
  device_copy_2D_to((void*)d_AplusT,
                    AplusT_size * sizeof(real),
                    (void*)data.localIntegration.nApNm1[0],
                    sizeof(LocalIntegrationData),
                    AplusT_size * sizeof(real),
                    num_cells);


  // retrieve indices from if-else statement
  std::vector<unsigned> bins[4];
  for(unsigned int face = 0; face < 4; face++) {
    unsigned counter = 0;
    for(unsigned int l_cell = 0; l_cell < num_cells; l_cell++) {
      if(data.cellInformation.faceTypes[face] != dynamicRupture) {
        bins[face].push_back(counter++);
      }
    }
  }


  kernel::localFlux lfKrnl = m_localFluxKernelPrototype;
  lfKrnl.Q = Q;
  lfKrnl.I = I;
  const unsigned num_faces = 4;
  unsigned *QI_strides[4];
  unsigned *AplusT_strides[4];
  for (unsigned face = 0; face < num_faces; ++face) {
    if (!bins[face].empty()) {
      // allocate indices in a device
      QI_strides[face] = (unsigned *) device_malloc(
          bins[face].size() * sizeof(unsigned));    // TODO: remove allocation of memory
      device_copy_to(QI_strides[face], bins[face].data(), bins[face].size() * sizeof(unsigned));
      device_scale_array(tensor::Q::Size, bins[face].size(), QI_strides[face]);

      AplusT_strides[face] = (unsigned *) device_malloc(
          bins[face].size() * sizeof(unsigned));    // TODO: remove allocation of memory
      device_copy_to(AplusT_strides[face], bins[face].data(), bins[face].size() * sizeof(unsigned));
      device_scale_array(AplusT_size, bins[face].size(), AplusT_strides[face]);

      // attach strides to the kernel
      tensor::Q::stride_ptr = QI_strides[face];
      tensor::I::stride_ptr = QI_strides[face];
      tensor::AplusT::stride_ptr = AplusT_strides[face];
      tensor::fMrT::stride_ptr[face] = seissol::tensor::device_zeros;
      tensor::rDivM::stride_ptr[face] = seissol::tensor::device_zeros;

      // attach a flux solver to the kernel
      lfKrnl.AplusT = d_AplusT + (face * tensor::AplusT::Size);

      // specify the number of faces to compute
      tensor::num_elements_in_cluster = bins[face].size();

      // execute the kernel for a particular face
      lfKrnl.execute(face);
    }
  }

  // copy dofs back to the host
  device_copy_from((void*)data.dofs,
                   (void*)Q,
                   num_cells * tensor::Q::Size * sizeof(real));

  // release resources
  device_free(Q);  // TODO: remove de-allocation of memory
  device_free(I);  // TODO: remove de-allocation of memory
  device_free(d_stars);  // TODO: remove de-allocation of memory
  device_free(d_AplusT);  // TODO: remove de-allocation of memory

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
