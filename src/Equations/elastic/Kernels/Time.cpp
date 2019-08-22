/******************************************************************************
** Copyright (c) 2014-2015, Intel Corporation                                **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
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
 * Time kernel of SeisSol.
 **/

#include "Kernels/TimeBase.h"
#include "Kernels/Time.h"

#ifndef NDEBUG
#pragma message "compiling time kernel with assertions"
extern long long libxsmm_num_total_flops;
#endif

#include <Kernels/denseMatrixOps.hpp>

#include <cstring>
#include <cassert>
#include <stdint.h>
#include <iostream>

#include <omp.h>

#include <yateto.h>

//DEBUGGING
#include "device_utils.h"
#include <algorithm>

seissol::kernels::TimeBase::TimeBase() {
  m_derivativesOffsets[0] = 0;
  for (int order = 0; order < CONVERGENCE_ORDER; ++order) {
    if (order > 0) {
      m_derivativesOffsets[order] = tensor::dQ::size(order-1) + m_derivativesOffsets[order-1];
    }
  }
}

void seissol::kernels::Time::setGlobalData(GlobalData const* global) {
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(0)) % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(1)) % ALIGNMENT == 0 );
  assert( ((uintptr_t)global->stiffnessMatricesTransposed(2)) % ALIGNMENT == 0 );

  m_krnlPrototype.kDivMT = global->stiffnessMatricesTransposed;
}

#define Q_VALUE 1.0
void seissol::kernels::Time::computeAderModified(double i_timeStepWidth,
                                                 LocalData::Loader& loader,
                                                 LocalTmp& tmp,
                                                 const unsigned num_cells,
                                                 real *o_timeIntegrated[tensor::I::size()],
                                                 real **o_timeDerivatives)
{
    // TODO: separate Teyalor expansion kernl into 2 parts, namely:
    // TODO:    1. A sublayer that have to hold derivatives
    // TODO:    2. A sublayer that doesn't have to store derivatives

    // TODO: allocate global data on gpu in advance. During allocation and initialization of the global data
    /*
     * assert alignments.
     */
    // TODO: enable
    /*
    assert( ((uintptr_t)data.dofs)              % ALIGNMENT == 0 );
    assert( ((uintptr_t)o_timeIntegrated )      % ALIGNMENT == 0 );
    assert( ((uintptr_t)o_timeDerivatives)      % ALIGNMENT == 0 || o_timeDerivatives == NULL );
    */

    /*
     * compute ADER scheme.
     */
    // temporary result

    // allocate a temporary storage for derivatives
    // in case of the user did not provide space to store derivatives

    kernel::derivative krnl = m_krnlPrototype;
    kernel::derivativeTaylorExpansion intKrnl;

    // allocate and init derivative buffers for all elements i.e. dQ(i)
    // Stitch Seissol and Yateto together
    real *d_temporaryBuffer[yateto::numFamilyMembers<tensor::dQ>()];
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
        device_malloc((void**)&d_temporaryBuffer[i], num_cells * tensor::dQ::Size[i] * sizeof(real));
    }


    // stream dofs to the derivative zero
    for (unsigned cell_indx = 0; cell_indx < num_cells; ++cell_indx) {
        auto data = loader.entry(cell_indx);

        device_copy_to((void*)(d_temporaryBuffer[0] + cell_indx * tensor::dQ::Size[0]),
                       (void*)data.dofs,
                       tensor::dQ::Size[0] * sizeof(real));
    }


    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
        krnl.dQ(i) = d_temporaryBuffer[i];
        intKrnl.dQ(i) = d_temporaryBuffer[i];
    }

    // allocate and init integrated unknowns tensors for all elements i.e. I
    // Stitch Seissol and Yateto together
    real *d_o_timeIntegrated;
    device_malloc((void**)&d_o_timeIntegrated, num_cells * tensor::I::Size * sizeof(real));
    device_copy_to((void*)d_o_timeIntegrated, (void*)o_timeIntegrated[0], num_cells * tensor::I::Size * sizeof(real));
    intKrnl.I = d_o_timeIntegrated;

    // allocate and init star tensors for all elements i.e. start(i)
    // Stitch Seissol and Yateto together
    real *d_stars[yateto::numFamilyMembers<tensor::star>()];
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
        device_malloc((void**)&d_stars[i], num_cells * tensor::star::Size[i] * sizeof(real));
    }

    // initialization of star matrices on GPU
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
        for (unsigned cell_indx = 0; cell_indx < num_cells; ++cell_indx) {
            auto data = loader.entry(cell_indx);

            device_copy_to((void*)(d_stars[i] + cell_indx * tensor::star::Size[i]),
                           (void*)data.localIntegration.starMatrices[i],
                           tensor::star::Size[i] * sizeof(real));
        }
    }

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
        krnl.star(i) = d_stars[i];
    }


    // allocate and init stiffness tensors of the reference element i.e. d_kDivMT(i)
    // Stitch Seissol and Yateto together
    real *d_kDivMT[yateto::numFamilyMembers<tensor::kDivMT>()];
    for (int i = 0; i < yateto::numFamilyMembers<tensor::kDivMT>(); ++i) {
        device_malloc((void**)&d_kDivMT[i], tensor::kDivMT::Size[i] * sizeof(real));


        device_copy_to((void*)(d_kDivMT[i]),
                       (void*)m_krnlPrototype.kDivMT(i),
                       tensor::kDivMT::Size[i] * sizeof(real));
    }

    for (int i = 0; i < yateto::numFamilyMembers<tensor::kDivMT>(); ++i) {
        krnl.kDivMT(i) = d_kDivMT[i];
    }


    // stitch a scalar which is used during Taylor expansion
    tensor::num_elements_in_cluster = num_cells;
    intKrnl.power = i_timeStepWidth;
    intKrnl.execute0();

    for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {
        krnl.execute(der);

        // update scalar for this derivative
        intKrnl.power *= i_timeStepWidth / real(der + 1);
        intKrnl.execute(der);
    }

    // TODO: check the assumption that elements within the linear space
    // NOTE: debugging showed that it worked.
    // copy results back to CPU
    device_copy_from((void*)o_timeIntegrated[0], (void*)d_o_timeIntegrated, num_cells * tensor::I::Size * sizeof(real));

    // free GPU memory
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
        device_free((void*)d_temporaryBuffer[i]);
    }

    device_free((void*)d_o_timeIntegrated);
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
        device_free((void*)d_stars[i]);
    }


    for (int i = 0; i < yateto::numFamilyMembers<tensor::kDivM>(); ++i) {
        device_free((void*)d_kDivMT[i]);
    }
}


void seissol::kernels::Time::computeAder(double i_timeStepWidth,
                                         LocalData& data,
                                         LocalTmp& tmp,
                                         real o_timeIntegrated[tensor::I::size()],
                                         real* o_timeDerivatives)
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)data.dofs)              % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated )      % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeDerivatives)      % ALIGNMENT == 0 || o_timeDerivatives == NULL );

  /*
   * compute ADER scheme.
   */
  // temporary result

  // allocate a temporary storage for derivatives
  // in case of the user did not provide space to store derivatives

  real temporaryBuffer[yateto::computeFamilySize<tensor::dQ>()] __attribute__((aligned(PAGESIZE_STACK)));

  // decide whether the store derivatives in the user's provided containers or locally
  // NOTE: we have to store derivatives because of the uniform interface of source code generator
  // TODO: ask whether it is done to avoid an extra layer??? I mean some interior cells must keep
  // TODO: their derivatives in order to update (or be updated by) other time clusters
  real* derivativesBuffer = (o_timeDerivatives != nullptr) ? o_timeDerivatives : temporaryBuffer;

  // prepare and init 'derivative' kernel generated by yateto
  kernel::derivative krnl = m_krnlPrototype;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
    krnl.star(i) = data.localIntegration.starMatrices[i];
  }

  krnl.dQ(0) = const_cast<real*>(data.dofs);  // <-- ask
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    krnl.dQ(i) = derivativesBuffer + m_derivativesOffsets[i];
  }



  // prepare and init 'derivativeTaylorExpansion' kernel generated by yateto
  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;
  intKrnl.dQ(0) = data.dofs;  // <-- ask
  for (unsigned i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = derivativesBuffer + m_derivativesOffsets[i];
  }


  // powers in the taylor-series expansion
  intKrnl.power = i_timeStepWidth;
  intKrnl.execute0();

  // stream out first derivative (order 0)
  if (o_timeDerivatives != nullptr) {
    // NOTE: both intKrnl.dQ(0) and krnl.dQ(0) hold pointers to data.dofs
    //       derivatives from 1 to <order> will be written to o_timeDerivatives
    //       but not the 0th one. That is the reason why we have to copy dofs in the line below
    streamstore(tensor::dQ::size(0), data.dofs, o_timeDerivatives);  // <-- ask
  }


  for (unsigned der = 1; der < CONVERGENCE_ORDER; ++der) {
    krnl.execute(der);

    // update scalar for this derivative
    intKrnl.power *= i_timeStepWidth / real(der+1);    
    intKrnl.execute(der);
  }
}

void seissol::kernels::Time::flopsAder( unsigned int        &o_nonZeroFlops,
                                        unsigned int        &o_hardwareFlops ) {
  // reset flops
  o_nonZeroFlops = 0; o_hardwareFlops =0;

  // initialization
  o_nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(0);
  o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(0);

  // interate over derivatives
  for( unsigned l_derivative = 1; l_derivative < CONVERGENCE_ORDER; l_derivative++ ) {
    o_nonZeroFlops  += kernel::derivative::nonZeroFlops(l_derivative);
    o_hardwareFlops += kernel::derivative::hardwareFlops(l_derivative);

    // update of time integrated DOFs
    o_nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(l_derivative);
    o_hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(l_derivative);
  }

}

unsigned seissol::kernels::Time::bytesAder()
{
  unsigned reals = 0;
  
  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q::size() + 2 * tensor::I::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>();
           
  /// \todo incorporate derivatives

  return reals * sizeof(real);
}

void seissol::kernels::Time::computeIntegral( double                            i_expansionPoint,
                                              double                            i_integrationStart,
                                              double                            i_integrationEnd,
                                              const real*                       i_timeDerivatives,
                                              real                              o_timeIntegrated[tensor::I::size()] )
{
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)i_timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)o_timeIntegrated)   % ALIGNMENT == 0 );

  // assert that this is a forwared integration in time
  assert( i_integrationStart + (real) 1.E-10 > i_expansionPoint   );
  assert( i_integrationEnd                   > i_integrationStart );

  /*
   * compute time integral.
   */
  // compute lengths of integration intervals
  real l_deltaTLower = i_integrationStart - i_expansionPoint;
  real l_deltaTUpper = i_integrationEnd   - i_expansionPoint;

  // initialization of scalars in the taylor series expansion (0th term)
  real l_firstTerm  = (real) 1;
  real l_secondTerm = (real) 1;
  real l_factorial  = (real) 1;
  
  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = o_timeIntegrated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = i_timeDerivatives + m_derivativesOffsets[i];
  }
 
  // iterate over time derivatives
  for(int der = 0; der < CONVERGENCE_ORDER; ++der ) {
    l_firstTerm  *= l_deltaTUpper;
    l_secondTerm *= l_deltaTLower;
    l_factorial  *= (real)(der+1);

    intKrnl.power  = l_firstTerm - l_secondTerm;
    intKrnl.power /= l_factorial;

    intKrnl.execute(der);
  }
}

void seissol::kernels::Time::computeTaylorExpansion( real         time,
                                                     real         expansionPoint,
                                                     real const*  timeDerivatives,
                                                     real         timeEvaluated[tensor::Q::size()] ) {
  /*
   * assert alignments.
   */
  assert( ((uintptr_t)timeDerivatives)  % ALIGNMENT == 0 );
  assert( ((uintptr_t)timeEvaluated)    % ALIGNMENT == 0 );

  // assert that this is a forward evaluation in time
  assert( time >= expansionPoint );

  real deltaT = time - expansionPoint;

  static_assert(tensor::I::size() == tensor::Q::size(), "Sizes of tensors I and Q must match");

  kernel::derivativeTaylorExpansion intKrnl;
  intKrnl.I = timeEvaluated;
  for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
    intKrnl.dQ(i) = timeDerivatives + m_derivativesOffsets[i];
  }
  intKrnl.power = 1.0;
 
  // iterate over time derivatives
  for(int derivative = 0; derivative < CONVERGENCE_ORDER; ++derivative) {
    intKrnl.execute(derivative);
    intKrnl.power *= deltaT / real(derivative+1);
  }
}

void seissol::kernels::Time::flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops) {
  // reset flops
  nonZeroFlops = 0; hardwareFlops = 0;

  // interate over derivatives
  for (unsigned der = 0; der < CONVERGENCE_ORDER; ++der) {
    nonZeroFlops  += kernel::derivativeTaylorExpansion::nonZeroFlops(der);
    hardwareFlops += kernel::derivativeTaylorExpansion::hardwareFlops(der);
  }
}
