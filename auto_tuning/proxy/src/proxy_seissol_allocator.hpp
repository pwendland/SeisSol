/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
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
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 **/

#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/DynamicRupture.h>
#include <Initializer/GlobalData.h>
#include <Solver/time_stepping/MiniSeisSol.cpp>
#include <yateto.h>

// DEBUGGING
#include "device_utils.h"
#include <generated_code/init.h>
#include <string>
#include <iostream>


seissol::initializers::LTSTree               m_ltsTree;
seissol::initializers::LTS                   m_lts;
seissol::initializers::LTSTree               m_dynRupTree;
seissol::initializers::DynamicRupture        m_dynRup;

GlobalData m_globalData;
GlobalData m_DeviceGlobalData;

real* m_fakeDerivatives = nullptr;

seissol::kernels::Time      m_timeKernel;
seissol::kernels::Local     m_localKernel;
seissol::kernels::Neighbor  m_neighborKernel;
seissol::kernels::DynamicRupture m_dynRupKernel;


real m_timeStepWidthSimulation = (real)1e-3;

namespace tensor = seissol::tensor;


/** Compares global data allocated on host with data copied on device.
 *
 * NOTE: it is a helper function for debugging. Which checks whether the data
 * between host and device are consistent.
 *
 * @param host Global data allocated on the host.
 * @param device Global data allocated on a device.
 * */
void compareGlobalData(const GlobalData &host, const GlobalData &device);


unsigned int init_data_structures(unsigned int i_cells,
                                  bool enableDynamicRupture,
                                  seissol::memory::ManagedAllocator &m_allocator)
{
  // init RNG
  srand48(i_cells);

  // 1. copy matrices from yateto to aligned memory allocated by seissol
  // 2. provide actual pointers to tensors and matrices back to yateto
  // 3. init integration-LTS-Buffer for all openmp threads
  seissol::initializers::initializeGlobalData(m_globalData, m_allocator, MEMKIND_GLOBAL);

#ifdef GPU
  // Do the same as above but with memory allocated on device
  seissol::initializers::initializeGlobalDataOnDevice(m_DeviceGlobalData,
                                                      m_allocator,
                                                      seissol::memory::Memkind::DeviceGlobalMemory);


  // DEBUGGING: make sure that arrays on host and device are the same
#endif

#ifdef GPU_DEBUGGIN
  compareGlobalData(m_globalData, m_DeviceGlobalData);
#endif

  // TODO: provide data to derivative kernel
  m_timeKernel.setGlobalData(&m_globalData);

  // TODO: provide data to volume and localFlux kernels
  m_localKernel.setGlobalData(&m_globalData);

  // TODO: provide data to localFlux, neighboringFlux and nodalFlux
  m_neighborKernel.setGlobalData(&m_globalData);
  
  m_lts.addTo(m_ltsTree);
  m_ltsTree.setNumberOfTimeClusters(1);
  m_ltsTree.fixate();
  
  seissol::initializers::TimeCluster& cluster = m_ltsTree.child(0);
  cluster.child<Ghost>().setNumberOfCells(0);
  cluster.child<Copy>().setNumberOfCells(0);
  cluster.child<Interior>().setNumberOfCells(i_cells);

#ifdef GPU
  seissol::initializers::prepareDeviceData(m_ltsTree, m_allocator);
#endif

  seissol::initializers::Layer& layer = cluster.child<Interior>();

  // allocate memory space to hold derivatives of Taylor expansion for each element
  layer.setBucketSize(m_lts.buffersDerivatives, sizeof(real) * tensor::I::size() * layer.getNumberOfCells());
  
  m_ltsTree.allocateVariables();
  m_ltsTree.touchVariables();
  m_ltsTree.allocateBuckets();
  
  if (enableDynamicRupture) {
    m_dynRup.addTo(m_dynRupTree);
    m_dynRupTree.setNumberOfTimeClusters(1);
    m_dynRupTree.fixate();
    
    seissol::initializers::TimeCluster& cluster = m_dynRupTree.child(0);
    cluster.child<Ghost>().setNumberOfCells(0);
    cluster.child<Copy>().setNumberOfCells(0);
    cluster.child<Interior>().setNumberOfCells(4*i_cells); /// Every face is a potential dynamic rupture face
  
    m_dynRupTree.allocateVariables();
    m_dynRupTree.touchVariables();
    
    m_fakeDerivatives = (real*) m_allocator.allocateMemory(i_cells * yateto::computeFamilySize<tensor::dQ>() * sizeof(real), PAGESIZE_HEAP, MEMKIND_TIMEDOFS);
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
    for (unsigned cell = 0; cell < i_cells; ++cell) {
      for (unsigned i = 0; i < yateto::computeFamilySize<tensor::dQ>(); i++) {
        m_fakeDerivatives[cell*yateto::computeFamilySize<tensor::dQ>() + i] = (real)drand48();
      }
    }
  }

  /* cell information and integration data*/
  seissol::fakeData(m_lts, layer, (enableDynamicRupture) ? dynamicRupture : regular);

  if (enableDynamicRupture) {
    // From lts tree
    CellDRMapping (*drMapping)[4] = m_ltsTree.var(m_lts.drMapping);

    // From dynamic rupture tree
    seissol::initializers::Layer& interior = m_dynRupTree.child(0).child<Interior>();
    real (*imposedStatePlus)[seissol::tensor::godunovState::size()] = interior.var(m_dynRup.imposedStatePlus);
    real (*fluxSolverPlus)[seissol::tensor::fluxSolver::size()]     = interior.var(m_dynRup.fluxSolverPlus);
    real** timeDerivativePlus = interior.var(m_dynRup.timeDerivativePlus);
    real** timeDerivativeMinus = interior.var(m_dynRup.timeDerivativeMinus);
    DRFaceInformation* faceInformation = interior.var(m_dynRup.faceInformation);
    
    /* init drMapping */
    for (unsigned cell = 0; cell < i_cells; ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        CellDRMapping& drm = drMapping[cell][face];
        unsigned side = (unsigned int)lrand48() % 4;
        unsigned orientation = (unsigned int)lrand48() % 3;
        unsigned drFace = (unsigned int)lrand48() % interior.getNumberOfCells();
        drm.side = side;
        drm.faceRelation = orientation;
        drm.godunov = imposedStatePlus[drFace];
        drm.fluxSolver = fluxSolverPlus[drFace];
      }
    }

    /* init dr godunov state */
    for (unsigned face = 0; face < interior.getNumberOfCells(); ++face) {
      unsigned plusCell = (unsigned int)lrand48() % i_cells;
      unsigned minusCell = (unsigned int)lrand48() % i_cells;
      timeDerivativePlus[face] = &m_fakeDerivatives[plusCell * yateto::computeFamilySize<tensor::dQ>()];
      timeDerivativeMinus[face] = &m_fakeDerivatives[minusCell * yateto::computeFamilySize<tensor::dQ>()];
      
      faceInformation[face].plusSide = (unsigned int)lrand48() % 4;
      faceInformation[face].minusSide = (unsigned int)lrand48() % 4;
      faceInformation[face].faceRelation = (unsigned int)lrand48() % 3;
    }
  }

  m_ltsTree.setDeviceVarInfo(m_lts);  // DEBUGGING
  
  seissol::initializers::DeviceVarInfo &manager = layer.getDeviceVarInfo();
  manager.allocateMemory(INTEGRATED_DOFS_ID);
  manager.allocateMemory(DERIVATIVES_ID);
  manager.allocateMemory(DOFS_ID);
  manager.allocateMemory(STARS_ID);
  manager.allocateMemory(APLUST_ID);
  manager.moveToDevice(m_lts.localIntegration, &LocalIntegrationData::starMatrices, STARS_ID);
  manager.moveToDevice(m_lts.localIntegration, &LocalIntegrationData::nApNm1, APLUST_ID);
  
  return i_cells;
}



void compareGlobalData(const GlobalData &host, const GlobalData &device) {
  for (unsigned i = 0; i < 3; ++i) {
    std::string array_name = "kDivM(" + std::to_string(i)+ ")";
    device_compare_with_host_array(host.stiffnessMatrices(i),
                                   device.stiffnessMatrices(i),
                                   seissol::init::kDivM::size(i),
                                   array_name.c_str());
  }


  for (unsigned i = 0; i < 3; ++i) {
    std::string array_name = "kDivMT(" + std::to_string(i)+ ")";
    device_compare_with_host_array(host.stiffnessMatricesTransposed(i),
                                   device.stiffnessMatricesTransposed(i),
                                   seissol::init::kDivMT::size(i),
                                   array_name.c_str());
  }

  for (unsigned i = 0; i < 4; ++i) {
    std::string array_name = "rDivM(" + std::to_string(i)+ ")";
    device_compare_with_host_array(host.changeOfBasisMatrices(i),
                                   device.changeOfBasisMatrices(i),
                                   seissol::init::rDivM::size(i),
                                   array_name.c_str());
  }

  for (unsigned i = 0; i < 4; ++i) {
    std::string array_name = "rT(" + std::to_string(i)+ ")";
    device_compare_with_host_array(host.neighbourChangeOfBasisMatricesTransposed(i),
                                   device.neighbourChangeOfBasisMatricesTransposed(i),
                                   seissol::init::rT::size(i),
                                   array_name.c_str());
  }

  for (unsigned i = 0; i < 4; ++i) {
    std::string array_name = "fMrT(" + std::to_string(i)+ ")";
    device_compare_with_host_array(host.localChangeOfBasisMatricesTransposed(i),
                                   device.localChangeOfBasisMatricesTransposed(i),
                                   seissol::init::fMrT::size(i),
                                   array_name.c_str());
  }

  for (unsigned i = 0; i < 3; ++i) {
    std::string array_name = "fP(" + std::to_string(i)+ ")";
    device_compare_with_host_array(host.neighbourFluxMatrices(i),
                                   device.neighbourFluxMatrices(i),
                                   seissol::init::fP::size(i),
                                   array_name.c_str());
  }

  device_compare_with_host_array(host.evalAtQPMatrix,
                                 device.evalAtQPMatrix,
                                 seissol::init::evalAtQP::size(),
                                 "evalAtQP");

  device_compare_with_host_array(host.projectQPMatrix,
                                 device.projectQPMatrix,
                                 seissol::init::projectQP::size(),
                                 "projectQP");


}