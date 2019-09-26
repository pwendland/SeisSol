/*
  created by Ravil
*/

#include <generated_code/tensor.h>

namespace tensor = seissol::tensor;
namespace kernels = seissol::kernels;


//--------------------------------------------------------------------------------------------------
void computeAderIntegrationModified() {

  auto& layer = m_ltsTree.child(0).child<Interior>();
  unsigned nrOfCells = layer.getNumberOfCells();
  real** buffers = layer.var(m_lts.buffers);
  real** derivatives = layer.var(m_lts.derivatives);

  kernels::LocalData::Loader loader;
  loader.load(m_lts, layer);

#ifdef _OPENMP
#pragma omp parallel
  {
    kernels::LocalTmp tmp;
#pragma omp for schedule(static)
#endif
    for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
      auto data = loader.entry(l_cell);
      m_timeKernel.computeAder(m_timeStepWidthSimulation,
                               data,
                               tmp,
                               buffers[l_cell],
                               derivatives[l_cell] );
    }
#ifdef _OPENMP
  }
#endif
}


//--------------------------------------------------------------------------------------------------
void computeLocalWithoutAderIntegrationModified() {
  auto&                 layer           = m_ltsTree.child(0).child<Interior>();
  unsigned              nrOfCells       = layer.getNumberOfCells();
  real**                buffers                       = layer.var(m_lts.buffers);

  kernels::LocalData::Loader loader;
  loader.load(m_lts, layer);

#ifdef _OPENMP
#pragma omp parallel
  {
    kernels::LocalTmp tmp;
#pragma omp for schedule(static)
#endif
    for( unsigned int l_cell = 0; l_cell < nrOfCells; l_cell++ ) {
      auto data = loader.entry(l_cell);
      m_localKernel.computeIntegral(  buffers[l_cell],
                                      data,
                                      tmp );
    }
#ifdef _OPENMP
  }
#endif
}


//--------------------------------------------------------------------------------------------------
/** Integrates local DOFs for a time cluster as well as volume and local flux integrals
 * NOTE: this function can run only with a GPU instance
 * */
#include "seissol_src/Initializer/LTS.h"
#include "seissol_src/Initializer/tree/DeviceVarInfo.hpp"

void computeLocalIntegrationModified() {
#ifndef GPU
#error "Cannot run computeLocalIntegrationModified without a GPU instance"
#endif

#pragma message("MESSAGE: 'computeLocalIntegration' procedure is switched to get running on GPU")

  auto& layer = m_ltsTree.child(0).child<Interior>();
  unsigned nrOfCells = layer.getNumberOfCells();
  real** buffers = layer.var(m_lts.buffers);
  real** derivatives = layer.var(m_lts.derivatives);

  kernels::LocalData::Loader loader;
  loader.load(m_lts, layer);
  kernels::LocalTmp tmp;

  m_timeKernel.setGlobalData(&m_DeviceGlobalData); // switch to the global data allocated on the device
  m_localKernel.setGlobalData(&m_DeviceGlobalData); // switch to the global data allocated on the device

  // TODO: analyse layer
  seissol::initializers::DeviceVarInfo manager = layer.getDeviceVarInfo();
  //container.moveToDevice(m_lts.localIntegration);

  manager.moveToDevice(m_lts.localIntegration, &LocalIntegrationData::starMatrices, STARS_ID);
  manager.moveToDevice(m_lts.dofs, DOF_ID);


  // TODO: enable assert statement
  //assert(container.ready());




  // TODO: move data on GPU

  m_timeKernel.computeAderModified((double)m_timeStepWidthSimulation,
                                   loader,
                                   manager,
                                   tmp,
                                   nrOfCells,
                                   buffers,
                                   derivatives);



  m_localKernel.computeIntegralModified(buffers,
                                        loader,
                                        manager,
                                        nrOfCells,
                                        tmp);

  m_timeKernel.setGlobalData(&m_globalData); // switch to back to the data allocated on the host
  m_localKernel.setGlobalData(&m_globalData); // switch to back to the data allocated on the host

  manager.freeAll();
}


//--------------------------------------------------------------------------------------------------
/** Mimics integration of neighbor's fluxes
 * NOTE: this function has the same functionality as computeNeighboringIntegration but performs integration
 *       face-wise. The function mimics unboxing element and face IDs from some data structure
 * */
void computeNeighboringIntegrationModified() {
  auto&                     layer                           = m_ltsTree.child(0).child<Interior>();
  unsigned                  nrOfCells                       = layer.getNumberOfCells();
  real*                     (*faceNeighbors)[4]             = layer.var(m_lts.faceNeighbors);
  CellDRMapping             (*drMapping)[4]                 = layer.var(m_lts.drMapping);
  CellLocalInformation*       cellInformation               = layer.var(m_lts.cellInformation);

  kernels::NeighborData::Loader loader;
  loader.load(m_lts, layer);

  real *l_timeIntegrated[4];
#ifdef ENABLE_MATRIX_PREFETCH
  real *l_faceNeighbors_prefetch[4];
#endif

#ifdef _OPENMP
#  ifdef ENABLE_MATRIX_PREFETCH
#pragma omp parallel private(l_timeIntegrated, l_faceNeighbors_prefetch)
#  else
#pragma omp parallel private(l_timeIntegrated)
#  endif
  {
#pragma omp for schedule(static)
#endif
    for(int l_cell = 0; l_cell < nrOfCells; l_cell++) {

      for (unsigned l_face_idx = 0; l_face_idx < 4; ++l_face_idx) {

        std::pair<unsigned int, unsigned int> element_face(l_cell, l_face_idx);
        unsigned int elementIdx = element_face.first;  // mimics unboxing from some data structure
        unsigned int faceIdx = element_face.second;  // mimics unboxing from some data structure

        auto data = loader.entry(elementIdx);

#ifdef _OPENMP
        real *tmp_storage = &(m_globalData.integrationBufferLTS[omp_get_thread_num() * tensor::I::size()]);
#else
        real *tmp_storage = m_globalData.integrationBufferLTS;
#endif

        // prepare DOFs of a neighbor for the subsequent face integral
        seissol::kernels::TimeCommon::computeIntegralsFacewise(m_timeKernel,
                                                               faceIdx,
                                                               cellInformation[elementIdx].ltsSetup,
                                                               cellInformation[elementIdx].faceTypes[faceIdx],
                                                               0.0,
                                                               (double) m_timeStepWidthSimulation,
                                                               faceNeighbors[elementIdx][faceIdx],
                                                               tmp_storage,
                                                               l_timeIntegrated[faceIdx]);

#ifdef ENABLE_MATRIX_PREFETCH
#pragma message("the current prefetch structure (flux matrices and tDOFs is tuned for higher order and shouldn't be harmful for lower orders")
        // prepare matrices for prefetch from the main memory
        if (faceIdx != 3) {
          l_faceNeighbors_prefetch[faceIdx] = (cellInformation[elementIdx].faceTypes[faceIdx + 1] != dynamicRupture)
                                               ? faceNeighbors[elementIdx][faceIdx + 1]
                                               : drMapping[elementIdx][faceIdx + 1].godunov;
        }
        else {
          if (elementIdx < (nrOfCells - 1)) {
            l_faceNeighbors_prefetch[faceIdx] = (cellInformation[elementIdx + 1].faceTypes[0] != dynamicRupture)
                                                 ? faceNeighbors[elementIdx + 1][0]
                                                 : drMapping[elementIdx + 1][0].godunov;
          }
          else {
            l_faceNeighbors_prefetch[faceIdx] = faceNeighbors[elementIdx][faceIdx];
          }
        }
#endif

        // compute the surface integral using neighbor's integrated DOFs
        m_neighborKernel.computeNeighborsIntegralFacewise(faceIdx,
                                                          data,
                                                          drMapping[elementIdx][faceIdx],
#ifdef ENABLE_MATRIX_PREFETCH
                                                          l_timeIntegrated[faceIdx], l_faceNeighbors_prefetch[faceIdx]
#else
                                                          l_timeIntegrated[faceIdx]
#endif
        );
      }
    }

#ifdef _OPENMP
  }
#endif
}
