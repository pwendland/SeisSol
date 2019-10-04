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
  kernels::LocalTmp tmp;

  // switch global data to a device
  m_timeKernel.setGlobalData(&m_DeviceGlobalData);

  seissol::initializers::DeviceVarInfo &manager = layer.getDeviceVarInfo();
  assert(manager.ready() && "a layer has not been analyzed for execution");

  manager.allocateMemory(INTEGRATED_DOFS_ID);
  manager.allocateMemory(DERIVATIVES_ID);
  manager.moveToDevice(m_lts.dofs, DOFS_ID);
  manager.copyDeviceToDevice(DOFS_ID, DERIVATIVES_ID);
  manager.moveToDevice(m_lts.localIntegration, &LocalIntegrationData::starMatrices, STARS_ID);

  // compute derivatives and time integrated dofs
  m_timeKernel.computeAderModified((double)m_timeStepWidthSimulation,
                                   manager,
                                   tmp,
                                   nrOfCells);


  device_copy_from(*(layer.var(m_lts.dofs)),
                   manager.getDevicePointer(INTEGRATED_DOFS_ID),
                   manager.getArraySize(INTEGRATED_DOFS_ID) * sizeof(real));



  // switch global data back to the host
  m_timeKernel.setGlobalData(&m_globalData);

  // free memory from a device
  manager.freeAll();
}


//--------------------------------------------------------------------------------------------------
void computeLocalWithoutAderIntegrationModified() {
  auto& layer = m_ltsTree.child(0).child<Interior>();
  unsigned nrOfCells = layer.getNumberOfCells();
  real** buffers = layer.var(m_lts.buffers);


  kernels::LocalData::Loader loader;
  loader.load(m_lts, layer);
  kernels::LocalTmp tmp;

  // switch global data to a device
  m_localKernel.setGlobalData(&m_DeviceGlobalData);

  seissol::initializers::DeviceVarInfo &manager = layer.getDeviceVarInfo();
  assert(manager.ready() && "a layer has not been analyzed for execution");

  manager.allocateMemory(INTEGRATED_DOFS_ID);
  device_copy_to(manager.getDevicePointer(INTEGRATED_DOFS_ID),
                 *buffers,
                 manager.getArraySize(INTEGRATED_DOFS_ID) * sizeof(real));


  manager.moveToDevice(m_lts.dofs, DOFS_ID);
  manager.moveToDevice(m_lts.localIntegration, &LocalIntegrationData::starMatrices, STARS_ID);
  manager.moveToDevice(m_lts.localIntegration, &LocalIntegrationData::nApNm1, APLUST_ID);


  m_localKernel.computeIntegralModified(buffers,
                                        loader,
                                        manager,
                                        nrOfCells,
                                        tmp);

  // copy updated dofs back to the host

  device_copy_from(*(layer.var(m_lts.dofs)),
                   manager.getDevicePointer(DOFS_ID),
                   manager.getArraySize(DOFS_ID) * sizeof(real));

  // switch global data back to the host
  m_localKernel.setGlobalData(&m_globalData);

  // free memory from a device
  manager.freeAll();
}


//--------------------------------------------------------------------------------------------------
/** Integrates local DOFs for a time cluster as well as volume and local flux integrals
 * NOTE: this function can run only with a GPU instance
 * */
#include "seissol_src/Initializer/LTS.h"
#include "seissol_src/Initializer/tree/DeviceVarInfo.hpp"

void computeLocalIntegrationModified() {
#ifndef GPU
//#error "Cannot run computeLocalIntegrationModified without a GPU instance"
#endif

#pragma message("MESSAGE: 'computeLocalIntegration' procedure is switched to get running on GPU")

  auto& layer = m_ltsTree.child(0).child<Interior>();
  unsigned nrOfCells = layer.getNumberOfCells();
  real** buffers = layer.var(m_lts.buffers);
  real** derivatives = layer.var(m_lts.derivatives);

  kernels::LocalData::Loader loader;
  loader.load(m_lts, layer);
  kernels::LocalTmp tmp;

  // switch global data to a device
  m_timeKernel.setGlobalData(&m_DeviceGlobalData);
  m_localKernel.setGlobalData(&m_DeviceGlobalData);

  seissol::initializers::DeviceVarInfo &manager = layer.getDeviceVarInfo();
  assert(manager.ready() && "a layer has not been analyzed for execution");

  // copy data to a device
  manager.allocateMemory(INTEGRATED_DOFS_ID);
  manager.allocateMemory(DERIVATIVES_ID);
  manager.moveToDevice(m_lts.dofs, DOFS_ID);
  manager.copyDeviceToDevice(DOFS_ID, DERIVATIVES_ID);
  manager.moveToDevice(m_lts.localIntegration, &LocalIntegrationData::starMatrices, STARS_ID);


  // compute derivatives and time integrated dofs
  m_timeKernel.computeAderModified((double)m_timeStepWidthSimulation,
                                   manager,
                                   tmp,
                                   nrOfCells);

  // prepare data for volume and local flux integral
  manager.moveToDevice(m_lts.localIntegration, &LocalIntegrationData::nApNm1, APLUST_ID);

  // copy integrated dofs back to the host
  device_copy_from((*buffers),
                   manager.getDevicePointer(INTEGRATED_DOFS_ID),
                   manager.getArraySize(INTEGRATED_DOFS_ID) * sizeof(real));

  // compute volume and local flux kernels
  m_localKernel.computeIntegralModified(buffers,
                                        loader,
                                        manager,
                                        nrOfCells,
                                        tmp);

  // copy updated dofs back to the host
  device_copy_from(*(layer.var(m_lts.dofs)),
                   manager.getDevicePointer(DOFS_ID),
                   manager.getArraySize(DOFS_ID) * sizeof(real));

  // switch global data back to the host
  m_timeKernel.setGlobalData(&m_globalData);
  m_localKernel.setGlobalData(&m_globalData);

  // free memory from a device
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
