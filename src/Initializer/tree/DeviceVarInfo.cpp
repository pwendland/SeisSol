#include <yateto.h>

#include "DeviceVarInfo.hpp"
#include "generated_code/tensor.h"
#include "generated_code/init.h"
#include "../typedefs.hpp"

#include <iostream>  // DEBUGGING

using namespace seissol::initializers;

DeviceVarInfo::DeviceVarInfo() : m_numberOfCells(0),
                                 m_vars(nullptr),
                                 m_is_ready(false) {
  m_DevicePointers.fill(nullptr);
}

DeviceVarInfo::~DeviceVarInfo() {}

real* DeviceVarInfo::getDevicePointer(TensorID ID) { return m_DevicePointers[ID]; }

inline
void DeviceVarInfo::free(TensorID ID) {
  if (m_DevicePointers[ID] != nullptr) {
    device_free(m_DevicePointers[ID]);
    m_DevicePointers[ID] = nullptr;
  }
}


void DeviceVarInfo::freeAll() {
  for (unsigned id = 0; id < NUM_IDS; ++id) {
    this->free((TensorID)id);
  }
}


void DeviceVarInfo::setData(unsigned numberOfCells, void **vars, LayerType layerType, LTS &tree_structure) {
  m_numberOfCells = numberOfCells;
  m_vars = vars;
  m_layerType = layerType;
  m_tree_structure = tree_structure;
}


void DeviceVarInfo::collectInfo() {
  computeArraySizes();
  computeIndices();

  m_is_ready = true;
}



void DeviceVarInfo::computeArraySizes() {
  m_ArraySizes[DOFS_ID] = m_numberOfCells * tensor::Q::size();
  m_ArraySizes[INTEGRATED_DOFS_ID] = m_numberOfCells * tensor::I::size();
  m_ArraySizes[DERIVATIVES_ID] = m_numberOfCells * yateto::computeFamilySize<init::dQ>();
  m_ArraySizes[STARS_ID] = m_numberOfCells * yateto::computeFamilySize<init::star>();
  m_ArraySizes[APLUST_ID] = m_numberOfCells * 4 * tensor::AplusT::size();
}




void DeviceVarInfo::computeIndices() {
  CellLocalInformation *l_cellInfo = (CellLocalInformation*)m_vars[m_tree_structure.cellInformation.index];

  // check whether a face of a cell belong to dynamic rupture
  for (unsigned face = 0; face < 4; ++face) {
    unsigned counter = 0;
    for (unsigned cell = 0; cell < m_numberOfCells; ++cell) {
      if (l_cellInfo[cell].faceTypes[face] == regular) {
        m_FaceElementBins[LOCAL_FLUX][face].push_back(counter++);
      }
    }
  }
}

