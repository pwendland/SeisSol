//
// Created by Ravil
//

#ifndef INITIALIZER_TREE_DEVICEVARINFO_HPP_
#define INITIALIZER_TREE_DEVICEVARINFO_HPP_

#include <iostream>  // DEBUGGING
#include <array>
#include "../LTS.h"


#include <type_traits>
enum TensorID {DOF_ID = 0, BUFFER_ID, DERIVATIVES_ID, STARS_ID, APLUST_ID, NUM_IDS};
//enum IndicesNames {DOF, DOF_DERIVATIVES, STARS, NUM_ARRAYS};


namespace seissol {
  namespace initializers {
    //struct LTS;
    class DeviceVarInfo;
  }
}

class seissol::initializers::DeviceVarInfo {
public:
  DeviceVarInfo();
  ~DeviceVarInfo();

  void setData(unsigned numberOfCells, void **vars, LayerType layerType, LTS &tree_structure);



  void collectInfo();

  /** Transfers data from a given structure
   * SOLUTION for offset calculation;
   *    https://en.cppreference.com/w/cpp/language/pointer
   *    https://gist.github.com/graphitemaster/494f21190bb2c63c5516
   * */
  template <typename H, typename T>
  void moveToDevice(seissol::initializers::Variable<H> const& handle, T H::*member, TensorID ID) {

    // compute offset and data_size to transfer
    H object{};  // TODO: rework to make a constant expression out of it
    size_t offset = size_t(&(object.*member)) - size_t(&object);

    real* d_ptr = (real*)device_malloc(m_ArraySizes[ID] * sizeof(real));
    device_copy_2D_to((void*)d_ptr,
                      sizeof(T),
                      (void*)(m_vars[handle.index] + offset),
                      sizeof(H),
                      sizeof(T),
                      m_numberOfCells);

    m_DevicePointers[ID] = d_ptr;
  }

  template <typename H>
  void moveToDevice(seissol::initializers::Variable<H> const& handle, TensorID ID) {
    real* d_ptr = (real*)device_malloc(m_ArraySizes[ID] * sizeof(real));
    device_copy_to((void*)d_ptr, (void*)m_vars[handle.index], m_ArraySizes[ID] * sizeof(real));
    m_DevicePointers[ID] = d_ptr;
  }

  real* getDevicePointer(TensorID ID);
  void free(TensorID ID);
  void freeAll();

  unsigned long getArraySize(TensorID ID) {return m_ArraySizes[ID]; }
  unsigned long getIndex(TensorID ID) {return m_Indices[ID]; }
  bool ready() {return m_is_ready;}

private:
  void computeArraySizes();
  void computeIndices();

  unsigned m_numberOfCells;
  void **m_vars;
  LayerType m_layerType{};
  LTS m_tree_structure{};


  std::array<unsigned long, NUM_IDS> m_ArraySizes{};
  std::array<unsigned long, NUM_IDS> m_Indices{};
  std::array<real*, NUM_IDS> m_DevicePointers{};
  bool m_is_ready;
};


#endif //INITIALIZER_TREE__DEVICEVARINFO_HPP_
