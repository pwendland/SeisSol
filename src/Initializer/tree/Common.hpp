#ifndef INITIALIZER_TREE_COMMON_HPP_
#define INITIALIZER_TREE_COMMON_HPP_

#include <bitset>
#include <limits>
#include <cstring>
#include <Initializer/MemoryAllocator.h>



enum LayerType {
  Ghost = (1 << 0),
  Copy = (1 << 1),
  Interior = (1 << 2),
  NUMBER_OF_LAYERS
};

namespace seissol {
  namespace initializers {
    typedef std::bitset<NUMBER_OF_LAYERS> LayerMask;
    template<typename T> class Variable;
    class Bucket;
    struct MemoryInfo;
  }
}



template<typename T>
struct seissol::initializers::Variable {
  unsigned index;
  LayerMask mask;
  unsigned count;

  Variable() : index(std::numeric_limits<unsigned>::max()), count(1) {}
};

struct seissol::initializers::Bucket {
  unsigned index;

  Bucket() : index(std::numeric_limits<unsigned>::max()) {}
};


struct seissol::initializers::MemoryInfo {
  size_t bytes;
  size_t alignment;
  LayerMask mask;
  memory::Memkind memkind;
};

#endif  // INITIALIZER_TREE_COMMON_HPP_