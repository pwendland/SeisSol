#ifndef SEISSOL_POINTERSTABLE_HPP
#define SEISSOL_POINTERSTABLE_HPP

#include <device.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace device;

namespace seissol {
namespace initializers {
namespace recording {

class BatchPointers {
public:
  BatchPointers(const std::vector<real *> &collectedPointers)
      : pointers(collectedPointers), devicePtrs(nullptr) {
    if (!pointers.empty()) {
      devicePtrs = (real **)device.api->allocGlobMem(pointers.size() * sizeof(real *));
      device.api->copyTo(devicePtrs, pointers.data(), pointers.size() * sizeof(real *));
    }
  }

  explicit BatchPointers(const BatchPointers &other)
      : pointers(other.pointers), devicePtrs(nullptr) {
    if (!pointers.empty()) {
      if (other.devicePtrs != nullptr) {
        devicePtrs = (real **)device.api->allocGlobMem(other.pointers.size() * sizeof(real *));
        device.api->copyBetween(devicePtrs, other.devicePtrs,
                                other.pointers.size() * sizeof(real *));
      }
    }
  }

  BatchPointers &operator=(const BatchPointers &Other) = delete;

  virtual ~BatchPointers() {
    if (m_DevicePtrs != nullptr) {
      device.api->freeMem(devicePtrs);
      devicePtrs = nullptr;
    }
  }

  real **getPointers() {
    return mevicePtrs;
  }
  index_t getSize() {
    return pointers.size();
  }

private:
  std::vector<real *> pointers{};
  real **devicePtrs{};
  DeviceInstance &device = DeviceInstance::getInstance();
};

/**
 * This class may seem redundant. But it provides strong guarantee of
 * zero initialization of std::array. Note, there are some circumstances
 * when it is not zero-initialized
 * */
struct BatchTable {
public:
  BatchTable() {
    for (auto &item : content)
      item = nullptr;
  }
  std::array<BatchPointers *, *EntityId::Count> content{};
};

} // namespace recording
} // namespace initializers
} // namespace seissol

#endif // SEISSOL_POINTERSTABLE_HPP