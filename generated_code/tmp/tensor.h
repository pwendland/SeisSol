#ifndef SEISSOL_TENSOR_H_
#define SEISSOL_TENSOR_H_
namespace seissol {
  namespace tensor {
    struct AminusT {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct AplusT {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct I {
      constexpr static unsigned const Shape[2] = {56, 9};
      constexpr static unsigned const Size = 504;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct Q {
      constexpr static unsigned const Shape[2] = {56, 9};
      constexpr static unsigned const Size = 504;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QAtPoint {
      constexpr static unsigned const Shape[1] = {9};
      constexpr static unsigned const Size = 9;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QFortran {
      constexpr static unsigned const Shape[2] = {56, 9};
      constexpr static unsigned const Size = 504;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QStress {
      constexpr static unsigned const Shape[2] = {56, 6};
      constexpr static unsigned const Size = 336;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QStressNodal {
      constexpr static unsigned const Shape[2] = {56, 6};
      constexpr static unsigned const Size = 336;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QgodLocal {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QgodNeighbor {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct T {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct Tinv {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct V3mTo2n {
      constexpr static unsigned const Shape[][2] = {{49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}};
      constexpr static unsigned const Size[] = {3136, 3136, 3136, 3136, 3136, 3136, 3136, 3136, 3136, 3136, 3136, 3136, 3136, 3136, 3136, 3136};
      constexpr static unsigned index(unsigned i0, unsigned i1) {
        return 1*i0 + 4*i1;
      }
      constexpr static unsigned size(unsigned i0, unsigned i1) {
        return Size[index(i0, i1)];
      }
      template<typename T>
      struct Container {
        T data[16];
        Container() : data{} {}
        inline T& operator()(unsigned i0, unsigned i1) {
          return data[index(i0, i1)];
        }
        inline T const& operator()(unsigned i0, unsigned i1) const {
          return data[index(i0, i1)];
        }
      };
    };
    struct V3mTo2nTWDivM {
      constexpr static unsigned const Shape[][2] = {{56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}};
      constexpr static unsigned const Size[] = {2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744};
      constexpr static unsigned index(unsigned i0, unsigned i1) {
        return 1*i0 + 4*i1;
      }
      constexpr static unsigned size(unsigned i0, unsigned i1) {
        return Size[index(i0, i1)];
      }
      template<typename T>
      struct Container {
        T data[16];
        Container() : data{} {}
        inline T& operator()(unsigned i0, unsigned i1) {
          return data[index(i0, i1)];
        }
        inline T const& operator()(unsigned i0, unsigned i1) const {
          return data[index(i0, i1)];
        }
      };
    };
    struct basisFunctions {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct dQ {
      constexpr static unsigned const Shape[][2] = {{56, 9}, {56, 9}, {56, 9}, {56, 9}, {56, 9}, {56, 9}};
      constexpr static unsigned const Size[] = {504, 360, 216, 144, 72, 72};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[6];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct displacement {
      constexpr static unsigned const Shape[2] = {56, 3};
      constexpr static unsigned const Size = 168;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct dofsQP {
      constexpr static unsigned const Shape[2] = {343, 9};
      constexpr static unsigned const Size = 3096;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct evalAtQP {
      constexpr static unsigned const Shape[2] = {343, 56};
      constexpr static unsigned const Size = 19264;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct fMrT {
      constexpr static unsigned const Shape[][2] = {{21, 56}, {21, 56}, {21, 56}, {21, 56}};
      constexpr static unsigned const Size[] = {1344, 1344, 1344, 1344};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct fP {
      constexpr static unsigned const Shape[][2] = {{21, 21}, {21, 21}, {21, 21}};
      constexpr static unsigned const Size[] = {504, 504, 504};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[3];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct fluxSolver {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct godunovMatrix {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct godunovState {
      constexpr static unsigned const Shape[2] = {49, 9};
      constexpr static unsigned const Size = 504;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct iniCond {
      constexpr static unsigned const Shape[2] = {343, 9};
      constexpr static unsigned const Size = 3096;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct initialLoading {
      constexpr static unsigned const Shape[1] = {6};
      constexpr static unsigned const Size = 6;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct kDivM {
      constexpr static unsigned const Shape[][2] = {{56, 56}, {56, 56}, {56, 56}};
      constexpr static unsigned const Size[] = {1960, 1960, 1960};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[3];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct kDivMT {
      constexpr static unsigned const Shape[][2] = {{56, 56}, {56, 56}, {56, 56}};
      constexpr static unsigned const Size[] = {2120, 2160, 2200};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[3];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct mInvJInvPhisAtSources {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct meanStress {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct momentFSRM {
      constexpr static unsigned const Shape[1] = {9};
      constexpr static unsigned const Size = 9;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct momentNRF {
      constexpr static unsigned const Shape[1] = {9};
      constexpr static unsigned const Size = 6;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct projectQP {
      constexpr static unsigned const Shape[2] = {56, 343};
      constexpr static unsigned const Size = 19208;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct rDivM {
      constexpr static unsigned const Shape[][2] = {{56, 21}, {56, 21}, {56, 21}, {56, 21}};
      constexpr static unsigned const Size[] = {1176, 1176, 1176, 1176};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct rT {
      constexpr static unsigned const Shape[][2] = {{21, 56}, {21, 56}, {21, 56}, {21, 56}};
      constexpr static unsigned const Size[] = {1344, 1344, 1344, 1344};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct replicateInitialLoading {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct secondInvariant {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct selectBulkAverage {
      constexpr static unsigned const Shape[1] = {6};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct selectBulkNegative {
      constexpr static unsigned const Shape[1] = {6};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct selectVelocity {
      constexpr static unsigned const Shape[2] = {9, 3};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct star {
      constexpr static unsigned const Shape[][2] = {{9, 9}, {9, 9}, {9, 9}};
      constexpr static unsigned const Size[] = {81, 81, 81};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[3];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct subTriangleDofs {
      constexpr static unsigned const Shape[][2] = {{1, 3}, {4, 3}, {16, 3}, {64, 3}};
      constexpr static unsigned const Size[] = {24, 24, 48, 192};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct subTriangleProjection {
      constexpr static unsigned const Shape[][2] = {{1, 56}, {4, 56}, {16, 56}, {64, 56}};
      constexpr static unsigned const Size[] = {448, 448, 896, 3584};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct v {
      constexpr static unsigned const Shape[2] = {56, 56};
      constexpr static unsigned const Size = 3136;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct vInv {
      constexpr static unsigned const Shape[2] = {56, 56};
      constexpr static unsigned const Size = 3136;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct weightSecondInvariant {
      constexpr static unsigned const Shape[1] = {6};
      constexpr static unsigned const Size = 6;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct yieldFactor {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
  }
}
#endif
