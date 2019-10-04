#ifndef SEISSOL_INIT_H_
#define SEISSOL_INIT_H_
#include "tensor.h"
#include "yateto.h"
namespace seissol {
  namespace init {
    struct AminusT : tensor::AminusT {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct AplusT : tensor::AplusT {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct I : tensor::I {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 9}, {0, 0}, {56, 9});
        }
      };
    };
    struct Q : tensor::Q {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 9}, {0, 0}, {56, 9});
        }
      };
    };
    struct QAtPoint : tensor::QAtPoint {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {9};

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {9}, {0}, {9});
        }
      };
    };
    struct QFortran : tensor::QFortran {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 9}, {0, 0}, {56, 9});
        }
      };
    };
    struct QStress : tensor::QStress {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 6};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 6}, {0, 0}, {56, 6});
        }
      };
    };
    struct QStressNodal : tensor::QStressNodal {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 6};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 6}, {0, 0}, {56, 6});
        }
      };
    };
    struct QgodLocal : tensor::QgodLocal {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct QgodNeighbor : tensor::QgodNeighbor {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct T : tensor::T {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct Tinv : tensor::Tinv {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct V3mTo2n : tensor::V3mTo2n {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {56, 56};
      constexpr static unsigned const Start4[] = {0, 0};
      constexpr static unsigned const Stop4[] = {56, 56};
      constexpr static unsigned const Start8[] = {0, 0};
      constexpr static unsigned const Stop8[] = {56, 56};
      constexpr static unsigned const Start12[] = {0, 0};
      constexpr static unsigned const Stop12[] = {56, 56};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {56, 56};
      constexpr static unsigned const Start5[] = {0, 0};
      constexpr static unsigned const Stop5[] = {56, 56};
      constexpr static unsigned const Start9[] = {0, 0};
      constexpr static unsigned const Stop9[] = {56, 56};
      constexpr static unsigned const Start13[] = {0, 0};
      constexpr static unsigned const Stop13[] = {56, 56};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {56, 56};
      constexpr static unsigned const Start6[] = {0, 0};
      constexpr static unsigned const Stop6[] = {56, 56};
      constexpr static unsigned const Start10[] = {0, 0};
      constexpr static unsigned const Stop10[] = {56, 56};
      constexpr static unsigned const Start14[] = {0, 0};
      constexpr static unsigned const Stop14[] = {56, 56};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {56, 56};
      constexpr static unsigned const Start7[] = {0, 0};
      constexpr static unsigned const Stop7[] = {56, 56};
      constexpr static unsigned const Start11[] = {0, 0};
      constexpr static unsigned const Stop11[] = {56, 56};
      constexpr static unsigned const Start15[] = {0, 0};
      constexpr static unsigned const Stop15[] = {56, 56};
      static float const Values0[] __attribute__((aligned(32)));
      static float const Values4[] __attribute__((aligned(32)));
      static float const Values8[] __attribute__((aligned(32)));
      static float const Values12[] __attribute__((aligned(32)));
      static float const Values1[] __attribute__((aligned(32)));
      static float const Values5[] __attribute__((aligned(32)));
      static float const Values9[] __attribute__((aligned(32)));
      static float const Values13[] __attribute__((aligned(32)));
      static float const Values2[] __attribute__((aligned(32)));
      static float const Values6[] __attribute__((aligned(32)));
      static float const Values10[] __attribute__((aligned(32)));
      static float const Values14[] __attribute__((aligned(32)));
      static float const Values3[] __attribute__((aligned(32)));
      static float const Values7[] __attribute__((aligned(32)));
      static float const Values11[] __attribute__((aligned(32)));
      static float const Values15[] __attribute__((aligned(32)));
      static float const* Values[];

      template<unsigned i0, unsigned i1> struct view {};
    };
    template<>
    struct V3mTo2n::view<0,0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<0,1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<0,2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<0,3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<1,0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<1,1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<1,2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<1,3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<2,0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<2,1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<2,2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<2,3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<3,0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<3,1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<3,2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    template<>
    struct V3mTo2n::view<3,3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 56}, {0, 0}, {56, 56});
      }
    };
    struct V3mTo2nTWDivM : tensor::V3mTo2nTWDivM {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {56, 49};
      constexpr static unsigned const Start4[] = {0, 0};
      constexpr static unsigned const Stop4[] = {56, 49};
      constexpr static unsigned const Start8[] = {0, 0};
      constexpr static unsigned const Stop8[] = {56, 49};
      constexpr static unsigned const Start12[] = {0, 0};
      constexpr static unsigned const Stop12[] = {56, 49};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {56, 49};
      constexpr static unsigned const Start5[] = {0, 0};
      constexpr static unsigned const Stop5[] = {56, 49};
      constexpr static unsigned const Start9[] = {0, 0};
      constexpr static unsigned const Stop9[] = {56, 49};
      constexpr static unsigned const Start13[] = {0, 0};
      constexpr static unsigned const Stop13[] = {56, 49};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {56, 49};
      constexpr static unsigned const Start6[] = {0, 0};
      constexpr static unsigned const Stop6[] = {56, 49};
      constexpr static unsigned const Start10[] = {0, 0};
      constexpr static unsigned const Stop10[] = {56, 49};
      constexpr static unsigned const Start14[] = {0, 0};
      constexpr static unsigned const Stop14[] = {56, 49};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {56, 49};
      constexpr static unsigned const Start7[] = {0, 0};
      constexpr static unsigned const Stop7[] = {56, 49};
      constexpr static unsigned const Start11[] = {0, 0};
      constexpr static unsigned const Stop11[] = {56, 49};
      constexpr static unsigned const Start15[] = {0, 0};
      constexpr static unsigned const Stop15[] = {56, 49};
      static float const Values0[] __attribute__((aligned(32)));
      static float const Values4[] __attribute__((aligned(32)));
      static float const Values8[] __attribute__((aligned(32)));
      static float const Values12[] __attribute__((aligned(32)));
      static float const Values1[] __attribute__((aligned(32)));
      static float const Values5[] __attribute__((aligned(32)));
      static float const Values9[] __attribute__((aligned(32)));
      static float const Values13[] __attribute__((aligned(32)));
      static float const Values2[] __attribute__((aligned(32)));
      static float const Values6[] __attribute__((aligned(32)));
      static float const Values10[] __attribute__((aligned(32)));
      static float const Values14[] __attribute__((aligned(32)));
      static float const Values3[] __attribute__((aligned(32)));
      static float const Values7[] __attribute__((aligned(32)));
      static float const Values11[] __attribute__((aligned(32)));
      static float const Values15[] __attribute__((aligned(32)));
      static float const* Values[];

      template<unsigned i0, unsigned i1> struct view {};
    };
    template<>
    struct V3mTo2nTWDivM::view<0,0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<0,1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<0,2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<0,3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<1,0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<1,1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<1,2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<1,3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<2,0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<2,1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<2,2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<2,3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<3,0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<3,1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<3,2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<3,3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    struct basisFunctions : tensor::basisFunctions {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
    struct dQ : tensor::dQ {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {56, 9};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {40, 9};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {24, 9};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {16, 9};
      constexpr static unsigned const Start4[] = {0, 0};
      constexpr static unsigned const Stop4[] = {8, 9};
      constexpr static unsigned const Start5[] = {0, 0};
      constexpr static unsigned const Stop5[] = {8, 9};

      template<unsigned i0> struct view {};
    };
    template<>
    struct dQ::view<0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 9}, {0, 0}, {56, 9});
      }
    };
    template<>
    struct dQ::view<1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 9}, {0, 0}, {40, 9});
      }
    };
    template<>
    struct dQ::view<2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 9}, {0, 0}, {24, 9});
      }
    };
    template<>
    struct dQ::view<3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 9}, {0, 0}, {16, 9});
      }
    };
    template<>
    struct dQ::view<4> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 9}, {0, 0}, {8, 9});
      }
    };
    template<>
    struct dQ::view<5> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 9}, {0, 0}, {8, 9});
      }
    };
    struct displacement : tensor::displacement {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 3};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 3}, {0, 0}, {56, 3});
        }
      };
    };
    struct dofsQP : tensor::dofsQP {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {344, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {343, 9}, {0, 0}, {344, 9});
        }
      };
    };
    struct evalAtQP : tensor::evalAtQP {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {344, 56};
      static float const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {343, 56}, {0, 0}, {344, 56});
        }
      };
    };
    struct fMrT : tensor::fMrT {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {24, 56};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {24, 56};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {24, 56};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {24, 56};
      static float const Values0[] __attribute__((aligned(32)));
      static float const Values1[] __attribute__((aligned(32)));
      static float const Values2[] __attribute__((aligned(32)));
      static float const Values3[] __attribute__((aligned(32)));
      static float const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct fMrT::view<0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct fMrT::view<1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct fMrT::view<2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct fMrT::view<3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    struct fP : tensor::fP {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {24, 21};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {24, 21};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {24, 21};
      static float const Values0[] __attribute__((aligned(32)));
      static float const Values1[] __attribute__((aligned(32)));
      static float const Values2[] __attribute__((aligned(32)));
      static float const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct fP::view<0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 21}, {0, 0}, {24, 21});
      }
    };
    template<>
    struct fP::view<1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 21}, {0, 0}, {24, 21});
      }
    };
    template<>
    struct fP::view<2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 21}, {0, 0}, {24, 21});
      }
    };
    struct fluxSolver : tensor::fluxSolver {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct godunovMatrix : tensor::godunovMatrix {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct godunovState : tensor::godunovState {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {49, 9}, {0, 0}, {56, 9});
        }
      };
    };
    struct iniCond : tensor::iniCond {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {344, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {343, 9}, {0, 0}, {344, 9});
        }
      };
    };
    struct initialLoading : tensor::initialLoading {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {6};

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {6}, {0}, {6});
        }
      };
    };
    struct kDivM : tensor::kDivM {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {56, 35};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {56, 35};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {56, 35};
      static float const Values0[] __attribute__((aligned(32)));
      static float const Values1[] __attribute__((aligned(32)));
      static float const Values2[] __attribute__((aligned(32)));
      static float const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct kDivM::view<0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 56}, {0, 0}, {56, 35});
      }
    };
    template<>
    struct kDivM::view<1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 56}, {0, 0}, {56, 35});
      }
    };
    template<>
    struct kDivM::view<2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 56}, {0, 0}, {56, 35});
      }
    };
    struct kDivMT : tensor::kDivMT {
      constexpr static unsigned const Start0[] = {0, 1};
      constexpr static unsigned const Stop0[] = {40, 54};
      constexpr static unsigned const Start1[] = {0, 1};
      constexpr static unsigned const Stop1[] = {40, 55};
      constexpr static unsigned const Start2[] = {0, 1};
      constexpr static unsigned const Stop2[] = {40, 56};
      static float const Values0[] __attribute__((aligned(32)));
      static float const Values1[] __attribute__((aligned(32)));
      static float const Values2[] __attribute__((aligned(32)));
      static float const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct kDivMT::view<0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 56}, {0, 1}, {40, 54});
      }
    };
    template<>
    struct kDivMT::view<1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 56}, {0, 1}, {40, 55});
      }
    };
    template<>
    struct kDivMT::view<2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 56}, {0, 1}, {40, 56});
      }
    };
    struct mInvJInvPhisAtSources : tensor::mInvJInvPhisAtSources {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
    struct meanStress : tensor::meanStress {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
    struct momentFSRM : tensor::momentFSRM {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {9};

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {9}, {0}, {9});
        }
      };
    };
    struct momentNRF : tensor::momentNRF {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {6};

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {9}, {0}, {6});
        }
      };
    };
    struct projectQP : tensor::projectQP {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 343};
      static float const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 343}, {0, 0}, {56, 343});
        }
      };
    };
    struct rDivM : tensor::rDivM {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {56, 21};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {56, 21};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {56, 21};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {56, 21};
      static float const Values0[] __attribute__((aligned(32)));
      static float const Values1[] __attribute__((aligned(32)));
      static float const Values2[] __attribute__((aligned(32)));
      static float const Values3[] __attribute__((aligned(32)));
      static float const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct rDivM::view<0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 21}, {0, 0}, {56, 21});
      }
    };
    template<>
    struct rDivM::view<1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 21}, {0, 0}, {56, 21});
      }
    };
    template<>
    struct rDivM::view<2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 21}, {0, 0}, {56, 21});
      }
    };
    template<>
    struct rDivM::view<3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 21}, {0, 0}, {56, 21});
      }
    };
    struct rT : tensor::rT {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {24, 56};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {24, 56};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {24, 56};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {24, 56};
      static float const Values0[] __attribute__((aligned(32)));
      static float const Values1[] __attribute__((aligned(32)));
      static float const Values2[] __attribute__((aligned(32)));
      static float const Values3[] __attribute__((aligned(32)));
      static float const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct rT::view<0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct rT::view<1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct rT::view<2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct rT::view<3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    struct replicateInitialLoading : tensor::replicateInitialLoading {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};
      static float const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
    struct secondInvariant : tensor::secondInvariant {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
    struct selectBulkAverage : tensor::selectBulkAverage {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {3};
      static float const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {6}, {0}, {3});
        }
      };
    };
    struct selectBulkNegative : tensor::selectBulkNegative {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {3};
      static float const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {6}, {0}, {3});
        }
      };
    };
    struct selectVelocity : tensor::selectVelocity {
      constexpr static unsigned const RowInd[] = {6, 7, 8};
      constexpr static unsigned const ColPtr[] = {0, 1, 2, 3};
      static float const Values[];

      struct view {
        typedef ::yateto::CSCMatrixView<float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::CSCMatrixView<float,unsigned>(values, {9, 3}, RowInd, ColPtr);
        }
      };
    };
    struct star : tensor::star {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {9, 9};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {9, 9};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {9, 9};

      template<unsigned i0> struct view {};
    };
    template<>
    struct star::view<0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
      }
    };
    template<>
    struct star::view<1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
      }
    };
    template<>
    struct star::view<2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
      }
    };
    struct subTriangleDofs : tensor::subTriangleDofs {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {8, 3};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {8, 3};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {16, 3};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {64, 3};

      template<unsigned i0> struct view {};
    };
    template<>
    struct subTriangleDofs::view<0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {1, 3}, {0, 0}, {8, 3});
      }
    };
    template<>
    struct subTriangleDofs::view<1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {4, 3}, {0, 0}, {8, 3});
      }
    };
    template<>
    struct subTriangleDofs::view<2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {16, 3}, {0, 0}, {16, 3});
      }
    };
    template<>
    struct subTriangleDofs::view<3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {64, 3}, {0, 0}, {64, 3});
      }
    };
    struct subTriangleProjection : tensor::subTriangleProjection {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {8, 56};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {8, 56};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {16, 56};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {64, 56};

      template<unsigned i0> struct view {};
    };
    template<>
    struct subTriangleProjection::view<0> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {1, 56}, {0, 0}, {8, 56});
      }
    };
    template<>
    struct subTriangleProjection::view<1> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {4, 56}, {0, 0}, {8, 56});
      }
    };
    template<>
    struct subTriangleProjection::view<2> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {16, 56}, {0, 0}, {16, 56});
      }
    };
    template<>
    struct subTriangleProjection::view<3> {
      typedef ::yateto::DenseTensorView<2,float,unsigned> type;
      static inline type create(float* values) {
        return ::yateto::DenseTensorView<2,float,unsigned>(values, {64, 56}, {0, 0}, {64, 56});
      }
    };
    struct v : tensor::v {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 56};
      static float const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 56}, {0, 0}, {56, 56});
        }
      };
    };
    struct vInv : tensor::vInv {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 56};
      static float const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<2,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<2,float,unsigned>(values, {56, 56}, {0, 0}, {56, 56});
        }
      };
    };
    struct weightSecondInvariant : tensor::weightSecondInvariant {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {6};
      static float const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {6}, {0}, {6});
        }
      };
    };
    struct yieldFactor : tensor::yieldFactor {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};

      struct view {
        typedef ::yateto::DenseTensorView<1,float,unsigned> type;
        static inline type create(float* values) {
          return ::yateto::DenseTensorView<1,float,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
  }
}
#endif
