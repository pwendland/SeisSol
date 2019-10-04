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
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct AplusT : tensor::AplusT {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct I : tensor::I {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 9}, {0, 0}, {56, 9});
        }
      };
    };
    struct Q : tensor::Q {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 9}, {0, 0}, {56, 9});
        }
      };
    };
    struct QAtPoint : tensor::QAtPoint {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {9};

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {9}, {0}, {9});
        }
      };
    };
    struct QFortran : tensor::QFortran {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 9}, {0, 0}, {56, 9});
        }
      };
    };
    struct QStress : tensor::QStress {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 6};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 6}, {0, 0}, {56, 6});
        }
      };
    };
    struct QStressNodal : tensor::QStressNodal {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 6};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 6}, {0, 0}, {56, 6});
        }
      };
    };
    struct QgodLocal : tensor::QgodLocal {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct QgodNeighbor : tensor::QgodNeighbor {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct T : tensor::T {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct Tinv : tensor::Tinv {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct V3mTo2n : tensor::V3mTo2n {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {52, 56};
      constexpr static unsigned const Start4[] = {0, 0};
      constexpr static unsigned const Stop4[] = {52, 56};
      constexpr static unsigned const Start8[] = {0, 0};
      constexpr static unsigned const Stop8[] = {52, 56};
      constexpr static unsigned const Start12[] = {0, 0};
      constexpr static unsigned const Stop12[] = {52, 56};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {52, 56};
      constexpr static unsigned const Start5[] = {0, 0};
      constexpr static unsigned const Stop5[] = {52, 56};
      constexpr static unsigned const Start9[] = {0, 0};
      constexpr static unsigned const Stop9[] = {52, 56};
      constexpr static unsigned const Start13[] = {0, 0};
      constexpr static unsigned const Stop13[] = {52, 56};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {52, 56};
      constexpr static unsigned const Start6[] = {0, 0};
      constexpr static unsigned const Stop6[] = {52, 56};
      constexpr static unsigned const Start10[] = {0, 0};
      constexpr static unsigned const Stop10[] = {52, 56};
      constexpr static unsigned const Start14[] = {0, 0};
      constexpr static unsigned const Stop14[] = {52, 56};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {52, 56};
      constexpr static unsigned const Start7[] = {0, 0};
      constexpr static unsigned const Stop7[] = {52, 56};
      constexpr static unsigned const Start11[] = {0, 0};
      constexpr static unsigned const Stop11[] = {52, 56};
      constexpr static unsigned const Start15[] = {0, 0};
      constexpr static unsigned const Stop15[] = {52, 56};
      static double const Values0[] __attribute__((aligned(32)));
      static double const Values4[] __attribute__((aligned(32)));
      static double const Values8[] __attribute__((aligned(32)));
      static double const Values12[] __attribute__((aligned(32)));
      static double const Values1[] __attribute__((aligned(32)));
      static double const Values5[] __attribute__((aligned(32)));
      static double const Values9[] __attribute__((aligned(32)));
      static double const Values13[] __attribute__((aligned(32)));
      static double const Values2[] __attribute__((aligned(32)));
      static double const Values6[] __attribute__((aligned(32)));
      static double const Values10[] __attribute__((aligned(32)));
      static double const Values14[] __attribute__((aligned(32)));
      static double const Values3[] __attribute__((aligned(32)));
      static double const Values7[] __attribute__((aligned(32)));
      static double const Values11[] __attribute__((aligned(32)));
      static double const Values15[] __attribute__((aligned(32)));
      static double const* Values[];

      template<unsigned i0, unsigned i1> struct view {};
    };
    template<>
    struct V3mTo2n::view<0,0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<0,1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<0,2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<0,3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<1,0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<1,1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<1,2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<1,3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<2,0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<2,1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<2,2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<2,3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<3,0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<3,1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<3,2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
      }
    };
    template<>
    struct V3mTo2n::view<3,3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 56}, {0, 0}, {52, 56});
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
      static double const Values0[] __attribute__((aligned(32)));
      static double const Values4[] __attribute__((aligned(32)));
      static double const Values8[] __attribute__((aligned(32)));
      static double const Values12[] __attribute__((aligned(32)));
      static double const Values1[] __attribute__((aligned(32)));
      static double const Values5[] __attribute__((aligned(32)));
      static double const Values9[] __attribute__((aligned(32)));
      static double const Values13[] __attribute__((aligned(32)));
      static double const Values2[] __attribute__((aligned(32)));
      static double const Values6[] __attribute__((aligned(32)));
      static double const Values10[] __attribute__((aligned(32)));
      static double const Values14[] __attribute__((aligned(32)));
      static double const Values3[] __attribute__((aligned(32)));
      static double const Values7[] __attribute__((aligned(32)));
      static double const Values11[] __attribute__((aligned(32)));
      static double const Values15[] __attribute__((aligned(32)));
      static double const* Values[];

      template<unsigned i0, unsigned i1> struct view {};
    };
    template<>
    struct V3mTo2nTWDivM::view<0,0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<0,1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<0,2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<0,3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<1,0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<1,1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<1,2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<1,3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<2,0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<2,1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<2,2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<2,3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<3,0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<3,1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<3,2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    template<>
    struct V3mTo2nTWDivM::view<3,3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 49}, {0, 0}, {56, 49});
      }
    };
    struct basisFunctions : tensor::basisFunctions {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
    struct dQ : tensor::dQ {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {56, 9};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {36, 9};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {20, 9};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {12, 9};
      constexpr static unsigned const Start4[] = {0, 0};
      constexpr static unsigned const Stop4[] = {4, 9};
      constexpr static unsigned const Start5[] = {0, 0};
      constexpr static unsigned const Stop5[] = {4, 9};

      template<unsigned i0> struct view {};
    };
    template<>
    struct dQ::view<0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 9}, {0, 0}, {56, 9});
      }
    };
    template<>
    struct dQ::view<1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 9}, {0, 0}, {36, 9});
      }
    };
    template<>
    struct dQ::view<2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 9}, {0, 0}, {20, 9});
      }
    };
    template<>
    struct dQ::view<3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 9}, {0, 0}, {12, 9});
      }
    };
    template<>
    struct dQ::view<4> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 9}, {0, 0}, {4, 9});
      }
    };
    template<>
    struct dQ::view<5> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 9}, {0, 0}, {4, 9});
      }
    };
    struct displacement : tensor::displacement {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 3};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 3}, {0, 0}, {56, 3});
        }
      };
    };
    struct dofsQP : tensor::dofsQP {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {344, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {343, 9}, {0, 0}, {344, 9});
        }
      };
    };
    struct evalAtQP : tensor::evalAtQP {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {344, 56};
      static double const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {343, 56}, {0, 0}, {344, 56});
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
      static double const Values0[] __attribute__((aligned(32)));
      static double const Values1[] __attribute__((aligned(32)));
      static double const Values2[] __attribute__((aligned(32)));
      static double const Values3[] __attribute__((aligned(32)));
      static double const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct fMrT::view<0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct fMrT::view<1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct fMrT::view<2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct fMrT::view<3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    struct fP : tensor::fP {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {24, 21};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {24, 21};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {24, 21};
      static double const Values0[] __attribute__((aligned(32)));
      static double const Values1[] __attribute__((aligned(32)));
      static double const Values2[] __attribute__((aligned(32)));
      static double const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct fP::view<0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 21}, {0, 0}, {24, 21});
      }
    };
    template<>
    struct fP::view<1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 21}, {0, 0}, {24, 21});
      }
    };
    template<>
    struct fP::view<2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 21}, {0, 0}, {24, 21});
      }
    };
    struct fluxSolver : tensor::fluxSolver {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct godunovMatrix : tensor::godunovMatrix {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {9, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
        }
      };
    };
    struct godunovState : tensor::godunovState {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {52, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {49, 9}, {0, 0}, {52, 9});
        }
      };
    };
    struct iniCond : tensor::iniCond {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {344, 9};

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {343, 9}, {0, 0}, {344, 9});
        }
      };
    };
    struct initialLoading : tensor::initialLoading {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {6};

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {6}, {0}, {6});
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
      static double const Values0[] __attribute__((aligned(32)));
      static double const Values1[] __attribute__((aligned(32)));
      static double const Values2[] __attribute__((aligned(32)));
      static double const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct kDivM::view<0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 56}, {0, 0}, {56, 35});
      }
    };
    template<>
    struct kDivM::view<1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 56}, {0, 0}, {56, 35});
      }
    };
    template<>
    struct kDivM::view<2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 56}, {0, 0}, {56, 35});
      }
    };
    struct kDivMT : tensor::kDivMT {
      constexpr static unsigned const Start0[] = {0, 1};
      constexpr static unsigned const Stop0[] = {36, 54};
      constexpr static unsigned const Start1[] = {0, 1};
      constexpr static unsigned const Stop1[] = {36, 55};
      constexpr static unsigned const Start2[] = {0, 1};
      constexpr static unsigned const Stop2[] = {36, 56};
      static double const Values0[] __attribute__((aligned(32)));
      static double const Values1[] __attribute__((aligned(32)));
      static double const Values2[] __attribute__((aligned(32)));
      static double const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct kDivMT::view<0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 56}, {0, 1}, {36, 54});
      }
    };
    template<>
    struct kDivMT::view<1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 56}, {0, 1}, {36, 55});
      }
    };
    template<>
    struct kDivMT::view<2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 56}, {0, 1}, {36, 56});
      }
    };
    struct mInvJInvPhisAtSources : tensor::mInvJInvPhisAtSources {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
    struct meanStress : tensor::meanStress {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
    struct momentFSRM : tensor::momentFSRM {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {9};

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {9}, {0}, {9});
        }
      };
    };
    struct momentNRF : tensor::momentNRF {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {6};

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {9}, {0}, {6});
        }
      };
    };
    struct projectQP : tensor::projectQP {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 343};
      static double const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 343}, {0, 0}, {56, 343});
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
      static double const Values0[] __attribute__((aligned(32)));
      static double const Values1[] __attribute__((aligned(32)));
      static double const Values2[] __attribute__((aligned(32)));
      static double const Values3[] __attribute__((aligned(32)));
      static double const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct rDivM::view<0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 21}, {0, 0}, {56, 21});
      }
    };
    template<>
    struct rDivM::view<1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 21}, {0, 0}, {56, 21});
      }
    };
    template<>
    struct rDivM::view<2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 21}, {0, 0}, {56, 21});
      }
    };
    template<>
    struct rDivM::view<3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 21}, {0, 0}, {56, 21});
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
      static double const Values0[] __attribute__((aligned(32)));
      static double const Values1[] __attribute__((aligned(32)));
      static double const Values2[] __attribute__((aligned(32)));
      static double const Values3[] __attribute__((aligned(32)));
      static double const* Values[];

      template<unsigned i0> struct view {};
    };
    template<>
    struct rT::view<0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct rT::view<1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct rT::view<2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    template<>
    struct rT::view<3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {21, 56}, {0, 0}, {24, 56});
      }
    };
    struct replicateInitialLoading : tensor::replicateInitialLoading {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};
      static double const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
    struct secondInvariant : tensor::secondInvariant {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
    struct selectBulkAverage : tensor::selectBulkAverage {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {3};
      static double const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {6}, {0}, {3});
        }
      };
    };
    struct selectBulkNegative : tensor::selectBulkNegative {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {3};
      static double const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {6}, {0}, {3});
        }
      };
    };
    struct selectVelocity : tensor::selectVelocity {
      constexpr static unsigned const RowInd[] = {6, 7, 8};
      constexpr static unsigned const ColPtr[] = {0, 1, 2, 3};
      static double const Values[];

      struct view {
        typedef ::yateto::CSCMatrixView<double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::CSCMatrixView<double,unsigned>(values, {9, 3}, RowInd, ColPtr);
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
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
      }
    };
    template<>
    struct star::view<1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
      }
    };
    template<>
    struct star::view<2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {9, 9}, {0, 0}, {9, 9});
      }
    };
    struct subTriangleDofs : tensor::subTriangleDofs {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {4, 3};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {4, 3};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {16, 3};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {64, 3};

      template<unsigned i0> struct view {};
    };
    template<>
    struct subTriangleDofs::view<0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {1, 3}, {0, 0}, {4, 3});
      }
    };
    template<>
    struct subTriangleDofs::view<1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {4, 3}, {0, 0}, {4, 3});
      }
    };
    template<>
    struct subTriangleDofs::view<2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {16, 3}, {0, 0}, {16, 3});
      }
    };
    template<>
    struct subTriangleDofs::view<3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {64, 3}, {0, 0}, {64, 3});
      }
    };
    struct subTriangleProjection : tensor::subTriangleProjection {
      constexpr static unsigned const Start0[] = {0, 0};
      constexpr static unsigned const Stop0[] = {4, 56};
      constexpr static unsigned const Start1[] = {0, 0};
      constexpr static unsigned const Stop1[] = {4, 56};
      constexpr static unsigned const Start2[] = {0, 0};
      constexpr static unsigned const Stop2[] = {16, 56};
      constexpr static unsigned const Start3[] = {0, 0};
      constexpr static unsigned const Stop3[] = {64, 56};

      template<unsigned i0> struct view {};
    };
    template<>
    struct subTriangleProjection::view<0> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {1, 56}, {0, 0}, {4, 56});
      }
    };
    template<>
    struct subTriangleProjection::view<1> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {4, 56}, {0, 0}, {4, 56});
      }
    };
    template<>
    struct subTriangleProjection::view<2> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {16, 56}, {0, 0}, {16, 56});
      }
    };
    template<>
    struct subTriangleProjection::view<3> {
      typedef ::yateto::DenseTensorView<2,double,unsigned> type;
      static inline type create(double* values) {
        return ::yateto::DenseTensorView<2,double,unsigned>(values, {64, 56}, {0, 0}, {64, 56});
      }
    };
    struct v : tensor::v {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 56};
      static double const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 56}, {0, 0}, {56, 56});
        }
      };
    };
    struct vInv : tensor::vInv {
      constexpr static unsigned const Start[] = {0, 0};
      constexpr static unsigned const Stop[] = {56, 56};
      static double const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<2,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<2,double,unsigned>(values, {56, 56}, {0, 0}, {56, 56});
        }
      };
    };
    struct weightSecondInvariant : tensor::weightSecondInvariant {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {6};
      static double const Values[] __attribute__((aligned(32)));

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {6}, {0}, {6});
        }
      };
    };
    struct yieldFactor : tensor::yieldFactor {
      constexpr static unsigned const Start[] = {0};
      constexpr static unsigned const Stop[] = {56};

      struct view {
        typedef ::yateto::DenseTensorView<1,double,unsigned> type;
        static inline type create(double* values) {
          return ::yateto::DenseTensorView<1,double,unsigned>(values, {56}, {0}, {56});
        }
      };
    };
  }
}
#endif
