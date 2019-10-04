#ifndef SEISSOL_KERNEL_H_
#define SEISSOL_KERNEL_H_
#include <cmath>
#include <limits>
#include "tensor.h"
namespace seissol {
  namespace kernel {
    struct computeFluxSolverLocal {
      constexpr static unsigned long const NonZeroFlops = 3186;
      constexpr static unsigned long const HardwareFlops = 4374;

      double fluxScale = std::numeric_limits<double>::signaling_NaN();
      double* AplusT{};
      double const* QgodLocal{};
      double const* T{};
      double const* Tinv{};
      tensor::star::Container<double const*> star;

      void execute();
    };
  }
  namespace kernel {
    struct computeFluxSolverNeighbor {
      constexpr static unsigned long const NonZeroFlops = 3186;
      constexpr static unsigned long const HardwareFlops = 4374;

      double fluxScale = std::numeric_limits<double>::signaling_NaN();
      double* AminusT{};
      double const* QgodNeighbor{};
      double const* T{};
      double const* Tinv{};
      tensor::star::Container<double const*> star;

      void execute();
    };
  }
  namespace kernel {
    struct copyQToQFortran {
      constexpr static unsigned long const NonZeroFlops = 0;
      constexpr static unsigned long const HardwareFlops = 0;

      double const* Q{};
      double* QFortran{};

      void execute();
    };
  }
  namespace kernel {
    struct projectIniCond {
      constexpr static unsigned long const NonZeroFlops = 345240;
      constexpr static unsigned long const HardwareFlops = 345744;

      double* Q{};
      double const* iniCond{};
      double const* projectQP{};

      void execute();
    };
  }
  namespace kernel {
    struct evalAtQP {
      constexpr static unsigned long const NonZeroFlops = 342657;
      constexpr static unsigned long const HardwareFlops = 346752;

      double const* Q{};
      double* dofsQP{};
      double const* evalAtQP{};

      void execute();
    };
  }
  namespace kernel {
    struct volume {
      constexpr static unsigned long const NonZeroFlops = 34839;
      constexpr static unsigned long const HardwareFlops = 123336;

      double const* I{};
      double* Q{};
      tensor::kDivM::Container<double const*> kDivM;
      tensor::star::Container<double const*> star;

      void execute();
    };
  }
  namespace kernel {
    struct rotateGodunovStateLocal {
      constexpr static unsigned long const NonZeroFlops = 1377;
      constexpr static unsigned long const HardwareFlops = 1458;

      double const* QgodLocal{};
      double const* Tinv{};
      double* godunovMatrix{};

      void execute();
    };
  }
  namespace kernel {
    struct rotateGodunovStateNeighbor {
      constexpr static unsigned long const NonZeroFlops = 1377;
      constexpr static unsigned long const HardwareFlops = 1458;

      double const* QgodNeighbor{};
      double const* Tinv{};
      double* godunovMatrix{};

      void execute();
    };
  }
  namespace kernel {
    struct rotateFluxMatrix {
      constexpr static unsigned long const NonZeroFlops = 432;
      constexpr static unsigned long const HardwareFlops = 1458;

      double fluxScale = std::numeric_limits<double>::signaling_NaN();
      double const* T{};
      double* fluxSolver{};
      tensor::star::Container<double const*> star;

      void execute();
    };
  }
  namespace kernel {
    struct plConvertToNodal {
      constexpr static unsigned long const NonZeroFlops = 33720;
      constexpr static unsigned long const HardwareFlops = 38304;

      double const* QStress{};
      double* QStressNodal{};
      double const* initialLoading{};
      double const* replicateInitialLoading{};
      double const* v{};

      void execute();
    };
  }
  namespace kernel {
    struct plComputeMean {
      constexpr static unsigned long const NonZeroFlops = 280;
      constexpr static unsigned long const HardwareFlops = 336;

      double const* QStressNodal{};
      double* meanStress{};
      double const* selectBulkAverage{};

      void execute();
    };
  }
  namespace kernel {
    struct plSubtractMean {
      constexpr static unsigned long const NonZeroFlops = 336;
      constexpr static unsigned long const HardwareFlops = 336;

      double* QStressNodal{};
      double const* meanStress{};
      double const* selectBulkNegative{};

      void execute();
    };
  }
  namespace kernel {
    struct plComputeSecondInvariant {
      constexpr static unsigned long const NonZeroFlops = 952;
      constexpr static unsigned long const HardwareFlops = 1008;

      double const* QStressNodal{};
      double* secondInvariant{};
      double const* weightSecondInvariant{};

      void execute();
    };
  }
  namespace kernel {
    struct plAdjustStresses {
      constexpr static unsigned long const NonZeroFlops = 31812;
      constexpr static unsigned long const HardwareFlops = 37968;

      double* QStress{};
      double const* QStressNodal{};
      double const* vInv{};
      double const* yieldFactor{};

      void execute();
    };
  }
  namespace kernel {
    struct addVelocity {
      constexpr static unsigned long const NonZeroFlops = 336;
      constexpr static unsigned long const HardwareFlops = 336;

      double const* I{};
      double* displacement{};
      double const* selectVelocity{};

      void execute();
    };
  }
  namespace kernel {
    struct sourceNRF {
      constexpr static unsigned long const NonZeroFlops = 672;
      constexpr static unsigned long const HardwareFlops = 672;

      double* Q{};
      double const* mInvJInvPhisAtSources{};
      double const* momentNRF{};

      void execute();
    };
  }
  namespace kernel {
    struct sourceFSRM {
      constexpr static unsigned long const NonZeroFlops = 1512;
      constexpr static unsigned long const HardwareFlops = 1512;

      double stfIntegral = std::numeric_limits<double>::signaling_NaN();
      double* Q{};
      double const* mInvJInvPhisAtSources{};
      double const* momentFSRM{};

      void execute();
    };
  }
  namespace kernel {
    struct evaluateDOFSAtPoint {
      constexpr static unsigned long const NonZeroFlops = 999;
      constexpr static unsigned long const HardwareFlops = 1008;

      double const* Q{};
      double* QAtPoint{};
      double const* basisFunctions{};

      void execute();
    };
  }
  namespace kernel {
    struct localFlux {
      constexpr static unsigned long const NonZeroFlops[] = {9936, 10080, 31968, 27216};
      constexpr static unsigned long const HardwareFlops[] = {49248, 49248, 49248, 49248};

      double const* AplusT{};
      double const* I{};
      double* Q{};
      tensor::fMrT::Container<double const*> fMrT;
      tensor::rDivM::Container<double const*> rDivM;

      struct Prefetch {
        double const* I{};
        double const* Q{};
      };
      Prefetch _prefetch;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      typedef void (localFlux::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {&localFlux::execute0, &localFlux::execute1, &localFlux::execute2, &localFlux::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
    };
  }
  namespace kernel {
    struct neighboringFlux {
      constexpr static unsigned long const NonZeroFlops[] = {11349, 10125, 11349, 11421, 10197, 11421, 22365, 21141, 22365, 19989, 18765, 19989, 11421, 10197, 11421, 11493, 10269, 11493, 22437, 21213, 22437, 20061, 18837, 20061, 22365, 21141, 22365, 22437, 21213, 22437, 33381, 32157, 33381, 31005, 29781, 31005, 19989, 18765, 19989, 20061, 18837, 20061, 31005, 29781, 31005, 28629, 27405, 28629};
      constexpr static unsigned long const HardwareFlops[] = {58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320, 58320};

      double const* AminusT{};
      double const* I{};
      double* Q{};
      tensor::fP::Container<double const*> fP;
      tensor::rDivM::Container<double const*> rDivM;
      tensor::rT::Container<double const*> rT;

      struct Prefetch {
        double const* I{};
      };
      Prefetch _prefetch;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      void execute4();
      void execute5();
      void execute6();
      void execute7();
      void execute8();
      void execute9();
      void execute10();
      void execute11();
      void execute12();
      void execute13();
      void execute14();
      void execute15();
      void execute16();
      void execute17();
      void execute18();
      void execute19();
      void execute20();
      void execute21();
      void execute22();
      void execute23();
      void execute24();
      void execute25();
      void execute26();
      void execute27();
      void execute28();
      void execute29();
      void execute30();
      void execute31();
      void execute32();
      void execute33();
      void execute34();
      void execute35();
      void execute36();
      void execute37();
      void execute38();
      void execute39();
      void execute40();
      void execute41();
      void execute42();
      void execute43();
      void execute44();
      void execute45();
      void execute46();
      void execute47();
      typedef void (neighboringFlux::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {&neighboringFlux::execute0, &neighboringFlux::execute1, &neighboringFlux::execute2, &neighboringFlux::execute3, &neighboringFlux::execute4, &neighboringFlux::execute5, &neighboringFlux::execute6, &neighboringFlux::execute7, &neighboringFlux::execute8, &neighboringFlux::execute9, &neighboringFlux::execute10, &neighboringFlux::execute11, &neighboringFlux::execute12, &neighboringFlux::execute13, &neighboringFlux::execute14, &neighboringFlux::execute15, &neighboringFlux::execute16, &neighboringFlux::execute17, &neighboringFlux::execute18, &neighboringFlux::execute19, &neighboringFlux::execute20, &neighboringFlux::execute21, &neighboringFlux::execute22, &neighboringFlux::execute23, &neighboringFlux::execute24, &neighboringFlux::execute25, &neighboringFlux::execute26, &neighboringFlux::execute27, &neighboringFlux::execute28, &neighboringFlux::execute29, &neighboringFlux::execute30, &neighboringFlux::execute31, &neighboringFlux::execute32, &neighboringFlux::execute33, &neighboringFlux::execute34, &neighboringFlux::execute35, &neighboringFlux::execute36, &neighboringFlux::execute37, &neighboringFlux::execute38, &neighboringFlux::execute39, &neighboringFlux::execute40, &neighboringFlux::execute41, &neighboringFlux::execute42, &neighboringFlux::execute43, &neighboringFlux::execute44, &neighboringFlux::execute45, &neighboringFlux::execute46, &neighboringFlux::execute47};
      constexpr static member_function_ptr findExecute(unsigned i0, unsigned i1, unsigned i2) {
        return ExecutePtrs[1*i0 + 3*i1 + 12*i2];
      }
      inline void execute(unsigned i0, unsigned i1, unsigned i2) {
        (this->*findExecute(i0, i1, i2))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0, unsigned i1, unsigned i2) {
        return NonZeroFlops[1*i0 + 3*i1 + 12*i2];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0, unsigned i1, unsigned i2) {
        return HardwareFlops[1*i0 + 3*i1 + 12*i2];
      }
    };
  }
  namespace kernel {
    struct derivativeTaylorExpansion {
      constexpr static unsigned long const NonZeroFlops[] = {504, 630, 360, 180, 72, 18};
      constexpr static unsigned long const HardwareFlops[] = {504, 630, 360, 180, 72, 18};

      double power = std::numeric_limits<double>::signaling_NaN();
      double* I{};
      tensor::dQ::Container<double const*> dQ;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      void execute4();
      void execute5();
      typedef void (derivativeTaylorExpansion::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {&derivativeTaylorExpansion::execute0, &derivativeTaylorExpansion::execute1, &derivativeTaylorExpansion::execute2, &derivativeTaylorExpansion::execute3, &derivativeTaylorExpansion::execute4, &derivativeTaylorExpansion::execute5};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
    };
  }
  namespace kernel {
    struct derivative {
      constexpr static unsigned long const NonZeroFlops[] = {0, 34524, 13806, 4716, 1260, 216};
      constexpr static unsigned long const HardwareFlops[] = {0, 122472, 45360, 17496, 3672, 2376};

      tensor::dQ::Container<double*> dQ;
      tensor::kDivMT::Container<double const*> kDivMT;
      tensor::star::Container<double const*> star;

      void execute1();
      void execute2();
      void execute3();
      void execute4();
      void execute5();
      typedef void (derivative::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {nullptr, &derivative::execute1, &derivative::execute2, &derivative::execute3, &derivative::execute4, &derivative::execute5};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
    };
  }
  namespace kernel {
    struct godunovState {
      constexpr static unsigned long const NonZeroFlops[] = {53676, 56448, 56448, 56448, 56889, 54117, 56889, 56889, 56889, 56889, 56889, 56889, 54117, 56889, 56889, 56889};
      constexpr static unsigned long const HardwareFlops[] = {60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840};

      double const* Q{};
      tensor::V3mTo2n::Container<double const*> V3mTo2n;
      double const* godunovMatrix{};
      double* godunovState{};

      struct Prefetch {
        double const* godunovState{};
      };
      Prefetch _prefetch;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      void execute4();
      void execute5();
      void execute6();
      void execute7();
      void execute8();
      void execute9();
      void execute10();
      void execute11();
      void execute12();
      void execute13();
      void execute14();
      void execute15();
      typedef void (godunovState::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {&godunovState::execute0, &godunovState::execute1, &godunovState::execute2, &godunovState::execute3, &godunovState::execute4, &godunovState::execute5, &godunovState::execute6, &godunovState::execute7, &godunovState::execute8, &godunovState::execute9, &godunovState::execute10, &godunovState::execute11, &godunovState::execute12, &godunovState::execute13, &godunovState::execute14, &godunovState::execute15};
      constexpr static member_function_ptr findExecute(unsigned i0, unsigned i1) {
        return ExecutePtrs[1*i0 + 4*i1];
      }
      inline void execute(unsigned i0, unsigned i1) {
        (this->*findExecute(i0, i1))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0, unsigned i1) {
        return NonZeroFlops[1*i0 + 4*i1];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0, unsigned i1) {
        return HardwareFlops[1*i0 + 4*i1];
      }
    };
  }
  namespace kernel {
    struct nodalFlux {
      constexpr static unsigned long const NonZeroFlops[] = {54117, 56889, 56889, 56889, 56889, 54117, 56889, 56889, 56889, 56889, 56889, 56889, 54117, 56889, 56889, 56889};
      constexpr static unsigned long const HardwareFlops[] = {57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816};

      double* Q{};
      tensor::V3mTo2nTWDivM::Container<double const*> V3mTo2nTWDivM;
      double const* fluxSolver{};
      double const* godunovState{};

      struct Prefetch {
        double const* I{};
      };
      Prefetch _prefetch;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      void execute4();
      void execute5();
      void execute6();
      void execute7();
      void execute8();
      void execute9();
      void execute10();
      void execute11();
      void execute12();
      void execute13();
      void execute14();
      void execute15();
      typedef void (nodalFlux::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {&nodalFlux::execute0, &nodalFlux::execute1, &nodalFlux::execute2, &nodalFlux::execute3, &nodalFlux::execute4, &nodalFlux::execute5, &nodalFlux::execute6, &nodalFlux::execute7, &nodalFlux::execute8, &nodalFlux::execute9, &nodalFlux::execute10, &nodalFlux::execute11, &nodalFlux::execute12, &nodalFlux::execute13, &nodalFlux::execute14, &nodalFlux::execute15};
      constexpr static member_function_ptr findExecute(unsigned i0, unsigned i1) {
        return ExecutePtrs[1*i0 + 4*i1];
      }
      inline void execute(unsigned i0, unsigned i1) {
        (this->*findExecute(i0, i1))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0, unsigned i1) {
        return NonZeroFlops[1*i0 + 4*i1];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0, unsigned i1) {
        return HardwareFlops[1*i0 + 4*i1];
      }
    };
  }
  namespace kernel {
    struct subTriangleDisplacement {
      constexpr static unsigned long const NonZeroFlops[] = {333, 1332, 5328, 21312};
      constexpr static unsigned long const HardwareFlops[] = {1344, 1344, 5376, 21504};

      double const* displacement{};
      tensor::subTriangleDofs::Container<double*> subTriangleDofs;
      tensor::subTriangleProjection::Container<double const*> subTriangleProjection;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      typedef void (subTriangleDisplacement::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {&subTriangleDisplacement::execute0, &subTriangleDisplacement::execute1, &subTriangleDisplacement::execute2, &subTriangleDisplacement::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
    };
  }
  namespace kernel {
    struct subTriangleVelocity {
      constexpr static unsigned long const NonZeroFlops[] = {336, 1344, 5376, 21480};
      constexpr static unsigned long const HardwareFlops[] = {1368, 1368, 5472, 21840};

      double const* Q{};
      double const* selectVelocity{};
      tensor::subTriangleDofs::Container<double*> subTriangleDofs;
      tensor::subTriangleProjection::Container<double const*> subTriangleProjection;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      typedef void (subTriangleVelocity::* const member_function_ptr)(void);
      constexpr static member_function_ptr ExecutePtrs[] = {&subTriangleVelocity::execute0, &subTriangleVelocity::execute1, &subTriangleVelocity::execute2, &subTriangleVelocity::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
    };
  }
}
#endif
