#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "kernel.h"
#include "tensor.h"
#include "Stopwatch.h"
#include "Util.h"
using namespace yateto;
void trashTheCache(double* trash, int size);
int main(int argc, char** argv) {
  int _fixedReps = (argc >= 2) ? atoi(argv[1]) : -1;
  int _reps, _error;
  Stopwatch _sw;
  double _time, _nzflops, _flops;
  printf("kernel,repetitions,time,numnzflop,numflop,nzgflops,gflops\n");
  {
    real* AplusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AplusT__), ALIGNMENT, tensor::AplusT::size()*sizeof(real));
    real* Tinv__;
    _error = posix_memalign(reinterpret_cast<void**>(&Tinv__), ALIGNMENT, tensor::Tinv::size()*sizeof(real));
    real* QgodLocal__;
    _error = posix_memalign(reinterpret_cast<void**>(&QgodLocal__), ALIGNMENT, tensor::QgodLocal::size()*sizeof(real));
    real* star__0;
    _error = posix_memalign(reinterpret_cast<void**>(&star__0), ALIGNMENT, tensor::star::size(0)*sizeof(real));
    real* T__;
    _error = posix_memalign(reinterpret_cast<void**>(&T__), ALIGNMENT, tensor::T::size()*sizeof(real));
    fillWithStuff(AplusT__, tensor::AplusT::size());
    fillWithStuff(Tinv__, tensor::Tinv::size());
    fillWithStuff(QgodLocal__, tensor::QgodLocal::size());
    fillWithStuff(star__0, tensor::star::size(0));
    fillWithStuff(T__, tensor::T::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::computeFluxSolverLocal::HardwareFlops);
    }
    kernel::computeFluxSolverLocal _kernel_computeFluxSolverLocal;
    _kernel_computeFluxSolverLocal.AplusT = AplusT__;
    _kernel_computeFluxSolverLocal.Tinv = Tinv__;
    _kernel_computeFluxSolverLocal.QgodLocal = QgodLocal__;
    _kernel_computeFluxSolverLocal.star(0) = star__0;
    _kernel_computeFluxSolverLocal.T = T__;
    _kernel_computeFluxSolverLocal.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel_computeFluxSolverLocal.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::computeFluxSolverLocal::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::computeFluxSolverLocal::HardwareFlops) * _reps / _time / 1.0e9;
    printf("computeFluxSolverLocal,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::computeFluxSolverLocal::NonZeroFlops, kernel::computeFluxSolverLocal::HardwareFlops, _nzflops, _flops);
    free(AplusT__);
    free(Tinv__);
    free(QgodLocal__);
    free(star__0);
    free(T__);
  }
  {
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    real* Tinv__;
    _error = posix_memalign(reinterpret_cast<void**>(&Tinv__), ALIGNMENT, tensor::Tinv::size()*sizeof(real));
    real* QgodNeighbor__;
    _error = posix_memalign(reinterpret_cast<void**>(&QgodNeighbor__), ALIGNMENT, tensor::QgodNeighbor::size()*sizeof(real));
    real* star__0;
    _error = posix_memalign(reinterpret_cast<void**>(&star__0), ALIGNMENT, tensor::star::size(0)*sizeof(real));
    real* T__;
    _error = posix_memalign(reinterpret_cast<void**>(&T__), ALIGNMENT, tensor::T::size()*sizeof(real));
    fillWithStuff(AminusT__, tensor::AminusT::size());
    fillWithStuff(Tinv__, tensor::Tinv::size());
    fillWithStuff(QgodNeighbor__, tensor::QgodNeighbor::size());
    fillWithStuff(star__0, tensor::star::size(0));
    fillWithStuff(T__, tensor::T::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::computeFluxSolverNeighbor::HardwareFlops);
    }
    kernel::computeFluxSolverNeighbor _kernel_computeFluxSolverNeighbor;
    _kernel_computeFluxSolverNeighbor.AminusT = AminusT__;
    _kernel_computeFluxSolverNeighbor.Tinv = Tinv__;
    _kernel_computeFluxSolverNeighbor.QgodNeighbor = QgodNeighbor__;
    _kernel_computeFluxSolverNeighbor.star(0) = star__0;
    _kernel_computeFluxSolverNeighbor.T = T__;
    _kernel_computeFluxSolverNeighbor.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel_computeFluxSolverNeighbor.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::computeFluxSolverNeighbor::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::computeFluxSolverNeighbor::HardwareFlops) * _reps / _time / 1.0e9;
    printf("computeFluxSolverNeighbor,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::computeFluxSolverNeighbor::NonZeroFlops, kernel::computeFluxSolverNeighbor::HardwareFlops, _nzflops, _flops);
    free(AminusT__);
    free(Tinv__);
    free(QgodNeighbor__);
    free(star__0);
    free(T__);
  }
  {
    real* QFortran__;
    _error = posix_memalign(reinterpret_cast<void**>(&QFortran__), ALIGNMENT, tensor::QFortran::size()*sizeof(real));
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    fillWithStuff(QFortran__, tensor::QFortran::size());
    fillWithStuff(Q__, tensor::Q::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::copyQToQFortran::HardwareFlops);
    }
    kernel::copyQToQFortran _kernel_copyQToQFortran;
    _kernel_copyQToQFortran.QFortran = QFortran__;
    _kernel_copyQToQFortran.Q = Q__;
    _kernel_copyQToQFortran.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel_copyQToQFortran.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::copyQToQFortran::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::copyQToQFortran::HardwareFlops) * _reps / _time / 1.0e9;
    printf("copyQToQFortran,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::copyQToQFortran::NonZeroFlops, kernel::copyQToQFortran::HardwareFlops, _nzflops, _flops);
    free(QFortran__);
    free(Q__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* projectQP__;
    _error = posix_memalign(reinterpret_cast<void**>(&projectQP__), ALIGNMENT, tensor::projectQP::size()*sizeof(real));
    real* iniCond__;
    _error = posix_memalign(reinterpret_cast<void**>(&iniCond__), ALIGNMENT, tensor::iniCond::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(projectQP__, tensor::projectQP::size());
    fillWithStuff(iniCond__, tensor::iniCond::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::projectIniCond::HardwareFlops);
    }
    kernel::projectIniCond _kernel_projectIniCond;
    _kernel_projectIniCond.Q = Q__;
    _kernel_projectIniCond.projectQP = projectQP__;
    _kernel_projectIniCond.iniCond = iniCond__;
    _kernel_projectIniCond.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel_projectIniCond.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::projectIniCond::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::projectIniCond::HardwareFlops) * _reps / _time / 1.0e9;
    printf("projectIniCond,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::projectIniCond::NonZeroFlops, kernel::projectIniCond::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(projectQP__);
    free(iniCond__);
  }
  {
    real* dofsQP__;
    _error = posix_memalign(reinterpret_cast<void**>(&dofsQP__), ALIGNMENT, tensor::dofsQP::size()*sizeof(real));
    real* evalAtQP__;
    _error = posix_memalign(reinterpret_cast<void**>(&evalAtQP__), ALIGNMENT, tensor::evalAtQP::size()*sizeof(real));
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    fillWithStuff(dofsQP__, tensor::dofsQP::size());
    fillWithStuff(evalAtQP__, tensor::evalAtQP::size());
    fillWithStuff(Q__, tensor::Q::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::evalAtQP::HardwareFlops);
    }
    kernel::evalAtQP _kernel_evalAtQP;
    _kernel_evalAtQP.dofsQP = dofsQP__;
    _kernel_evalAtQP.evalAtQP = evalAtQP__;
    _kernel_evalAtQP.Q = Q__;
    _kernel_evalAtQP.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel_evalAtQP.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::evalAtQP::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::evalAtQP::HardwareFlops) * _reps / _time / 1.0e9;
    printf("evalAtQP,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::evalAtQP::NonZeroFlops, kernel::evalAtQP::HardwareFlops, _nzflops, _flops);
    free(dofsQP__);
    free(evalAtQP__);
    free(Q__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* kDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivM__0), ALIGNMENT, tensor::kDivM::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* star__0;
    _error = posix_memalign(reinterpret_cast<void**>(&star__0), ALIGNMENT, tensor::star::size(0)*sizeof(real));
    real* kDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivM__1), ALIGNMENT, tensor::kDivM::size(1)*sizeof(real));
    real* star__1;
    _error = posix_memalign(reinterpret_cast<void**>(&star__1), ALIGNMENT, tensor::star::size(1)*sizeof(real));
    real* kDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivM__2), ALIGNMENT, tensor::kDivM::size(2)*sizeof(real));
    real* star__2;
    _error = posix_memalign(reinterpret_cast<void**>(&star__2), ALIGNMENT, tensor::star::size(2)*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(kDivM__0, tensor::kDivM::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(star__0, tensor::star::size(0));
    fillWithStuff(kDivM__1, tensor::kDivM::size(1));
    fillWithStuff(star__1, tensor::star::size(1));
    fillWithStuff(kDivM__2, tensor::kDivM::size(2));
    fillWithStuff(star__2, tensor::star::size(2));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::volume::HardwareFlops);
    }
    kernel::volume _kernel_volume;
    _kernel_volume.Q = Q__;
    _kernel_volume.kDivM(0) = kDivM__0;
    _kernel_volume.I = I__;
    _kernel_volume.star(0) = star__0;
    _kernel_volume.kDivM(1) = kDivM__1;
    _kernel_volume.star(1) = star__1;
    _kernel_volume.kDivM(2) = kDivM__2;
    _kernel_volume.star(2) = star__2;
    _kernel_volume.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel_volume.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::volume::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::volume::HardwareFlops) * _reps / _time / 1.0e9;
    printf("volume,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::volume::NonZeroFlops, kernel::volume::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(kDivM__0);
    free(I__);
    free(star__0);
    free(kDivM__1);
    free(star__1);
    free(kDivM__2);
    free(star__2);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fMrT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fMrT__0), ALIGNMENT, tensor::fMrT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AplusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AplusT__), ALIGNMENT, tensor::AplusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fMrT__0, tensor::fMrT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AplusT__, tensor::AplusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_localFlux_0::HardwareFlops);
    }
    kernel::_localFlux_0 _kernel__localFlux_0;
    _kernel__localFlux_0.Q = Q__;
    _kernel__localFlux_0.rDivM(0) = rDivM__0;
    _kernel__localFlux_0.fMrT(0) = fMrT__0;
    _kernel__localFlux_0.I = I__;
    _kernel__localFlux_0.AplusT = AplusT__;
    _kernel__localFlux_0.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__localFlux_0.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_localFlux_0::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_localFlux_0::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_localFlux_0,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_localFlux_0::NonZeroFlops, kernel::_localFlux_0::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fMrT__0);
    free(I__);
    free(AplusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fMrT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fMrT__1), ALIGNMENT, tensor::fMrT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AplusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AplusT__), ALIGNMENT, tensor::AplusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fMrT__1, tensor::fMrT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AplusT__, tensor::AplusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_localFlux_1::HardwareFlops);
    }
    kernel::_localFlux_1 _kernel__localFlux_1;
    _kernel__localFlux_1.Q = Q__;
    _kernel__localFlux_1.rDivM(1) = rDivM__1;
    _kernel__localFlux_1.fMrT(1) = fMrT__1;
    _kernel__localFlux_1.I = I__;
    _kernel__localFlux_1.AplusT = AplusT__;
    _kernel__localFlux_1.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__localFlux_1.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_localFlux_1::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_localFlux_1::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_localFlux_1,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_localFlux_1::NonZeroFlops, kernel::_localFlux_1::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fMrT__1);
    free(I__);
    free(AplusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fMrT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fMrT__2), ALIGNMENT, tensor::fMrT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AplusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AplusT__), ALIGNMENT, tensor::AplusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fMrT__2, tensor::fMrT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AplusT__, tensor::AplusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_localFlux_2::HardwareFlops);
    }
    kernel::_localFlux_2 _kernel__localFlux_2;
    _kernel__localFlux_2.Q = Q__;
    _kernel__localFlux_2.rDivM(2) = rDivM__2;
    _kernel__localFlux_2.fMrT(2) = fMrT__2;
    _kernel__localFlux_2.I = I__;
    _kernel__localFlux_2.AplusT = AplusT__;
    _kernel__localFlux_2.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__localFlux_2.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_localFlux_2::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_localFlux_2::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_localFlux_2,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_localFlux_2::NonZeroFlops, kernel::_localFlux_2::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fMrT__2);
    free(I__);
    free(AplusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fMrT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&fMrT__3), ALIGNMENT, tensor::fMrT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AplusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AplusT__), ALIGNMENT, tensor::AplusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fMrT__3, tensor::fMrT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AplusT__, tensor::AplusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_localFlux_3::HardwareFlops);
    }
    kernel::_localFlux_3 _kernel__localFlux_3;
    _kernel__localFlux_3.Q = Q__;
    _kernel__localFlux_3.rDivM(3) = rDivM__3;
    _kernel__localFlux_3.fMrT(3) = fMrT__3;
    _kernel__localFlux_3.I = I__;
    _kernel__localFlux_3.AplusT = AplusT__;
    _kernel__localFlux_3.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__localFlux_3.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_localFlux_3::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_localFlux_3::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_localFlux_3,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_localFlux_3::NonZeroFlops, kernel::_localFlux_3::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fMrT__3);
    free(I__);
    free(AplusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_0::HardwareFlops);
    }
    kernel::_neighboringFlux_0 _kernel__neighboringFlux_0;
    _kernel__neighboringFlux_0.Q = Q__;
    _kernel__neighboringFlux_0.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_0.fP(0) = fP__0;
    _kernel__neighboringFlux_0.rT(0) = rT__0;
    _kernel__neighboringFlux_0.I = I__;
    _kernel__neighboringFlux_0.AminusT = AminusT__;
    _kernel__neighboringFlux_0.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_0.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_0::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_0::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_0,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_0::NonZeroFlops, kernel::_neighboringFlux_0::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__0);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_12::HardwareFlops);
    }
    kernel::_neighboringFlux_12 _kernel__neighboringFlux_12;
    _kernel__neighboringFlux_12.Q = Q__;
    _kernel__neighboringFlux_12.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_12.fP(0) = fP__0;
    _kernel__neighboringFlux_12.rT(0) = rT__0;
    _kernel__neighboringFlux_12.I = I__;
    _kernel__neighboringFlux_12.AminusT = AminusT__;
    _kernel__neighboringFlux_12.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_12.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_12::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_12::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_12,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_12::NonZeroFlops, kernel::_neighboringFlux_12::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__0);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_24::HardwareFlops);
    }
    kernel::_neighboringFlux_24 _kernel__neighboringFlux_24;
    _kernel__neighboringFlux_24.Q = Q__;
    _kernel__neighboringFlux_24.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_24.fP(0) = fP__0;
    _kernel__neighboringFlux_24.rT(0) = rT__0;
    _kernel__neighboringFlux_24.I = I__;
    _kernel__neighboringFlux_24.AminusT = AminusT__;
    _kernel__neighboringFlux_24.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_24.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_24::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_24::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_24,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_24::NonZeroFlops, kernel::_neighboringFlux_24::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__0);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_36::HardwareFlops);
    }
    kernel::_neighboringFlux_36 _kernel__neighboringFlux_36;
    _kernel__neighboringFlux_36.Q = Q__;
    _kernel__neighboringFlux_36.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_36.fP(0) = fP__0;
    _kernel__neighboringFlux_36.rT(0) = rT__0;
    _kernel__neighboringFlux_36.I = I__;
    _kernel__neighboringFlux_36.AminusT = AminusT__;
    _kernel__neighboringFlux_36.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_36.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_36::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_36::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_36,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_36::NonZeroFlops, kernel::_neighboringFlux_36::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__0);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_3::HardwareFlops);
    }
    kernel::_neighboringFlux_3 _kernel__neighboringFlux_3;
    _kernel__neighboringFlux_3.Q = Q__;
    _kernel__neighboringFlux_3.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_3.fP(0) = fP__0;
    _kernel__neighboringFlux_3.rT(1) = rT__1;
    _kernel__neighboringFlux_3.I = I__;
    _kernel__neighboringFlux_3.AminusT = AminusT__;
    _kernel__neighboringFlux_3.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_3.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_3::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_3::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_3,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_3::NonZeroFlops, kernel::_neighboringFlux_3::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__0);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_15::HardwareFlops);
    }
    kernel::_neighboringFlux_15 _kernel__neighboringFlux_15;
    _kernel__neighboringFlux_15.Q = Q__;
    _kernel__neighboringFlux_15.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_15.fP(0) = fP__0;
    _kernel__neighboringFlux_15.rT(1) = rT__1;
    _kernel__neighboringFlux_15.I = I__;
    _kernel__neighboringFlux_15.AminusT = AminusT__;
    _kernel__neighboringFlux_15.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_15.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_15::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_15::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_15,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_15::NonZeroFlops, kernel::_neighboringFlux_15::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__0);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_27::HardwareFlops);
    }
    kernel::_neighboringFlux_27 _kernel__neighboringFlux_27;
    _kernel__neighboringFlux_27.Q = Q__;
    _kernel__neighboringFlux_27.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_27.fP(0) = fP__0;
    _kernel__neighboringFlux_27.rT(1) = rT__1;
    _kernel__neighboringFlux_27.I = I__;
    _kernel__neighboringFlux_27.AminusT = AminusT__;
    _kernel__neighboringFlux_27.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_27.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_27::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_27::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_27,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_27::NonZeroFlops, kernel::_neighboringFlux_27::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__0);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_39::HardwareFlops);
    }
    kernel::_neighboringFlux_39 _kernel__neighboringFlux_39;
    _kernel__neighboringFlux_39.Q = Q__;
    _kernel__neighboringFlux_39.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_39.fP(0) = fP__0;
    _kernel__neighboringFlux_39.rT(1) = rT__1;
    _kernel__neighboringFlux_39.I = I__;
    _kernel__neighboringFlux_39.AminusT = AminusT__;
    _kernel__neighboringFlux_39.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_39.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_39::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_39::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_39,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_39::NonZeroFlops, kernel::_neighboringFlux_39::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__0);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_6::HardwareFlops);
    }
    kernel::_neighboringFlux_6 _kernel__neighboringFlux_6;
    _kernel__neighboringFlux_6.Q = Q__;
    _kernel__neighboringFlux_6.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_6.fP(0) = fP__0;
    _kernel__neighboringFlux_6.rT(2) = rT__2;
    _kernel__neighboringFlux_6.I = I__;
    _kernel__neighboringFlux_6.AminusT = AminusT__;
    _kernel__neighboringFlux_6.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_6.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_6::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_6::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_6,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_6::NonZeroFlops, kernel::_neighboringFlux_6::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__0);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_18::HardwareFlops);
    }
    kernel::_neighboringFlux_18 _kernel__neighboringFlux_18;
    _kernel__neighboringFlux_18.Q = Q__;
    _kernel__neighboringFlux_18.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_18.fP(0) = fP__0;
    _kernel__neighboringFlux_18.rT(2) = rT__2;
    _kernel__neighboringFlux_18.I = I__;
    _kernel__neighboringFlux_18.AminusT = AminusT__;
    _kernel__neighboringFlux_18.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_18.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_18::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_18::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_18,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_18::NonZeroFlops, kernel::_neighboringFlux_18::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__0);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_30::HardwareFlops);
    }
    kernel::_neighboringFlux_30 _kernel__neighboringFlux_30;
    _kernel__neighboringFlux_30.Q = Q__;
    _kernel__neighboringFlux_30.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_30.fP(0) = fP__0;
    _kernel__neighboringFlux_30.rT(2) = rT__2;
    _kernel__neighboringFlux_30.I = I__;
    _kernel__neighboringFlux_30.AminusT = AminusT__;
    _kernel__neighboringFlux_30.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_30.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_30::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_30::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_30,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_30::NonZeroFlops, kernel::_neighboringFlux_30::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__0);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_42::HardwareFlops);
    }
    kernel::_neighboringFlux_42 _kernel__neighboringFlux_42;
    _kernel__neighboringFlux_42.Q = Q__;
    _kernel__neighboringFlux_42.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_42.fP(0) = fP__0;
    _kernel__neighboringFlux_42.rT(2) = rT__2;
    _kernel__neighboringFlux_42.I = I__;
    _kernel__neighboringFlux_42.AminusT = AminusT__;
    _kernel__neighboringFlux_42.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_42.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_42::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_42::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_42,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_42::NonZeroFlops, kernel::_neighboringFlux_42::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__0);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_9::HardwareFlops);
    }
    kernel::_neighboringFlux_9 _kernel__neighboringFlux_9;
    _kernel__neighboringFlux_9.Q = Q__;
    _kernel__neighboringFlux_9.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_9.fP(0) = fP__0;
    _kernel__neighboringFlux_9.rT(3) = rT__3;
    _kernel__neighboringFlux_9.I = I__;
    _kernel__neighboringFlux_9.AminusT = AminusT__;
    _kernel__neighboringFlux_9.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_9.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_9::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_9::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_9,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_9::NonZeroFlops, kernel::_neighboringFlux_9::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__0);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_21::HardwareFlops);
    }
    kernel::_neighboringFlux_21 _kernel__neighboringFlux_21;
    _kernel__neighboringFlux_21.Q = Q__;
    _kernel__neighboringFlux_21.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_21.fP(0) = fP__0;
    _kernel__neighboringFlux_21.rT(3) = rT__3;
    _kernel__neighboringFlux_21.I = I__;
    _kernel__neighboringFlux_21.AminusT = AminusT__;
    _kernel__neighboringFlux_21.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_21.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_21::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_21::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_21,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_21::NonZeroFlops, kernel::_neighboringFlux_21::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__0);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_33::HardwareFlops);
    }
    kernel::_neighboringFlux_33 _kernel__neighboringFlux_33;
    _kernel__neighboringFlux_33.Q = Q__;
    _kernel__neighboringFlux_33.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_33.fP(0) = fP__0;
    _kernel__neighboringFlux_33.rT(3) = rT__3;
    _kernel__neighboringFlux_33.I = I__;
    _kernel__neighboringFlux_33.AminusT = AminusT__;
    _kernel__neighboringFlux_33.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_33.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_33::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_33::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_33,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_33::NonZeroFlops, kernel::_neighboringFlux_33::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__0);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__0;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__0), ALIGNMENT, tensor::fP::size(0)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__0, tensor::fP::size(0));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_45::HardwareFlops);
    }
    kernel::_neighboringFlux_45 _kernel__neighboringFlux_45;
    _kernel__neighboringFlux_45.Q = Q__;
    _kernel__neighboringFlux_45.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_45.fP(0) = fP__0;
    _kernel__neighboringFlux_45.rT(3) = rT__3;
    _kernel__neighboringFlux_45.I = I__;
    _kernel__neighboringFlux_45.AminusT = AminusT__;
    _kernel__neighboringFlux_45.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_45.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_45::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_45::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_45,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_45::NonZeroFlops, kernel::_neighboringFlux_45::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__0);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_1::HardwareFlops);
    }
    kernel::_neighboringFlux_1 _kernel__neighboringFlux_1;
    _kernel__neighboringFlux_1.Q = Q__;
    _kernel__neighboringFlux_1.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_1.fP(1) = fP__1;
    _kernel__neighboringFlux_1.rT(0) = rT__0;
    _kernel__neighboringFlux_1.I = I__;
    _kernel__neighboringFlux_1.AminusT = AminusT__;
    _kernel__neighboringFlux_1.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_1.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_1::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_1::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_1,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_1::NonZeroFlops, kernel::_neighboringFlux_1::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__1);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_13::HardwareFlops);
    }
    kernel::_neighboringFlux_13 _kernel__neighboringFlux_13;
    _kernel__neighboringFlux_13.Q = Q__;
    _kernel__neighboringFlux_13.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_13.fP(1) = fP__1;
    _kernel__neighboringFlux_13.rT(0) = rT__0;
    _kernel__neighboringFlux_13.I = I__;
    _kernel__neighboringFlux_13.AminusT = AminusT__;
    _kernel__neighboringFlux_13.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_13.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_13::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_13::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_13,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_13::NonZeroFlops, kernel::_neighboringFlux_13::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__1);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_25::HardwareFlops);
    }
    kernel::_neighboringFlux_25 _kernel__neighboringFlux_25;
    _kernel__neighboringFlux_25.Q = Q__;
    _kernel__neighboringFlux_25.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_25.fP(1) = fP__1;
    _kernel__neighboringFlux_25.rT(0) = rT__0;
    _kernel__neighboringFlux_25.I = I__;
    _kernel__neighboringFlux_25.AminusT = AminusT__;
    _kernel__neighboringFlux_25.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_25.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_25::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_25::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_25,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_25::NonZeroFlops, kernel::_neighboringFlux_25::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__1);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_37::HardwareFlops);
    }
    kernel::_neighboringFlux_37 _kernel__neighboringFlux_37;
    _kernel__neighboringFlux_37.Q = Q__;
    _kernel__neighboringFlux_37.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_37.fP(1) = fP__1;
    _kernel__neighboringFlux_37.rT(0) = rT__0;
    _kernel__neighboringFlux_37.I = I__;
    _kernel__neighboringFlux_37.AminusT = AminusT__;
    _kernel__neighboringFlux_37.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_37.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_37::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_37::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_37,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_37::NonZeroFlops, kernel::_neighboringFlux_37::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__1);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_4::HardwareFlops);
    }
    kernel::_neighboringFlux_4 _kernel__neighboringFlux_4;
    _kernel__neighboringFlux_4.Q = Q__;
    _kernel__neighboringFlux_4.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_4.fP(1) = fP__1;
    _kernel__neighboringFlux_4.rT(1) = rT__1;
    _kernel__neighboringFlux_4.I = I__;
    _kernel__neighboringFlux_4.AminusT = AminusT__;
    _kernel__neighboringFlux_4.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_4.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_4::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_4::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_4,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_4::NonZeroFlops, kernel::_neighboringFlux_4::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__1);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_16::HardwareFlops);
    }
    kernel::_neighboringFlux_16 _kernel__neighboringFlux_16;
    _kernel__neighboringFlux_16.Q = Q__;
    _kernel__neighboringFlux_16.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_16.fP(1) = fP__1;
    _kernel__neighboringFlux_16.rT(1) = rT__1;
    _kernel__neighboringFlux_16.I = I__;
    _kernel__neighboringFlux_16.AminusT = AminusT__;
    _kernel__neighboringFlux_16.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_16.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_16::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_16::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_16,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_16::NonZeroFlops, kernel::_neighboringFlux_16::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__1);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_28::HardwareFlops);
    }
    kernel::_neighboringFlux_28 _kernel__neighboringFlux_28;
    _kernel__neighboringFlux_28.Q = Q__;
    _kernel__neighboringFlux_28.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_28.fP(1) = fP__1;
    _kernel__neighboringFlux_28.rT(1) = rT__1;
    _kernel__neighboringFlux_28.I = I__;
    _kernel__neighboringFlux_28.AminusT = AminusT__;
    _kernel__neighboringFlux_28.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_28.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_28::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_28::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_28,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_28::NonZeroFlops, kernel::_neighboringFlux_28::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__1);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_40::HardwareFlops);
    }
    kernel::_neighboringFlux_40 _kernel__neighboringFlux_40;
    _kernel__neighboringFlux_40.Q = Q__;
    _kernel__neighboringFlux_40.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_40.fP(1) = fP__1;
    _kernel__neighboringFlux_40.rT(1) = rT__1;
    _kernel__neighboringFlux_40.I = I__;
    _kernel__neighboringFlux_40.AminusT = AminusT__;
    _kernel__neighboringFlux_40.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_40.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_40::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_40::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_40,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_40::NonZeroFlops, kernel::_neighboringFlux_40::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__1);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_7::HardwareFlops);
    }
    kernel::_neighboringFlux_7 _kernel__neighboringFlux_7;
    _kernel__neighboringFlux_7.Q = Q__;
    _kernel__neighboringFlux_7.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_7.fP(1) = fP__1;
    _kernel__neighboringFlux_7.rT(2) = rT__2;
    _kernel__neighboringFlux_7.I = I__;
    _kernel__neighboringFlux_7.AminusT = AminusT__;
    _kernel__neighboringFlux_7.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_7.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_7::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_7::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_7,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_7::NonZeroFlops, kernel::_neighboringFlux_7::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__1);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_19::HardwareFlops);
    }
    kernel::_neighboringFlux_19 _kernel__neighboringFlux_19;
    _kernel__neighboringFlux_19.Q = Q__;
    _kernel__neighboringFlux_19.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_19.fP(1) = fP__1;
    _kernel__neighboringFlux_19.rT(2) = rT__2;
    _kernel__neighboringFlux_19.I = I__;
    _kernel__neighboringFlux_19.AminusT = AminusT__;
    _kernel__neighboringFlux_19.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_19.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_19::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_19::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_19,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_19::NonZeroFlops, kernel::_neighboringFlux_19::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__1);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_31::HardwareFlops);
    }
    kernel::_neighboringFlux_31 _kernel__neighboringFlux_31;
    _kernel__neighboringFlux_31.Q = Q__;
    _kernel__neighboringFlux_31.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_31.fP(1) = fP__1;
    _kernel__neighboringFlux_31.rT(2) = rT__2;
    _kernel__neighboringFlux_31.I = I__;
    _kernel__neighboringFlux_31.AminusT = AminusT__;
    _kernel__neighboringFlux_31.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_31.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_31::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_31::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_31,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_31::NonZeroFlops, kernel::_neighboringFlux_31::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__1);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_43::HardwareFlops);
    }
    kernel::_neighboringFlux_43 _kernel__neighboringFlux_43;
    _kernel__neighboringFlux_43.Q = Q__;
    _kernel__neighboringFlux_43.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_43.fP(1) = fP__1;
    _kernel__neighboringFlux_43.rT(2) = rT__2;
    _kernel__neighboringFlux_43.I = I__;
    _kernel__neighboringFlux_43.AminusT = AminusT__;
    _kernel__neighboringFlux_43.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_43.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_43::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_43::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_43,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_43::NonZeroFlops, kernel::_neighboringFlux_43::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__1);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_10::HardwareFlops);
    }
    kernel::_neighboringFlux_10 _kernel__neighboringFlux_10;
    _kernel__neighboringFlux_10.Q = Q__;
    _kernel__neighboringFlux_10.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_10.fP(1) = fP__1;
    _kernel__neighboringFlux_10.rT(3) = rT__3;
    _kernel__neighboringFlux_10.I = I__;
    _kernel__neighboringFlux_10.AminusT = AminusT__;
    _kernel__neighboringFlux_10.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_10.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_10::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_10::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_10,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_10::NonZeroFlops, kernel::_neighboringFlux_10::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__1);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_22::HardwareFlops);
    }
    kernel::_neighboringFlux_22 _kernel__neighboringFlux_22;
    _kernel__neighboringFlux_22.Q = Q__;
    _kernel__neighboringFlux_22.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_22.fP(1) = fP__1;
    _kernel__neighboringFlux_22.rT(3) = rT__3;
    _kernel__neighboringFlux_22.I = I__;
    _kernel__neighboringFlux_22.AminusT = AminusT__;
    _kernel__neighboringFlux_22.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_22.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_22::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_22::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_22,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_22::NonZeroFlops, kernel::_neighboringFlux_22::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__1);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_34::HardwareFlops);
    }
    kernel::_neighboringFlux_34 _kernel__neighboringFlux_34;
    _kernel__neighboringFlux_34.Q = Q__;
    _kernel__neighboringFlux_34.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_34.fP(1) = fP__1;
    _kernel__neighboringFlux_34.rT(3) = rT__3;
    _kernel__neighboringFlux_34.I = I__;
    _kernel__neighboringFlux_34.AminusT = AminusT__;
    _kernel__neighboringFlux_34.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_34.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_34::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_34::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_34,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_34::NonZeroFlops, kernel::_neighboringFlux_34::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__1);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__1;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__1), ALIGNMENT, tensor::fP::size(1)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__1, tensor::fP::size(1));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_46::HardwareFlops);
    }
    kernel::_neighboringFlux_46 _kernel__neighboringFlux_46;
    _kernel__neighboringFlux_46.Q = Q__;
    _kernel__neighboringFlux_46.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_46.fP(1) = fP__1;
    _kernel__neighboringFlux_46.rT(3) = rT__3;
    _kernel__neighboringFlux_46.I = I__;
    _kernel__neighboringFlux_46.AminusT = AminusT__;
    _kernel__neighboringFlux_46.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_46.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_46::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_46::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_46,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_46::NonZeroFlops, kernel::_neighboringFlux_46::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__1);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_2::HardwareFlops);
    }
    kernel::_neighboringFlux_2 _kernel__neighboringFlux_2;
    _kernel__neighboringFlux_2.Q = Q__;
    _kernel__neighboringFlux_2.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_2.fP(2) = fP__2;
    _kernel__neighboringFlux_2.rT(0) = rT__0;
    _kernel__neighboringFlux_2.I = I__;
    _kernel__neighboringFlux_2.AminusT = AminusT__;
    _kernel__neighboringFlux_2.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_2.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_2::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_2::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_2,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_2::NonZeroFlops, kernel::_neighboringFlux_2::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__2);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_14::HardwareFlops);
    }
    kernel::_neighboringFlux_14 _kernel__neighboringFlux_14;
    _kernel__neighboringFlux_14.Q = Q__;
    _kernel__neighboringFlux_14.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_14.fP(2) = fP__2;
    _kernel__neighboringFlux_14.rT(0) = rT__0;
    _kernel__neighboringFlux_14.I = I__;
    _kernel__neighboringFlux_14.AminusT = AminusT__;
    _kernel__neighboringFlux_14.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_14.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_14::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_14::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_14,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_14::NonZeroFlops, kernel::_neighboringFlux_14::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__2);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_26::HardwareFlops);
    }
    kernel::_neighboringFlux_26 _kernel__neighboringFlux_26;
    _kernel__neighboringFlux_26.Q = Q__;
    _kernel__neighboringFlux_26.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_26.fP(2) = fP__2;
    _kernel__neighboringFlux_26.rT(0) = rT__0;
    _kernel__neighboringFlux_26.I = I__;
    _kernel__neighboringFlux_26.AminusT = AminusT__;
    _kernel__neighboringFlux_26.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_26.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_26::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_26::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_26,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_26::NonZeroFlops, kernel::_neighboringFlux_26::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__2);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__0), ALIGNMENT, tensor::rT::size(0)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__0, tensor::rT::size(0));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_38::HardwareFlops);
    }
    kernel::_neighboringFlux_38 _kernel__neighboringFlux_38;
    _kernel__neighboringFlux_38.Q = Q__;
    _kernel__neighboringFlux_38.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_38.fP(2) = fP__2;
    _kernel__neighboringFlux_38.rT(0) = rT__0;
    _kernel__neighboringFlux_38.I = I__;
    _kernel__neighboringFlux_38.AminusT = AminusT__;
    _kernel__neighboringFlux_38.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_38.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_38::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_38::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_38,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_38::NonZeroFlops, kernel::_neighboringFlux_38::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__2);
    free(rT__0);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_5::HardwareFlops);
    }
    kernel::_neighboringFlux_5 _kernel__neighboringFlux_5;
    _kernel__neighboringFlux_5.Q = Q__;
    _kernel__neighboringFlux_5.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_5.fP(2) = fP__2;
    _kernel__neighboringFlux_5.rT(1) = rT__1;
    _kernel__neighboringFlux_5.I = I__;
    _kernel__neighboringFlux_5.AminusT = AminusT__;
    _kernel__neighboringFlux_5.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_5.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_5::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_5::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_5,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_5::NonZeroFlops, kernel::_neighboringFlux_5::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__2);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_17::HardwareFlops);
    }
    kernel::_neighboringFlux_17 _kernel__neighboringFlux_17;
    _kernel__neighboringFlux_17.Q = Q__;
    _kernel__neighboringFlux_17.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_17.fP(2) = fP__2;
    _kernel__neighboringFlux_17.rT(1) = rT__1;
    _kernel__neighboringFlux_17.I = I__;
    _kernel__neighboringFlux_17.AminusT = AminusT__;
    _kernel__neighboringFlux_17.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_17.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_17::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_17::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_17,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_17::NonZeroFlops, kernel::_neighboringFlux_17::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__2);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_29::HardwareFlops);
    }
    kernel::_neighboringFlux_29 _kernel__neighboringFlux_29;
    _kernel__neighboringFlux_29.Q = Q__;
    _kernel__neighboringFlux_29.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_29.fP(2) = fP__2;
    _kernel__neighboringFlux_29.rT(1) = rT__1;
    _kernel__neighboringFlux_29.I = I__;
    _kernel__neighboringFlux_29.AminusT = AminusT__;
    _kernel__neighboringFlux_29.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_29.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_29::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_29::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_29,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_29::NonZeroFlops, kernel::_neighboringFlux_29::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__2);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__1), ALIGNMENT, tensor::rT::size(1)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__1, tensor::rT::size(1));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_41::HardwareFlops);
    }
    kernel::_neighboringFlux_41 _kernel__neighboringFlux_41;
    _kernel__neighboringFlux_41.Q = Q__;
    _kernel__neighboringFlux_41.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_41.fP(2) = fP__2;
    _kernel__neighboringFlux_41.rT(1) = rT__1;
    _kernel__neighboringFlux_41.I = I__;
    _kernel__neighboringFlux_41.AminusT = AminusT__;
    _kernel__neighboringFlux_41.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_41.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_41::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_41::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_41,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_41::NonZeroFlops, kernel::_neighboringFlux_41::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__2);
    free(rT__1);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_8::HardwareFlops);
    }
    kernel::_neighboringFlux_8 _kernel__neighboringFlux_8;
    _kernel__neighboringFlux_8.Q = Q__;
    _kernel__neighboringFlux_8.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_8.fP(2) = fP__2;
    _kernel__neighboringFlux_8.rT(2) = rT__2;
    _kernel__neighboringFlux_8.I = I__;
    _kernel__neighboringFlux_8.AminusT = AminusT__;
    _kernel__neighboringFlux_8.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_8.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_8::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_8::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_8,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_8::NonZeroFlops, kernel::_neighboringFlux_8::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__2);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_20::HardwareFlops);
    }
    kernel::_neighboringFlux_20 _kernel__neighboringFlux_20;
    _kernel__neighboringFlux_20.Q = Q__;
    _kernel__neighboringFlux_20.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_20.fP(2) = fP__2;
    _kernel__neighboringFlux_20.rT(2) = rT__2;
    _kernel__neighboringFlux_20.I = I__;
    _kernel__neighboringFlux_20.AminusT = AminusT__;
    _kernel__neighboringFlux_20.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_20.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_20::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_20::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_20,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_20::NonZeroFlops, kernel::_neighboringFlux_20::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__2);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_32::HardwareFlops);
    }
    kernel::_neighboringFlux_32 _kernel__neighboringFlux_32;
    _kernel__neighboringFlux_32.Q = Q__;
    _kernel__neighboringFlux_32.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_32.fP(2) = fP__2;
    _kernel__neighboringFlux_32.rT(2) = rT__2;
    _kernel__neighboringFlux_32.I = I__;
    _kernel__neighboringFlux_32.AminusT = AminusT__;
    _kernel__neighboringFlux_32.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_32.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_32::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_32::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_32,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_32::NonZeroFlops, kernel::_neighboringFlux_32::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__2);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__2), ALIGNMENT, tensor::rT::size(2)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__2, tensor::rT::size(2));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_44::HardwareFlops);
    }
    kernel::_neighboringFlux_44 _kernel__neighboringFlux_44;
    _kernel__neighboringFlux_44.Q = Q__;
    _kernel__neighboringFlux_44.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_44.fP(2) = fP__2;
    _kernel__neighboringFlux_44.rT(2) = rT__2;
    _kernel__neighboringFlux_44.I = I__;
    _kernel__neighboringFlux_44.AminusT = AminusT__;
    _kernel__neighboringFlux_44.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_44.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_44::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_44::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_44,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_44::NonZeroFlops, kernel::_neighboringFlux_44::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__2);
    free(rT__2);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__0;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__0), ALIGNMENT, tensor::rDivM::size(0)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__0, tensor::rDivM::size(0));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_11::HardwareFlops);
    }
    kernel::_neighboringFlux_11 _kernel__neighboringFlux_11;
    _kernel__neighboringFlux_11.Q = Q__;
    _kernel__neighboringFlux_11.rDivM(0) = rDivM__0;
    _kernel__neighboringFlux_11.fP(2) = fP__2;
    _kernel__neighboringFlux_11.rT(3) = rT__3;
    _kernel__neighboringFlux_11.I = I__;
    _kernel__neighboringFlux_11.AminusT = AminusT__;
    _kernel__neighboringFlux_11.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_11.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_11::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_11::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_11,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_11::NonZeroFlops, kernel::_neighboringFlux_11::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__0);
    free(fP__2);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__1;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__1), ALIGNMENT, tensor::rDivM::size(1)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__1, tensor::rDivM::size(1));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_23::HardwareFlops);
    }
    kernel::_neighboringFlux_23 _kernel__neighboringFlux_23;
    _kernel__neighboringFlux_23.Q = Q__;
    _kernel__neighboringFlux_23.rDivM(1) = rDivM__1;
    _kernel__neighboringFlux_23.fP(2) = fP__2;
    _kernel__neighboringFlux_23.rT(3) = rT__3;
    _kernel__neighboringFlux_23.I = I__;
    _kernel__neighboringFlux_23.AminusT = AminusT__;
    _kernel__neighboringFlux_23.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_23.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_23::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_23::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_23,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_23::NonZeroFlops, kernel::_neighboringFlux_23::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__1);
    free(fP__2);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__2;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__2), ALIGNMENT, tensor::rDivM::size(2)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__2, tensor::rDivM::size(2));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_35::HardwareFlops);
    }
    kernel::_neighboringFlux_35 _kernel__neighboringFlux_35;
    _kernel__neighboringFlux_35.Q = Q__;
    _kernel__neighboringFlux_35.rDivM(2) = rDivM__2;
    _kernel__neighboringFlux_35.fP(2) = fP__2;
    _kernel__neighboringFlux_35.rT(3) = rT__3;
    _kernel__neighboringFlux_35.I = I__;
    _kernel__neighboringFlux_35.AminusT = AminusT__;
    _kernel__neighboringFlux_35.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_35.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_35::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_35::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_35,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_35::NonZeroFlops, kernel::_neighboringFlux_35::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__2);
    free(fP__2);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* Q__;
    _error = posix_memalign(reinterpret_cast<void**>(&Q__), ALIGNMENT, tensor::Q::size()*sizeof(real));
    real* rDivM__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rDivM__3), ALIGNMENT, tensor::rDivM::size(3)*sizeof(real));
    real* fP__2;
    _error = posix_memalign(reinterpret_cast<void**>(&fP__2), ALIGNMENT, tensor::fP::size(2)*sizeof(real));
    real* rT__3;
    _error = posix_memalign(reinterpret_cast<void**>(&rT__3), ALIGNMENT, tensor::rT::size(3)*sizeof(real));
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* AminusT__;
    _error = posix_memalign(reinterpret_cast<void**>(&AminusT__), ALIGNMENT, tensor::AminusT::size()*sizeof(real));
    fillWithStuff(Q__, tensor::Q::size());
    fillWithStuff(rDivM__3, tensor::rDivM::size(3));
    fillWithStuff(fP__2, tensor::fP::size(2));
    fillWithStuff(rT__3, tensor::rT::size(3));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(AminusT__, tensor::AminusT::size());
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_neighboringFlux_47::HardwareFlops);
    }
    kernel::_neighboringFlux_47 _kernel__neighboringFlux_47;
    _kernel__neighboringFlux_47.Q = Q__;
    _kernel__neighboringFlux_47.rDivM(3) = rDivM__3;
    _kernel__neighboringFlux_47.fP(2) = fP__2;
    _kernel__neighboringFlux_47.rT(3) = rT__3;
    _kernel__neighboringFlux_47.I = I__;
    _kernel__neighboringFlux_47.AminusT = AminusT__;
    _kernel__neighboringFlux_47.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__neighboringFlux_47.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_neighboringFlux_47::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_neighboringFlux_47::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_neighboringFlux_47,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_neighboringFlux_47::NonZeroFlops, kernel::_neighboringFlux_47::HardwareFlops, _nzflops, _flops);
    free(Q__);
    free(rDivM__3);
    free(fP__2);
    free(rT__3);
    free(I__);
    free(AminusT__);
  }
  {
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* dQ__0;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__0), ALIGNMENT, tensor::dQ::size(0)*sizeof(real));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(dQ__0, tensor::dQ::size(0));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivativeTaylorExpansion_0::HardwareFlops);
    }
    kernel::_derivativeTaylorExpansion_0 _kernel__derivativeTaylorExpansion_0;
    _kernel__derivativeTaylorExpansion_0.I = I__;
    _kernel__derivativeTaylorExpansion_0.dQ(0) = dQ__0;
    _kernel__derivativeTaylorExpansion_0.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivativeTaylorExpansion_0.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivativeTaylorExpansion_0::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivativeTaylorExpansion_0::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivativeTaylorExpansion_0,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivativeTaylorExpansion_0::NonZeroFlops, kernel::_derivativeTaylorExpansion_0::HardwareFlops, _nzflops, _flops);
    free(I__);
    free(dQ__0);
  }
  {
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* dQ__1;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__1), ALIGNMENT, tensor::dQ::size(1)*sizeof(real));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(dQ__1, tensor::dQ::size(1));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivativeTaylorExpansion_1::HardwareFlops);
    }
    kernel::_derivativeTaylorExpansion_1 _kernel__derivativeTaylorExpansion_1;
    _kernel__derivativeTaylorExpansion_1.I = I__;
    _kernel__derivativeTaylorExpansion_1.dQ(1) = dQ__1;
    _kernel__derivativeTaylorExpansion_1.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivativeTaylorExpansion_1.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivativeTaylorExpansion_1::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivativeTaylorExpansion_1::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivativeTaylorExpansion_1,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivativeTaylorExpansion_1::NonZeroFlops, kernel::_derivativeTaylorExpansion_1::HardwareFlops, _nzflops, _flops);
    free(I__);
    free(dQ__1);
  }
  {
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* dQ__2;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__2), ALIGNMENT, tensor::dQ::size(2)*sizeof(real));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(dQ__2, tensor::dQ::size(2));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivativeTaylorExpansion_2::HardwareFlops);
    }
    kernel::_derivativeTaylorExpansion_2 _kernel__derivativeTaylorExpansion_2;
    _kernel__derivativeTaylorExpansion_2.I = I__;
    _kernel__derivativeTaylorExpansion_2.dQ(2) = dQ__2;
    _kernel__derivativeTaylorExpansion_2.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivativeTaylorExpansion_2.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivativeTaylorExpansion_2::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivativeTaylorExpansion_2::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivativeTaylorExpansion_2,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivativeTaylorExpansion_2::NonZeroFlops, kernel::_derivativeTaylorExpansion_2::HardwareFlops, _nzflops, _flops);
    free(I__);
    free(dQ__2);
  }
  {
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* dQ__3;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__3), ALIGNMENT, tensor::dQ::size(3)*sizeof(real));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(dQ__3, tensor::dQ::size(3));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivativeTaylorExpansion_3::HardwareFlops);
    }
    kernel::_derivativeTaylorExpansion_3 _kernel__derivativeTaylorExpansion_3;
    _kernel__derivativeTaylorExpansion_3.I = I__;
    _kernel__derivativeTaylorExpansion_3.dQ(3) = dQ__3;
    _kernel__derivativeTaylorExpansion_3.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivativeTaylorExpansion_3.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivativeTaylorExpansion_3::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivativeTaylorExpansion_3::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivativeTaylorExpansion_3,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivativeTaylorExpansion_3::NonZeroFlops, kernel::_derivativeTaylorExpansion_3::HardwareFlops, _nzflops, _flops);
    free(I__);
    free(dQ__3);
  }
  {
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* dQ__4;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__4), ALIGNMENT, tensor::dQ::size(4)*sizeof(real));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(dQ__4, tensor::dQ::size(4));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivativeTaylorExpansion_4::HardwareFlops);
    }
    kernel::_derivativeTaylorExpansion_4 _kernel__derivativeTaylorExpansion_4;
    _kernel__derivativeTaylorExpansion_4.I = I__;
    _kernel__derivativeTaylorExpansion_4.dQ(4) = dQ__4;
    _kernel__derivativeTaylorExpansion_4.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivativeTaylorExpansion_4.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivativeTaylorExpansion_4::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivativeTaylorExpansion_4::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivativeTaylorExpansion_4,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivativeTaylorExpansion_4::NonZeroFlops, kernel::_derivativeTaylorExpansion_4::HardwareFlops, _nzflops, _flops);
    free(I__);
    free(dQ__4);
  }
  {
    real* I__;
    _error = posix_memalign(reinterpret_cast<void**>(&I__), ALIGNMENT, tensor::I::size()*sizeof(real));
    real* dQ__5;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__5), ALIGNMENT, tensor::dQ::size(5)*sizeof(real));
    fillWithStuff(I__, tensor::I::size());
    fillWithStuff(dQ__5, tensor::dQ::size(5));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivativeTaylorExpansion_5::HardwareFlops);
    }
    kernel::_derivativeTaylorExpansion_5 _kernel__derivativeTaylorExpansion_5;
    _kernel__derivativeTaylorExpansion_5.I = I__;
    _kernel__derivativeTaylorExpansion_5.dQ(5) = dQ__5;
    _kernel__derivativeTaylorExpansion_5.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivativeTaylorExpansion_5.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivativeTaylorExpansion_5::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivativeTaylorExpansion_5::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivativeTaylorExpansion_5,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivativeTaylorExpansion_5::NonZeroFlops, kernel::_derivativeTaylorExpansion_5::HardwareFlops, _nzflops, _flops);
    free(I__);
    free(dQ__5);
  }
  {
    real* dQ__1;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__1), ALIGNMENT, tensor::dQ::size(1)*sizeof(real));
    real* kDivMT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__0), ALIGNMENT, tensor::kDivMT::size(0)*sizeof(real));
    real* dQ__0;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__0), ALIGNMENT, tensor::dQ::size(0)*sizeof(real));
    real* star__0;
    _error = posix_memalign(reinterpret_cast<void**>(&star__0), ALIGNMENT, tensor::star::size(0)*sizeof(real));
    real* kDivMT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__1), ALIGNMENT, tensor::kDivMT::size(1)*sizeof(real));
    real* star__1;
    _error = posix_memalign(reinterpret_cast<void**>(&star__1), ALIGNMENT, tensor::star::size(1)*sizeof(real));
    real* kDivMT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__2), ALIGNMENT, tensor::kDivMT::size(2)*sizeof(real));
    real* star__2;
    _error = posix_memalign(reinterpret_cast<void**>(&star__2), ALIGNMENT, tensor::star::size(2)*sizeof(real));
    fillWithStuff(dQ__1, tensor::dQ::size(1));
    fillWithStuff(kDivMT__0, tensor::kDivMT::size(0));
    fillWithStuff(dQ__0, tensor::dQ::size(0));
    fillWithStuff(star__0, tensor::star::size(0));
    fillWithStuff(kDivMT__1, tensor::kDivMT::size(1));
    fillWithStuff(star__1, tensor::star::size(1));
    fillWithStuff(kDivMT__2, tensor::kDivMT::size(2));
    fillWithStuff(star__2, tensor::star::size(2));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivative_1::HardwareFlops);
    }
    kernel::_derivative_1 _kernel__derivative_1;
    _kernel__derivative_1.dQ(1) = dQ__1;
    _kernel__derivative_1.kDivMT(0) = kDivMT__0;
    _kernel__derivative_1.dQ(0) = dQ__0;
    _kernel__derivative_1.star(0) = star__0;
    _kernel__derivative_1.kDivMT(1) = kDivMT__1;
    _kernel__derivative_1.star(1) = star__1;
    _kernel__derivative_1.kDivMT(2) = kDivMT__2;
    _kernel__derivative_1.star(2) = star__2;
    _kernel__derivative_1.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivative_1.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivative_1::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivative_1::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivative_1,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivative_1::NonZeroFlops, kernel::_derivative_1::HardwareFlops, _nzflops, _flops);
    free(dQ__1);
    free(kDivMT__0);
    free(dQ__0);
    free(star__0);
    free(kDivMT__1);
    free(star__1);
    free(kDivMT__2);
    free(star__2);
  }
  {
    real* dQ__2;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__2), ALIGNMENT, tensor::dQ::size(2)*sizeof(real));
    real* kDivMT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__0), ALIGNMENT, tensor::kDivMT::size(0)*sizeof(real));
    real* dQ__1;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__1), ALIGNMENT, tensor::dQ::size(1)*sizeof(real));
    real* star__0;
    _error = posix_memalign(reinterpret_cast<void**>(&star__0), ALIGNMENT, tensor::star::size(0)*sizeof(real));
    real* kDivMT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__1), ALIGNMENT, tensor::kDivMT::size(1)*sizeof(real));
    real* star__1;
    _error = posix_memalign(reinterpret_cast<void**>(&star__1), ALIGNMENT, tensor::star::size(1)*sizeof(real));
    real* kDivMT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__2), ALIGNMENT, tensor::kDivMT::size(2)*sizeof(real));
    real* star__2;
    _error = posix_memalign(reinterpret_cast<void**>(&star__2), ALIGNMENT, tensor::star::size(2)*sizeof(real));
    fillWithStuff(dQ__2, tensor::dQ::size(2));
    fillWithStuff(kDivMT__0, tensor::kDivMT::size(0));
    fillWithStuff(dQ__1, tensor::dQ::size(1));
    fillWithStuff(star__0, tensor::star::size(0));
    fillWithStuff(kDivMT__1, tensor::kDivMT::size(1));
    fillWithStuff(star__1, tensor::star::size(1));
    fillWithStuff(kDivMT__2, tensor::kDivMT::size(2));
    fillWithStuff(star__2, tensor::star::size(2));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivative_2::HardwareFlops);
    }
    kernel::_derivative_2 _kernel__derivative_2;
    _kernel__derivative_2.dQ(2) = dQ__2;
    _kernel__derivative_2.kDivMT(0) = kDivMT__0;
    _kernel__derivative_2.dQ(1) = dQ__1;
    _kernel__derivative_2.star(0) = star__0;
    _kernel__derivative_2.kDivMT(1) = kDivMT__1;
    _kernel__derivative_2.star(1) = star__1;
    _kernel__derivative_2.kDivMT(2) = kDivMT__2;
    _kernel__derivative_2.star(2) = star__2;
    _kernel__derivative_2.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivative_2.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivative_2::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivative_2::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivative_2,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivative_2::NonZeroFlops, kernel::_derivative_2::HardwareFlops, _nzflops, _flops);
    free(dQ__2);
    free(kDivMT__0);
    free(dQ__1);
    free(star__0);
    free(kDivMT__1);
    free(star__1);
    free(kDivMT__2);
    free(star__2);
  }
  {
    real* dQ__3;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__3), ALIGNMENT, tensor::dQ::size(3)*sizeof(real));
    real* kDivMT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__0), ALIGNMENT, tensor::kDivMT::size(0)*sizeof(real));
    real* dQ__2;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__2), ALIGNMENT, tensor::dQ::size(2)*sizeof(real));
    real* star__0;
    _error = posix_memalign(reinterpret_cast<void**>(&star__0), ALIGNMENT, tensor::star::size(0)*sizeof(real));
    real* kDivMT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__1), ALIGNMENT, tensor::kDivMT::size(1)*sizeof(real));
    real* star__1;
    _error = posix_memalign(reinterpret_cast<void**>(&star__1), ALIGNMENT, tensor::star::size(1)*sizeof(real));
    real* kDivMT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__2), ALIGNMENT, tensor::kDivMT::size(2)*sizeof(real));
    real* star__2;
    _error = posix_memalign(reinterpret_cast<void**>(&star__2), ALIGNMENT, tensor::star::size(2)*sizeof(real));
    fillWithStuff(dQ__3, tensor::dQ::size(3));
    fillWithStuff(kDivMT__0, tensor::kDivMT::size(0));
    fillWithStuff(dQ__2, tensor::dQ::size(2));
    fillWithStuff(star__0, tensor::star::size(0));
    fillWithStuff(kDivMT__1, tensor::kDivMT::size(1));
    fillWithStuff(star__1, tensor::star::size(1));
    fillWithStuff(kDivMT__2, tensor::kDivMT::size(2));
    fillWithStuff(star__2, tensor::star::size(2));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivative_3::HardwareFlops);
    }
    kernel::_derivative_3 _kernel__derivative_3;
    _kernel__derivative_3.dQ(3) = dQ__3;
    _kernel__derivative_3.kDivMT(0) = kDivMT__0;
    _kernel__derivative_3.dQ(2) = dQ__2;
    _kernel__derivative_3.star(0) = star__0;
    _kernel__derivative_3.kDivMT(1) = kDivMT__1;
    _kernel__derivative_3.star(1) = star__1;
    _kernel__derivative_3.kDivMT(2) = kDivMT__2;
    _kernel__derivative_3.star(2) = star__2;
    _kernel__derivative_3.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivative_3.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivative_3::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivative_3::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivative_3,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivative_3::NonZeroFlops, kernel::_derivative_3::HardwareFlops, _nzflops, _flops);
    free(dQ__3);
    free(kDivMT__0);
    free(dQ__2);
    free(star__0);
    free(kDivMT__1);
    free(star__1);
    free(kDivMT__2);
    free(star__2);
  }
  {
    real* dQ__4;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__4), ALIGNMENT, tensor::dQ::size(4)*sizeof(real));
    real* kDivMT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__0), ALIGNMENT, tensor::kDivMT::size(0)*sizeof(real));
    real* dQ__3;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__3), ALIGNMENT, tensor::dQ::size(3)*sizeof(real));
    real* star__0;
    _error = posix_memalign(reinterpret_cast<void**>(&star__0), ALIGNMENT, tensor::star::size(0)*sizeof(real));
    real* kDivMT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__1), ALIGNMENT, tensor::kDivMT::size(1)*sizeof(real));
    real* star__1;
    _error = posix_memalign(reinterpret_cast<void**>(&star__1), ALIGNMENT, tensor::star::size(1)*sizeof(real));
    real* kDivMT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__2), ALIGNMENT, tensor::kDivMT::size(2)*sizeof(real));
    real* star__2;
    _error = posix_memalign(reinterpret_cast<void**>(&star__2), ALIGNMENT, tensor::star::size(2)*sizeof(real));
    fillWithStuff(dQ__4, tensor::dQ::size(4));
    fillWithStuff(kDivMT__0, tensor::kDivMT::size(0));
    fillWithStuff(dQ__3, tensor::dQ::size(3));
    fillWithStuff(star__0, tensor::star::size(0));
    fillWithStuff(kDivMT__1, tensor::kDivMT::size(1));
    fillWithStuff(star__1, tensor::star::size(1));
    fillWithStuff(kDivMT__2, tensor::kDivMT::size(2));
    fillWithStuff(star__2, tensor::star::size(2));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivative_4::HardwareFlops);
    }
    kernel::_derivative_4 _kernel__derivative_4;
    _kernel__derivative_4.dQ(4) = dQ__4;
    _kernel__derivative_4.kDivMT(0) = kDivMT__0;
    _kernel__derivative_4.dQ(3) = dQ__3;
    _kernel__derivative_4.star(0) = star__0;
    _kernel__derivative_4.kDivMT(1) = kDivMT__1;
    _kernel__derivative_4.star(1) = star__1;
    _kernel__derivative_4.kDivMT(2) = kDivMT__2;
    _kernel__derivative_4.star(2) = star__2;
    _kernel__derivative_4.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivative_4.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivative_4::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivative_4::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivative_4,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivative_4::NonZeroFlops, kernel::_derivative_4::HardwareFlops, _nzflops, _flops);
    free(dQ__4);
    free(kDivMT__0);
    free(dQ__3);
    free(star__0);
    free(kDivMT__1);
    free(star__1);
    free(kDivMT__2);
    free(star__2);
  }
  {
    real* dQ__5;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__5), ALIGNMENT, tensor::dQ::size(5)*sizeof(real));
    real* kDivMT__0;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__0), ALIGNMENT, tensor::kDivMT::size(0)*sizeof(real));
    real* dQ__4;
    _error = posix_memalign(reinterpret_cast<void**>(&dQ__4), ALIGNMENT, tensor::dQ::size(4)*sizeof(real));
    real* star__0;
    _error = posix_memalign(reinterpret_cast<void**>(&star__0), ALIGNMENT, tensor::star::size(0)*sizeof(real));
    real* kDivMT__1;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__1), ALIGNMENT, tensor::kDivMT::size(1)*sizeof(real));
    real* star__1;
    _error = posix_memalign(reinterpret_cast<void**>(&star__1), ALIGNMENT, tensor::star::size(1)*sizeof(real));
    real* kDivMT__2;
    _error = posix_memalign(reinterpret_cast<void**>(&kDivMT__2), ALIGNMENT, tensor::kDivMT::size(2)*sizeof(real));
    real* star__2;
    _error = posix_memalign(reinterpret_cast<void**>(&star__2), ALIGNMENT, tensor::star::size(2)*sizeof(real));
    fillWithStuff(dQ__5, tensor::dQ::size(5));
    fillWithStuff(kDivMT__0, tensor::kDivMT::size(0));
    fillWithStuff(dQ__4, tensor::dQ::size(4));
    fillWithStuff(star__0, tensor::star::size(0));
    fillWithStuff(kDivMT__1, tensor::kDivMT::size(1));
    fillWithStuff(star__1, tensor::star::size(1));
    fillWithStuff(kDivMT__2, tensor::kDivMT::size(2));
    fillWithStuff(star__2, tensor::star::size(2));
    _reps = _fixedReps;
    if (_reps < 0) {
      _reps = ceil(40000000000.0/kernel::_derivative_5::HardwareFlops);
    }
    kernel::_derivative_5 _kernel__derivative_5;
    _kernel__derivative_5.dQ(5) = dQ__5;
    _kernel__derivative_5.kDivMT(0) = kDivMT__0;
    _kernel__derivative_5.dQ(4) = dQ__4;
    _kernel__derivative_5.star(0) = star__0;
    _kernel__derivative_5.kDivMT(1) = kDivMT__1;
    _kernel__derivative_5.star(1) = star__1;
    _kernel__derivative_5.kDivMT(2) = kDivMT__2;
    _kernel__derivative_5.star(2) = star__2;
    _kernel__derivative_5.execute();
    _sw.start();
    for (int i = 0; i < _reps; ++i) {
      _kernel__derivative_5.execute();
    }
    _time = _sw.stop();
    _nzflops = static_cast<double>(kernel::_derivative_5::NonZeroFlops) * _reps / _time / 1.0e9;
    _flops = static_cast<double>(kernel::_derivative_5::HardwareFlops) * _reps / _time / 1.0e9;
    printf("_derivative_5,%u,%lf,%lu,%lu,%lf,%lf\n", _reps, _time, kernel::_derivative_5::NonZeroFlops, kernel::_derivative_5::HardwareFlops, _nzflops, _flops);
    free(dQ__5);
    free(kDivMT__0);
    free(dQ__4);
    free(star__0);
    free(kDivMT__1);
    free(star__1);
    free(kDivMT__2);
    free(star__2);
  }
  return 0;
}
