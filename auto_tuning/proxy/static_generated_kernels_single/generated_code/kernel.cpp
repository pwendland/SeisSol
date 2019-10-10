#include "mkl_cblas.h"
#include "device_utils.h"
#include <cassert>
#include <cstring>
#include <cstdlib>
#include "subroutine.h"
#include "kernel.h"
namespace seissol {
  constexpr unsigned long const kernel::computeFluxSolverLocal::NonZeroFlops;
  constexpr unsigned long const kernel::computeFluxSolverLocal::HardwareFlops;
  void kernel::computeFluxSolverLocal::execute() {
    assert(!std::isnan(fluxScale));
    assert(AplusT != nullptr);
    assert(QgodLocal != nullptr);
    assert(T != nullptr);
    assert(Tinv != nullptr);
    assert(star(0) != nullptr);
    float *_tmp0, *_tmp1;
    float _buffer0[81] __attribute__((aligned(32)));
    float _buffer1[81] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, 9, 9, 9, 1.0, Tinv, 9, QgodLocal, 9, 0.0, _tmp0, 9);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, 9, 9, 9, 1.0, star(0), 9, T, 9, 0.0, _tmp1, 9);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 9, 9, 9, fluxScale, _tmp0, 9, _tmp1, 9, 0.0, AplusT, 9);
  }
  constexpr unsigned long const kernel::computeFluxSolverNeighbor::NonZeroFlops;
  constexpr unsigned long const kernel::computeFluxSolverNeighbor::HardwareFlops;
  void kernel::computeFluxSolverNeighbor::execute() {
    assert(!std::isnan(fluxScale));
    assert(AminusT != nullptr);
    assert(QgodNeighbor != nullptr);
    assert(T != nullptr);
    assert(Tinv != nullptr);
    assert(star(0) != nullptr);
    float *_tmp0, *_tmp1;
    float _buffer0[81] __attribute__((aligned(32)));
    float _buffer1[81] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, 9, 9, 9, 1.0, Tinv, 9, QgodNeighbor, 9, 0.0, _tmp0, 9);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, 9, 9, 9, 1.0, star(0), 9, T, 9, 0.0, _tmp1, 9);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 9, 9, 9, fluxScale, _tmp0, 9, _tmp1, 9, 0.0, AminusT, 9);
  }
  constexpr unsigned long const kernel::copyQToQFortran::NonZeroFlops;
  constexpr unsigned long const kernel::copyQToQFortran::HardwareFlops;
  void kernel::copyQToQFortran::execute() {
    assert(Q != nullptr);
    assert(QFortran != nullptr);
    for (int b = 0; b < 9; ++b) {
      #pragma omp simd
      for (int a = 0; a < 56; ++a) {
        QFortran[1*a + 56*b] = Q[1*a + 56*b];
      }
    }
  }
  constexpr unsigned long const kernel::projectIniCond::NonZeroFlops;
  constexpr unsigned long const kernel::projectIniCond::HardwareFlops;
  void kernel::projectIniCond::execute() {
    assert(Q != nullptr);
    assert(iniCond != nullptr);
    assert(projectQP != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 343, 1.0, projectQP, 56, iniCond, 344, 0.0, Q, 56);
  }
  constexpr unsigned long const kernel::evalAtQP::NonZeroFlops;
  constexpr unsigned long const kernel::evalAtQP::HardwareFlops;
  void kernel::evalAtQP::execute() {
    assert(Q != nullptr);
    assert(dofsQP != nullptr);
    assert(evalAtQP != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 344, 9, 56, 1.0, evalAtQP, 344, Q, 56, 0.0, dofsQP, 344);
  }
  #ifdef GPU
constexpr unsigned long const kernel::volume::NonZeroFlops;
  constexpr unsigned long const kernel::volume::HardwareFlops;
  void kernel::volume::execute() {
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(kDivM(2) != nullptr);
    assert(kDivM(0) != nullptr);
    assert(kDivM(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    // TODO: allocate all buffers at the main entry point of the program
    // allocating temp memory only on gpu
    float *d_buffer0 = (float*)DeviceScratchMem::get_instance().get_mem(504 * tensor::num_elements_in_cluster );
    _tmp0 = d_buffer0;
    {
    const unsigned jump_to_next_I = tensor::I::jump_to_next;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[0];
    const unsigned jump_to_next__tmp0 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, I, 56, star(0), 9, 0.0, _tmp0, 40, jump_to_next_I, jump_to_next_star, jump_to_next__tmp0, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next_kDivM = tensor::kDivM::jump_to_next[0];
    const unsigned jump_to_next__tmp0 = 504;
    const unsigned jump_to_next_Q = tensor::Q::jump_to_next;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(0), 56, _tmp0, 40, 1.0, Q, 56, jump_to_next_kDivM, jump_to_next__tmp0, jump_to_next_Q, tensor::num_elements_in_cluster);
    }
    _tmp2 = d_buffer0;
    {
    const unsigned jump_to_next_I = tensor::I::jump_to_next;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[1];
    const unsigned jump_to_next__tmp2 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, I, 56, star(1), 9, 0.0, _tmp2, 40, jump_to_next_I, jump_to_next_star, jump_to_next__tmp2, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next_kDivM = tensor::kDivM::jump_to_next[1];
    const unsigned jump_to_next__tmp2 = 504;
    const unsigned jump_to_next_Q = tensor::Q::jump_to_next;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(1), 56, _tmp2, 40, 1.0, Q, 56, jump_to_next_kDivM, jump_to_next__tmp2, jump_to_next_Q, tensor::num_elements_in_cluster);
    }
    _tmp4 = d_buffer0;
    {
    const unsigned jump_to_next_I = tensor::I::jump_to_next;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[2];
    const unsigned jump_to_next__tmp4 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, I, 56, star(2), 9, 0.0, _tmp4, 40, jump_to_next_I, jump_to_next_star, jump_to_next__tmp4, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next_kDivM = tensor::kDivM::jump_to_next[2];
    const unsigned jump_to_next__tmp4 = 504;
    const unsigned jump_to_next_Q = tensor::Q::jump_to_next;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(2), 56, _tmp4, 40, 1.0, Q, 56, jump_to_next_kDivM, jump_to_next__tmp4, jump_to_next_Q, tensor::num_elements_in_cluster);
    }
    DeviceScratchMem::get_instance().free();
  }
  #else
  constexpr unsigned long const kernel::volume::NonZeroFlops;
  constexpr unsigned long const kernel::volume::HardwareFlops;
  void kernel::volume::execute() {
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(kDivM(2) != nullptr);
    assert(kDivM(0) != nullptr);
    assert(kDivM(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    float _buffer0[360] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, I, 56, star(0), 9, 0.0, _tmp0, 40);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(0), 56, _tmp0, 40, 1.0, Q, 56);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, I, 56, star(1), 9, 0.0, _tmp2, 40);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(1), 56, _tmp2, 40, 1.0, Q, 56);
    _tmp4 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, I, 56, star(2), 9, 0.0, _tmp4, 40);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 35, 1.0, kDivM(2), 56, _tmp4, 40, 1.0, Q, 56);
  }
  #endif
  constexpr unsigned long const kernel::rotateGodunovStateLocal::NonZeroFlops;
  constexpr unsigned long const kernel::rotateGodunovStateLocal::HardwareFlops;
  void kernel::rotateGodunovStateLocal::execute() {
    assert(QgodLocal != nullptr);
    assert(Tinv != nullptr);
    assert(godunovMatrix != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, 9, 9, 9, 1.0, Tinv, 9, QgodLocal, 9, 0.0, godunovMatrix, 9);
  }
  constexpr unsigned long const kernel::rotateGodunovStateNeighbor::NonZeroFlops;
  constexpr unsigned long const kernel::rotateGodunovStateNeighbor::HardwareFlops;
  void kernel::rotateGodunovStateNeighbor::execute() {
    assert(QgodNeighbor != nullptr);
    assert(Tinv != nullptr);
    assert(godunovMatrix != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, 9, 9, 9, 1.0, Tinv, 9, QgodNeighbor, 9, 0.0, godunovMatrix, 9);
  }
  constexpr unsigned long const kernel::rotateFluxMatrix::NonZeroFlops;
  constexpr unsigned long const kernel::rotateFluxMatrix::HardwareFlops;
  void kernel::rotateFluxMatrix::execute() {
    assert(!std::isnan(fluxScale));
    assert(T != nullptr);
    assert(fluxSolver != nullptr);
    assert(star(0) != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, 9, 9, 9, fluxScale, star(0), 9, T, 9, 0.0, fluxSolver, 9);
  }
  constexpr unsigned long const kernel::plConvertToNodal::NonZeroFlops;
  constexpr unsigned long const kernel::plConvertToNodal::HardwareFlops;
  void kernel::plConvertToNodal::execute() {
    assert(QStress != nullptr);
    assert(QStressNodal != nullptr);
    assert(initialLoading != nullptr);
    assert(replicateInitialLoading != nullptr);
    assert(v != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 6, 56, 1.0, v, 56, QStress, 56, 0.0, QStressNodal, 56);
    for (int p = 0; p < 6; ++p) {
      #pragma omp simd
      for (int k = 0; k < 56; ++k) {
        QStressNodal[1*k + 56*p] += replicateInitialLoading[1*k] * initialLoading[1*p];
      }
    }
  }
  constexpr unsigned long const kernel::plComputeMean::NonZeroFlops;
  constexpr unsigned long const kernel::plComputeMean::HardwareFlops;
  void kernel::plComputeMean::execute() {
    assert(QStressNodal != nullptr);
    assert(meanStress != nullptr);
    assert(selectBulkAverage != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 1, 3, 1.0, QStressNodal, 56, selectBulkAverage, 3, 0.0, meanStress, 56);
  }
  constexpr unsigned long const kernel::plSubtractMean::NonZeroFlops;
  constexpr unsigned long const kernel::plSubtractMean::HardwareFlops;
  void kernel::plSubtractMean::execute() {
    assert(QStressNodal != nullptr);
    assert(meanStress != nullptr);
    assert(selectBulkNegative != nullptr);
    for (int p = 0; p < 3; ++p) {
      #pragma omp simd
      for (int k = 0; k < 56; ++k) {
        QStressNodal[1*k + 56*p] += meanStress[1*k] * selectBulkNegative[1*p];
      }
    }
  }
  constexpr unsigned long const kernel::plComputeSecondInvariant::NonZeroFlops;
  constexpr unsigned long const kernel::plComputeSecondInvariant::HardwareFlops;
  void kernel::plComputeSecondInvariant::execute() {
    assert(QStressNodal != nullptr);
    assert(secondInvariant != nullptr);
    assert(weightSecondInvariant != nullptr);
    float *_tmp0;
    float _buffer0[336] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    for (int q = 0; q < 6; ++q) {
      #pragma omp simd
      for (int k = 0; k < 56; ++k) {
        _tmp0[1*k + 56*q] = QStressNodal[1*k + 56*q] * QStressNodal[1*k + 56*q];
      }
    }
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 1, 6, 1.0, _tmp0, 56, weightSecondInvariant, 6, 0.0, secondInvariant, 56);
  }
  constexpr unsigned long const kernel::plAdjustStresses::NonZeroFlops;
  constexpr unsigned long const kernel::plAdjustStresses::HardwareFlops;
  void kernel::plAdjustStresses::execute() {
    assert(QStress != nullptr);
    assert(QStressNodal != nullptr);
    assert(vInv != nullptr);
    assert(yieldFactor != nullptr);
    float *_tmp0;
    float _buffer0[336] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    for (int p = 0; p < 6; ++p) {
      #pragma omp simd
      for (int l = 0; l < 56; ++l) {
        _tmp0[1*l + 56*p] = QStressNodal[1*l + 56*p] * yieldFactor[1*l];
      }
    }
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 6, 56, 1.0, vInv, 56, _tmp0, 56, 1.0, QStress, 56);
  }
  constexpr unsigned long const kernel::addVelocity::NonZeroFlops;
  constexpr unsigned long const kernel::addVelocity::HardwareFlops;
  void kernel::addVelocity::execute() {
    assert(I != nullptr);
    assert(displacement != nullptr);
    assert(selectVelocity != nullptr);
    for (int o = 0; o < 56; ++o) {
      displacement[0 + 1*o + 56*0] += 1.0 * I[336 + 1*o + 56*0] * selectVelocity[0 + 0];
      displacement[0 + 1*o + 56*1] += 1.0 * I[336 + 1*o + 56*1] * selectVelocity[0 + 1];
      displacement[0 + 1*o + 56*2] += 1.0 * I[336 + 1*o + 56*2] * selectVelocity[0 + 2];
    }
  }
  constexpr unsigned long const kernel::sourceNRF::NonZeroFlops;
  constexpr unsigned long const kernel::sourceNRF::HardwareFlops;
  void kernel::sourceNRF::execute() {
    assert(Q != nullptr);
    assert(mInvJInvPhisAtSources != nullptr);
    assert(momentNRF != nullptr);
    for (int p = 0; p < 6; ++p) {
      #pragma omp simd
      for (int k = 0; k < 56; ++k) {
        Q[1*k + 56*p] += -1.0 * mInvJInvPhisAtSources[1*k] * momentNRF[1*p];
      }
    }
  }
  constexpr unsigned long const kernel::sourceFSRM::NonZeroFlops;
  constexpr unsigned long const kernel::sourceFSRM::HardwareFlops;
  void kernel::sourceFSRM::execute() {
    assert(!std::isnan(stfIntegral));
    assert(Q != nullptr);
    assert(mInvJInvPhisAtSources != nullptr);
    assert(momentFSRM != nullptr);
    for (int p = 0; p < 9; ++p) {
      #pragma omp simd
      for (int k = 0; k < 56; ++k) {
        Q[1*k + 56*p] += stfIntegral * mInvJInvPhisAtSources[1*k] * momentFSRM[1*p];
      }
    }
  }
  constexpr unsigned long const kernel::evaluateDOFSAtPoint::NonZeroFlops;
  constexpr unsigned long const kernel::evaluateDOFSAtPoint::HardwareFlops;
  void kernel::evaluateDOFSAtPoint::execute() {
    assert(Q != nullptr);
    assert(QAtPoint != nullptr);
    assert(basisFunctions != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, 9, 1, 56, 1.0, Q, 56, basisFunctions, 56, 0.0, QAtPoint, 9);
  }
  #ifdef GPU
 constexpr unsigned long const kernel::localFlux::NonZeroFlops[];
  constexpr unsigned long const kernel::localFlux::HardwareFlops[];
  constexpr kernel::localFlux::member_function_ptr kernel::localFlux::ExecutePtrs[];
  void kernel::localFlux::execute0() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(0) != nullptr);
    assert(rDivM(0) != nullptr);
    float *_tmp0, *_tmp1;
    // TODO: allocate all buffers at the main entry point of the program
    // allocating temp memory only on gpu
    float *d_buffer0 = (float*)DeviceScratchMem::get_instance().get_mem(216 * tensor::num_elements_in_cluster );
    // allocating temp memory only on gpu
    float *d_buffer1 = (float*)DeviceScratchMem::get_instance().get_mem(216 * tensor::num_elements_in_cluster );
    _tmp0 = d_buffer0;
    {
      unsigned *jump_to_next_fMrT = tensor::fMrT::stride_ptr[0];
      unsigned *jump_to_next_I = tensor::I::stride_ptr;
      const unsigned jump_to_next__tmp0 = 216;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(0), 24, I, 56, 0.0, _tmp0, 24, jump_to_next_fMrT, jump_to_next_I, jump_to_next__tmp0, tensor::num_elements_in_cluster);
    }
    _tmp1 = d_buffer1;
    {
      const unsigned jump_to_next__tmp0 = 216;
      unsigned *jump_to_next_AplusT = tensor::AplusT::stride_ptr;
      const unsigned jump_to_next__tmp1 = 216;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, jump_to_next__tmp0, jump_to_next_AplusT, jump_to_next__tmp1, tensor::num_elements_in_cluster);
    }
    {
      unsigned *jump_to_next_rDivM = tensor::rDivM::stride_ptr[0];
      const unsigned jump_to_next__tmp1 = 216;
      unsigned *jump_to_next_Q = tensor::Q::stride_ptr;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp1, 24, 1.0, Q, 56, jump_to_next_rDivM, jump_to_next__tmp1, jump_to_next_Q, tensor::num_elements_in_cluster);
    }
    DeviceScratchMem::get_instance().free();
    DeviceScratchMem::get_instance().free();
  }
  void kernel::localFlux::execute1() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(1) != nullptr);
    assert(rDivM(1) != nullptr);
    float *_tmp0, *_tmp1;
    // TODO: allocate all buffers at the main entry point of the program
    // allocating temp memory only on gpu
    float *d_buffer0 = (float*)DeviceScratchMem::get_instance().get_mem(216 * tensor::num_elements_in_cluster );
    // allocating temp memory only on gpu
    float *d_buffer1 = (float*)DeviceScratchMem::get_instance().get_mem(216 * tensor::num_elements_in_cluster );
    _tmp0 = d_buffer0;
    {
      unsigned* jump_to_next_fMrT = tensor::fMrT::stride_ptr[1];
      unsigned* jump_to_next_I = tensor::I::stride_ptr;
      const unsigned jump_to_next__tmp0 = 216;
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(1), 24, I, 56, 0.0, _tmp0, 24, jump_to_next_fMrT, jump_to_next_I, jump_to_next__tmp0, tensor::num_elements_in_cluster);
    }
    _tmp1 = d_buffer1;
    {
      const unsigned jump_to_next__tmp0 = 216;
      unsigned* jump_to_next_AplusT = tensor::AplusT::stride_ptr;
      const unsigned jump_to_next__tmp1 = 216;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, jump_to_next__tmp0, jump_to_next_AplusT, jump_to_next__tmp1, tensor::num_elements_in_cluster);
    }
    {
      unsigned* jump_to_next_rDivM = tensor::rDivM::stride_ptr[1];
      const unsigned jump_to_next__tmp1 = 216;
      unsigned* jump_to_next_Q = tensor::Q::stride_ptr;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp1, 24, 1.0, Q, 56, jump_to_next_rDivM, jump_to_next__tmp1, jump_to_next_Q, tensor::num_elements_in_cluster);
    }
    DeviceScratchMem::get_instance().free();
    DeviceScratchMem::get_instance().free();
  }
  void kernel::localFlux::execute2() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(2) != nullptr);
    assert(rDivM(2) != nullptr);
    float *_tmp0, *_tmp1;
    // TODO: allocate all buffers at the main entry point of the program
    // allocating temp memory only on gpu
    float *d_buffer0 = (float*)DeviceScratchMem::get_instance().get_mem(216 * tensor::num_elements_in_cluster );
    // allocating temp memory only on gpu
    float *d_buffer1 = (float*)DeviceScratchMem::get_instance().get_mem(216 * tensor::num_elements_in_cluster );
    _tmp0 = d_buffer0;
    {
      unsigned* jump_to_next_fMrT = tensor::fMrT::stride_ptr[2];
      unsigned* jump_to_next_I = tensor::I::stride_ptr;
      const unsigned jump_to_next__tmp0 = 216;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(2), 24, I, 56, 0.0, _tmp0, 24, jump_to_next_fMrT, jump_to_next_I, jump_to_next__tmp0, tensor::num_elements_in_cluster);
    }
    _tmp1 = d_buffer1;
    {
      const unsigned jump_to_next__tmp0 = 216;
      unsigned* jump_to_next_AplusT = tensor::AplusT::stride_ptr;
      const unsigned jump_to_next__tmp1 = 216;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, jump_to_next__tmp0, jump_to_next_AplusT, jump_to_next__tmp1, tensor::num_elements_in_cluster);
    }
    {
      unsigned* jump_to_next_rDivM = tensor::rDivM::stride_ptr[2];
      const unsigned jump_to_next__tmp1 = 216;
      //unsigned* jump_to_next_Q = tensor::Q::stride_ptr;
      unsigned* jump_to_next_Q = tensor::Q::stride_ptr;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp1, 24, 1.0, Q, 56, jump_to_next_rDivM, jump_to_next__tmp1, jump_to_next_Q, tensor::num_elements_in_cluster);
    }
    DeviceScratchMem::get_instance().free();
    DeviceScratchMem::get_instance().free();
  }
  void kernel::localFlux::execute3() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(3) != nullptr);
    assert(rDivM(3) != nullptr);
    float *_tmp0, *_tmp1;
    // TODO: allocate all buffers at the main entry point of the program
    // allocating temp memory only on gpu
    float *d_buffer0 = (float*)DeviceScratchMem::get_instance().get_mem(216 * tensor::num_elements_in_cluster );
    // allocating temp memory only on gpu
    float *d_buffer1 = (float*)DeviceScratchMem::get_instance().get_mem(216 * tensor::num_elements_in_cluster );
    _tmp0 = d_buffer0;
    {
      unsigned* jump_to_next_fMrT = tensor::fMrT::stride_ptr[3];
      unsigned* jump_to_next_I = tensor::I::stride_ptr;
      const unsigned jump_to_next__tmp0 = 216;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(3), 24, I, 56, 0.0, _tmp0, 24, jump_to_next_fMrT, jump_to_next_I, jump_to_next__tmp0, tensor::num_elements_in_cluster);
    }
    _tmp1 = d_buffer1;
    {
      const unsigned jump_to_next__tmp0 = 216;
      unsigned* jump_to_next_AplusT = tensor::AplusT::stride_ptr;
      const unsigned jump_to_next__tmp1 = 216;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24, jump_to_next__tmp0, jump_to_next_AplusT, jump_to_next__tmp1, tensor::num_elements_in_cluster);
    }
    {
      unsigned* jump_to_next_rDivM = tensor::rDivM::stride_ptr[3];
      const unsigned jump_to_next__tmp1 = 216;
      unsigned* jump_to_next_Q = tensor::Q::stride_ptr;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp1, 24, 1.0, Q, 56, jump_to_next_rDivM, jump_to_next__tmp1, jump_to_next_Q, tensor::num_elements_in_cluster);
    }
    DeviceScratchMem::get_instance().free();
    DeviceScratchMem::get_instance().free();
  }
  #else
  constexpr unsigned long const kernel::localFlux::NonZeroFlops[];
  constexpr unsigned long const kernel::localFlux::HardwareFlops[];
  constexpr kernel::localFlux::member_function_ptr kernel::localFlux::ExecutePtrs[];
  void kernel::localFlux::execute0() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(0) != nullptr);
    assert(rDivM(0) != nullptr);
    float *_tmp0, *_tmp1;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp1, 24, 1.0, Q, 56);
  }
  void kernel::localFlux::execute1() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(1) != nullptr);
    assert(rDivM(1) != nullptr);
    float *_tmp0, *_tmp1;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp1, 24, 1.0, Q, 56);
  }
  void kernel::localFlux::execute2() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(2) != nullptr);
    assert(rDivM(2) != nullptr);
    float *_tmp0, *_tmp1;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp1, 24, 1.0, Q, 56);
  }
  void kernel::localFlux::execute3() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(3) != nullptr);
    assert(rDivM(3) != nullptr);
    float *_tmp0, *_tmp1;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, fMrT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, AplusT, 9, 0.0, _tmp1, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp1, 24, 1.0, Q, 56);
  }
  #endif
  constexpr unsigned long const kernel::neighboringFlux::NonZeroFlops[];
  constexpr unsigned long const kernel::neighboringFlux::HardwareFlops[];
  constexpr kernel::neighboringFlux::member_function_ptr kernel::neighboringFlux::ExecutePtrs[];
  void kernel::neighboringFlux::execute0() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute1() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute2() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute3() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute4() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute5() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute6() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute7() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute8() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute9() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute10() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute11() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(0), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute12() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute13() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute14() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute15() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute16() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute17() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute18() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute19() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute20() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute21() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute22() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute23() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(1), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute24() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute25() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute26() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute27() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute28() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute29() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute30() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute31() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute32() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute33() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute34() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute35() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(2), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute36() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute37() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute38() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(0), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute39() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute40() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute41() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(1), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute42() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute43() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute44() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(2), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute45() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(0), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute46() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(1), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  void kernel::neighboringFlux::execute47() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    float *_tmp0, *_tmp1, *_tmp2;
    float _buffer0[216] __attribute__((aligned(32)));
    float _buffer1[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 56, 1.0, rT(3), 24, I, 56, 0.0, _tmp0, 24);
    _tmp1 = _buffer1;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 21, 1.0, fP(2), 24, _tmp0, 24, 0.0, _tmp1, 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp1, 24, AminusT, 9, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 21, 1.0, rDivM(3), 56, _tmp2, 24, 1.0, Q, 56);
  }
  #ifdef GPU
  constexpr unsigned long const kernel::derivativeTaylorExpansion::NonZeroFlops[];
  constexpr unsigned long const kernel::derivativeTaylorExpansion::HardwareFlops[];
  constexpr kernel::derivativeTaylorExpansion::member_function_ptr kernel::derivativeTaylorExpansion::ExecutePtrs[];
  void kernel::derivativeTaylorExpansion::execute0() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(0) != nullptr);
    // TODO: allocate all buffers at the main entry point of the program
    cuda_copy_add_scale(56, 9, power, dQ(0) + 0, 56, 0.0, I + 0, 56, tensor::dQ::jump_to_next[0], tensor::I::jump_to_next, tensor::num_elements_in_cluster);
  }
  void kernel::derivativeTaylorExpansion::execute1() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(1) != nullptr);
    // TODO: allocate all buffers at the main entry point of the program
    cuda_copy_add_scale(35, 9, power, dQ(1) + 0, 40, 1.0, I + 0, 56, tensor::dQ::jump_to_next[1], tensor::I::jump_to_next, tensor::num_elements_in_cluster);
  }
  void kernel::derivativeTaylorExpansion::execute2() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(2) != nullptr);
    // TODO: allocate all buffers at the main entry point of the program
    cuda_copy_add_scale(20, 9, power, dQ(2) + 0, 24, 1.0, I + 0, 56, tensor::dQ::jump_to_next[2], tensor::I::jump_to_next, tensor::num_elements_in_cluster);
  }
  void kernel::derivativeTaylorExpansion::execute3() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(3) != nullptr);
    // TODO: allocate all buffers at the main entry point of the program
    cuda_copy_add_scale(10, 9, power, dQ(3) + 0, 16, 1.0, I + 0, 56, tensor::dQ::jump_to_next[3], tensor::I::jump_to_next, tensor::num_elements_in_cluster);
  }
  void kernel::derivativeTaylorExpansion::execute4() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(4) != nullptr);
    // TODO: allocate all buffers at the main entry point of the program
    cuda_copy_add_scale(4, 9, power, dQ(4) + 0, 8, 1.0, I + 0, 56, tensor::dQ::jump_to_next[4], tensor::I::jump_to_next, tensor::num_elements_in_cluster);
  }
  void kernel::derivativeTaylorExpansion::execute5() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(5) != nullptr);
    // TODO: allocate all buffers at the main entry point of the program
    cuda_copy_add_scale(1, 9, power, dQ(5) + 0, 8, 1.0, I + 0, 56, tensor::dQ::jump_to_next[5], tensor::I::jump_to_next, tensor::num_elements_in_cluster);
  }
  #else
  constexpr unsigned long const kernel::derivativeTaylorExpansion::NonZeroFlops[];
  constexpr unsigned long const kernel::derivativeTaylorExpansion::HardwareFlops[];
  constexpr kernel::derivativeTaylorExpansion::member_function_ptr kernel::derivativeTaylorExpansion::ExecutePtrs[];
  void kernel::derivativeTaylorExpansion::execute0() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(0) != nullptr);
    for (int b = 0; b < 9; ++b) {
      #pragma omp simd
      for (int a = 0; a < 56; ++a) {
        I[1*a + 56*b] = power * dQ(0)[1*a + 56*b];
      }
    }
  }
  void kernel::derivativeTaylorExpansion::execute1() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(1) != nullptr);
    for (int b = 0; b < 9; ++b) {
      #pragma omp simd
      for (int a = 0; a < 35; ++a) {
        I[1*a + 56*b] += power * dQ(1)[1*a + 40*b];
      }
    }
  }
  void kernel::derivativeTaylorExpansion::execute2() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(2) != nullptr);
    for (int b = 0; b < 9; ++b) {
      #pragma omp simd
      for (int a = 0; a < 20; ++a) {
        I[1*a + 56*b] += power * dQ(2)[1*a + 24*b];
      }
    }
  }
  void kernel::derivativeTaylorExpansion::execute3() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(3) != nullptr);
    for (int b = 0; b < 9; ++b) {
      #pragma omp simd
      for (int a = 0; a < 10; ++a) {
        I[1*a + 56*b] += power * dQ(3)[1*a + 16*b];
      }
    }
  }
  void kernel::derivativeTaylorExpansion::execute4() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(4) != nullptr);
    for (int b = 0; b < 9; ++b) {
      #pragma omp simd
      for (int a = 0; a < 4; ++a) {
        I[1*a + 56*b] += power * dQ(4)[1*a + 8*b];
      }
    }
  }
  void kernel::derivativeTaylorExpansion::execute5() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(5) != nullptr);
    for (int b = 0; b < 9; ++b) {
      #pragma omp simd
      for (int a = 0; a < 1; ++a) {
        I[1*a + 56*b] += power * dQ(5)[1*a + 8*b];
      }
    }
  }
  #endif

  #ifdef GPU
  constexpr unsigned long const kernel::derivative::NonZeroFlops[];
  constexpr unsigned long const kernel::derivative::HardwareFlops[];
  constexpr kernel::derivative::member_function_ptr kernel::derivative::ExecutePtrs[];
  void kernel::derivative::execute1() {
    assert(dQ(0) != nullptr);
    assert(dQ(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    // TODO: allocate all buffers at the main entry point of the program
    // allocating temp memory only on gpu
    float *d_buffer0 = (float*)DeviceScratchMem::get_instance().get_mem(504 * tensor::num_elements_in_cluster );
    _tmp0 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[0];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[0];
    const unsigned jump_to_next__tmp0 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 53, 1.0, kDivMT(0), 40, dQ(0) + 1, 56, 0.0, _tmp0, 40, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp0, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp0 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[0];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[1];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, _tmp0, 40, star(0), 9, 0.0, dQ(1), 40, jump_to_next__tmp0, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    _tmp2 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[1];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[0];
    const unsigned jump_to_next__tmp2 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 54, 1.0, kDivMT(1), 40, dQ(0) + 1, 56, 0.0, _tmp2, 40, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp2, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp2 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[1];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[1];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, _tmp2, 40, star(1), 9, 1.0, dQ(1), 40, jump_to_next__tmp2, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    _tmp4 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[2];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[0];
    const unsigned jump_to_next__tmp4 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 55, 1.0, kDivMT(2), 40, dQ(0) + 1, 56, 0.0, _tmp4, 40, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp4, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp4 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[2];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[1];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, _tmp4, 40, star(2), 9, 1.0, dQ(1), 40, jump_to_next__tmp4, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    DeviceScratchMem::get_instance().free();
  }
  void kernel::derivative::execute2() {
    assert(dQ(2) != nullptr);
    assert(dQ(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    // TODO: allocate all buffers at the main entry point of the program
    // allocating temp memory only on gpu
    float *d_buffer0 = (float*)DeviceScratchMem::get_instance().get_mem(504 * tensor::num_elements_in_cluster );
    _tmp0 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[0];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[1];
    const unsigned jump_to_next__tmp0 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 32, 1.0, kDivMT(0), 40, dQ(1) + 1, 40, 0.0, _tmp0, 24, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp0, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp0 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[0];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[2];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, star(0), 9, 0.0, dQ(2), 24, jump_to_next__tmp0, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    _tmp2 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[1];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[1];
    const unsigned jump_to_next__tmp2 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 33, 1.0, kDivMT(1), 40, dQ(1) + 1, 40, 0.0, _tmp2, 24, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp2, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp2 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[1];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[2];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp2, 24, star(1), 9, 1.0, dQ(2), 24, jump_to_next__tmp2, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    _tmp4 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[2];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[1];
    const unsigned jump_to_next__tmp4 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 34, 1.0, kDivMT(2), 40, dQ(1) + 1, 40, 0.0, _tmp4, 24, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp4, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp4 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[2];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[2];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp4, 24, star(2), 9, 1.0, dQ(2), 24, jump_to_next__tmp4, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    DeviceScratchMem::get_instance().free();
  }
  void kernel::derivative::execute3() {
    assert(dQ(2) != nullptr);
    assert(dQ(3) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    // TODO: allocate all buffers at the main entry point of the program
    // allocating temp memory only on gpu
    float *d_buffer0 = (float*)DeviceScratchMem::get_instance().get_mem(504 * tensor::num_elements_in_cluster );
    _tmp0 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[0];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[2];
    const unsigned jump_to_next__tmp0 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 17, 1.0, kDivMT(0), 40, dQ(2) + 1, 24, 0.0, _tmp0, 16, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp0, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp0 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[0];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[3];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 9, 1.0, _tmp0, 16, star(0), 9, 0.0, dQ(3), 16, jump_to_next__tmp0, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    _tmp2 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[1];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[2];
    const unsigned jump_to_next__tmp2 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 18, 1.0, kDivMT(1), 40, dQ(2) + 1, 24, 0.0, _tmp2, 16, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp2, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp2 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[1];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[3];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 9, 1.0, _tmp2, 16, star(1), 9, 1.0, dQ(3), 16, jump_to_next__tmp2, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    _tmp4 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[2];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[2];
    const unsigned jump_to_next__tmp4 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 19, 1.0, kDivMT(2), 40, dQ(2) + 1, 24, 0.0, _tmp4, 16, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp4, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp4 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[2];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[3];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 9, 1.0, _tmp4, 16, star(2), 9, 1.0, dQ(3), 16, jump_to_next__tmp4, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    DeviceScratchMem::get_instance().free();
  }
  void kernel::derivative::execute4() {
    assert(dQ(3) != nullptr);
    assert(dQ(4) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    // TODO: allocate all buffers at the main entry point of the program
    // allocating temp memory only on gpu
    float *d_buffer0 = (float*)DeviceScratchMem::get_instance().get_mem(504 * tensor::num_elements_in_cluster );
    _tmp0 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[0];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[3];
    const unsigned jump_to_next__tmp0 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 7, 1.0, kDivMT(0), 40, dQ(3) + 1, 16, 0.0, _tmp0, 8, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp0, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp0 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[0];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[4];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp0, 8, star(0), 9, 0.0, dQ(4), 8, jump_to_next__tmp0, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    _tmp2 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[1];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[3];
    const unsigned jump_to_next__tmp2 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 8, 1.0, kDivMT(1), 40, dQ(3) + 1, 16, 0.0, _tmp2, 8, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp2, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp2 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[1];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[4];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp2, 8, star(1), 9, 1.0, dQ(4), 8, jump_to_next__tmp2, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    _tmp4 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[2];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[3];
    const unsigned jump_to_next__tmp4 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, kDivMT(2), 40, dQ(3) + 1, 16, 0.0, _tmp4, 8, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp4, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp4 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[2];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[4];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp4, 8, star(2), 9, 1.0, dQ(4), 8, jump_to_next__tmp4, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    DeviceScratchMem::get_instance().free();
  }
  void kernel::derivative::execute5() {
    assert(dQ(5) != nullptr);
    assert(dQ(4) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    // TODO: allocate all buffers at the main entry point of the program
    // allocating temp memory only on gpu
    float *d_buffer0 = (float*)DeviceScratchMem::get_instance().get_mem(504 * tensor::num_elements_in_cluster );
    _tmp0 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[0];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[4];
    const unsigned jump_to_next__tmp0 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 1, 1.0, kDivMT(0), 40, dQ(4) + 1, 8, 0.0, _tmp0, 8, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp0, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp0 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[0];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[5];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp0, 8, star(0), 9, 0.0, dQ(5), 8, jump_to_next__tmp0, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    _tmp2 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[1];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[4];
    const unsigned jump_to_next__tmp2 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 2, 1.0, kDivMT(1), 40, dQ(4) + 1, 8, 0.0, _tmp2, 8, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp2, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp2 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[1];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[5];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp2, 8, star(1), 9, 1.0, dQ(5), 8, jump_to_next__tmp2, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    _tmp4 = d_buffer0;
    {
    const unsigned jump_to_next_kDivMT = tensor::kDivMT::jump_to_next[2];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[4];
    const unsigned jump_to_next__tmp4 = 504;
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 3, 1.0, kDivMT(2), 40, dQ(4) + 1, 8, 0.0, _tmp4, 8, jump_to_next_kDivMT, jump_to_next_dQ, jump_to_next__tmp4, tensor::num_elements_in_cluster);
    }
    {
    const unsigned jump_to_next__tmp4 = 504;
    const unsigned jump_to_next_star = tensor::star::jump_to_next[2];
    const unsigned jump_to_next_dQ = tensor::dQ::jump_to_next[5];
    
    
    cuda_blas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp4, 8, star(2), 9, 1.0, dQ(5), 8, jump_to_next__tmp4, jump_to_next_star, jump_to_next_dQ, tensor::num_elements_in_cluster);
    }
    DeviceScratchMem::get_instance().free();
  }
  #else
  constexpr unsigned long const kernel::derivative::NonZeroFlops[];
  constexpr unsigned long const kernel::derivative::HardwareFlops[];
  constexpr kernel::derivative::member_function_ptr kernel::derivative::ExecutePtrs[];
  void kernel::derivative::execute1() {
    assert(dQ(0) != nullptr);
    assert(dQ(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    float _buffer0[360] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 53, 1.0, kDivMT(0), 40, dQ(0) + 1, 56, 0.0, _tmp0, 40);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, _tmp0, 40, star(0), 9, 0.0, dQ(1), 40);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 54, 1.0, kDivMT(1), 40, dQ(0) + 1, 56, 0.0, _tmp2, 40);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, _tmp2, 40, star(1), 9, 1.0, dQ(1), 40);
    _tmp4 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 55, 1.0, kDivMT(2), 40, dQ(0) + 1, 56, 0.0, _tmp4, 40);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 40, 9, 9, 1.0, _tmp4, 40, star(2), 9, 1.0, dQ(1), 40);
  }
  void kernel::derivative::execute2() {
    assert(dQ(2) != nullptr);
    assert(dQ(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    float _buffer0[216] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 32, 1.0, kDivMT(0), 40, dQ(1) + 1, 40, 0.0, _tmp0, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp0, 24, star(0), 9, 0.0, dQ(2), 24);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 33, 1.0, kDivMT(1), 40, dQ(1) + 1, 40, 0.0, _tmp2, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp2, 24, star(1), 9, 1.0, dQ(2), 24);
    _tmp4 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 34, 1.0, kDivMT(2), 40, dQ(1) + 1, 40, 0.0, _tmp4, 24);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 24, 9, 9, 1.0, _tmp4, 24, star(2), 9, 1.0, dQ(2), 24);
  }
  void kernel::derivative::execute3() {
    assert(dQ(2) != nullptr);
    assert(dQ(3) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    float _buffer0[144] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 17, 1.0, kDivMT(0), 40, dQ(2) + 1, 24, 0.0, _tmp0, 16);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 9, 1.0, _tmp0, 16, star(0), 9, 0.0, dQ(3), 16);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 18, 1.0, kDivMT(1), 40, dQ(2) + 1, 24, 0.0, _tmp2, 16);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 9, 1.0, _tmp2, 16, star(1), 9, 1.0, dQ(3), 16);
    _tmp4 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 19, 1.0, kDivMT(2), 40, dQ(2) + 1, 24, 0.0, _tmp4, 16);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 9, 9, 1.0, _tmp4, 16, star(2), 9, 1.0, dQ(3), 16);
  }
  void kernel::derivative::execute4() {
    assert(dQ(3) != nullptr);
    assert(dQ(4) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    float _buffer0[72] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 7, 1.0, kDivMT(0), 40, dQ(3) + 1, 16, 0.0, _tmp0, 8);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp0, 8, star(0), 9, 0.0, dQ(4), 8);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 8, 1.0, kDivMT(1), 40, dQ(3) + 1, 16, 0.0, _tmp2, 8);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp2, 8, star(1), 9, 1.0, dQ(4), 8);
    _tmp4 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, kDivMT(2), 40, dQ(3) + 1, 16, 0.0, _tmp4, 8);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp4, 8, star(2), 9, 1.0, dQ(4), 8);
  }
  void kernel::derivative::execute5() {
    assert(dQ(5) != nullptr);
    assert(dQ(4) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(star(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    float *_tmp0, *_tmp2, *_tmp4;
    float _buffer0[72] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 1, 1.0, kDivMT(0), 40, dQ(4) + 1, 8, 0.0, _tmp0, 8);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp0, 8, star(0), 9, 0.0, dQ(5), 8);
    _tmp2 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 2, 1.0, kDivMT(1), 40, dQ(4) + 1, 8, 0.0, _tmp2, 8);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp2, 8, star(1), 9, 1.0, dQ(5), 8);
    _tmp4 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 3, 1.0, kDivMT(2), 40, dQ(4) + 1, 8, 0.0, _tmp4, 8);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 9, 9, 1.0, _tmp4, 8, star(2), 9, 1.0, dQ(5), 8);
  }
  #endif

  constexpr unsigned long const kernel::godunovState::NonZeroFlops[];
  constexpr unsigned long const kernel::godunovState::HardwareFlops[];
  constexpr kernel::godunovState::member_function_ptr kernel::godunovState::ExecutePtrs[];
  void kernel::godunovState::execute0() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(0,0), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 0.0, godunovState, 56);
  }
  void kernel::godunovState::execute1() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(1,0), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 0.0, godunovState, 56);
  }
  void kernel::godunovState::execute2() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(2,0), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 0.0, godunovState, 56);
  }
  void kernel::godunovState::execute3() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,0) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(3,0), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 0.0, godunovState, 56);
  }
  void kernel::godunovState::execute4() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(0,1), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute5() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(1,1), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute6() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(2,1), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute7() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,1) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(3,1), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute8() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(0,2), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute9() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(1,2), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute10() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(2,2), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute11() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,2) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(3,2), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute12() {
    assert(Q != nullptr);
    assert(V3mTo2n(0,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(0,3), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute13() {
    assert(Q != nullptr);
    assert(V3mTo2n(1,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(1,3), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute14() {
    assert(Q != nullptr);
    assert(V3mTo2n(2,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(2,3), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  void kernel::godunovState::execute15() {
    assert(Q != nullptr);
    assert(V3mTo2n(3,3) != nullptr);
    assert(godunovMatrix != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 56, 1.0, V3mTo2n(3,3), 56, Q, 56, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, _tmp0, 56, godunovMatrix, 9, 1.0, godunovState, 56);
  }
  constexpr unsigned long const kernel::nodalFlux::NonZeroFlops[];
  constexpr unsigned long const kernel::nodalFlux::HardwareFlops[];
  constexpr kernel::nodalFlux::member_function_ptr kernel::nodalFlux::ExecutePtrs[];
  void kernel::nodalFlux::execute0() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,0), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute1() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,0), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute2() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,0), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute3() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,0) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,0), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute4() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,1), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute5() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,1), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute6() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,1), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute7() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,1) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,1), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute8() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,2), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute9() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,2), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute10() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,2), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute11() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,2) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,2), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute12() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(0,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(0,3), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute13() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(1,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(1,3), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute14() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(2,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(2,3), 56, _tmp0, 56, 1.0, Q, 56);
  }
  void kernel::nodalFlux::execute15() {
    assert(Q != nullptr);
    assert(V3mTo2nTWDivM(3,3) != nullptr);
    assert(fluxSolver != nullptr);
    assert(godunovState != nullptr);
    float *_tmp0;
    float _buffer0[504] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 9, 1.0, godunovState, 56, fluxSolver, 9, 0.0, _tmp0, 56);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 56, 9, 49, 1.0, V3mTo2nTWDivM(3,3), 56, _tmp0, 56, 1.0, Q, 56);
  }
  constexpr unsigned long const kernel::subTriangleDisplacement::NonZeroFlops[];
  constexpr unsigned long const kernel::subTriangleDisplacement::HardwareFlops[];
  constexpr kernel::subTriangleDisplacement::member_function_ptr kernel::subTriangleDisplacement::ExecutePtrs[];
  void kernel::subTriangleDisplacement::execute0() {
    assert(displacement != nullptr);
    assert(subTriangleDofs(0) != nullptr);
    assert(subTriangleProjection(0) != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 3, 56, 1.0, subTriangleProjection(0), 8, displacement, 56, 0.0, subTriangleDofs(0), 8);
  }
  void kernel::subTriangleDisplacement::execute1() {
    assert(displacement != nullptr);
    assert(subTriangleDofs(1) != nullptr);
    assert(subTriangleProjection(1) != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 3, 56, 1.0, subTriangleProjection(1), 8, displacement, 56, 0.0, subTriangleDofs(1), 8);
  }
  void kernel::subTriangleDisplacement::execute2() {
    assert(displacement != nullptr);
    assert(subTriangleDofs(2) != nullptr);
    assert(subTriangleProjection(2) != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 3, 56, 1.0, subTriangleProjection(2), 16, displacement, 56, 0.0, subTriangleDofs(2), 16);
  }
  void kernel::subTriangleDisplacement::execute3() {
    assert(displacement != nullptr);
    assert(subTriangleDofs(3) != nullptr);
    assert(subTriangleProjection(3) != nullptr);
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 64, 3, 56, 1.0, subTriangleProjection(3), 64, displacement, 56, 0.0, subTriangleDofs(3), 64);
  }
  constexpr unsigned long const kernel::subTriangleVelocity::NonZeroFlops[];
  constexpr unsigned long const kernel::subTriangleVelocity::HardwareFlops[];
  constexpr kernel::subTriangleVelocity::member_function_ptr kernel::subTriangleVelocity::ExecutePtrs[];
  void kernel::subTriangleVelocity::execute0() {
    assert(Q != nullptr);
    assert(selectVelocity != nullptr);
    assert(subTriangleDofs(0) != nullptr);
    assert(subTriangleProjection(0) != nullptr);
    float *_tmp0;
    float _buffer0[24] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 3, 56, 1.0, subTriangleProjection(0), 8, Q + 336, 56, 0.0, _tmp0, 8);
    for (int o = 0; o < 8; ++o) {
      for (int i = 0; i < 3; ++i) {
        subTriangleDofs(0)[0 + 1*o + 8*i] = 0.0;
      }
      subTriangleDofs(0)[0 + 1*o + 8*0] += 1.0 * _tmp0[0 + 1*o + 8*0] * selectVelocity[0 + 0];
      subTriangleDofs(0)[0 + 1*o + 8*1] += 1.0 * _tmp0[0 + 1*o + 8*1] * selectVelocity[0 + 1];
      subTriangleDofs(0)[0 + 1*o + 8*2] += 1.0 * _tmp0[0 + 1*o + 8*2] * selectVelocity[0 + 2];
    }
  }
  void kernel::subTriangleVelocity::execute1() {
    assert(Q != nullptr);
    assert(selectVelocity != nullptr);
    assert(subTriangleDofs(1) != nullptr);
    assert(subTriangleProjection(1) != nullptr);
    float *_tmp0;
    float _buffer0[24] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 8, 3, 56, 1.0, subTriangleProjection(1), 8, Q + 336, 56, 0.0, _tmp0, 8);
    for (int o = 0; o < 8; ++o) {
      for (int i = 0; i < 3; ++i) {
        subTriangleDofs(1)[0 + 1*o + 8*i] = 0.0;
      }
      subTriangleDofs(1)[0 + 1*o + 8*0] += 1.0 * _tmp0[0 + 1*o + 8*0] * selectVelocity[0 + 0];
      subTriangleDofs(1)[0 + 1*o + 8*1] += 1.0 * _tmp0[0 + 1*o + 8*1] * selectVelocity[0 + 1];
      subTriangleDofs(1)[0 + 1*o + 8*2] += 1.0 * _tmp0[0 + 1*o + 8*2] * selectVelocity[0 + 2];
    }
  }
  void kernel::subTriangleVelocity::execute2() {
    assert(Q != nullptr);
    assert(selectVelocity != nullptr);
    assert(subTriangleDofs(2) != nullptr);
    assert(subTriangleProjection(2) != nullptr);
    float *_tmp0;
    float _buffer0[48] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 16, 3, 56, 1.0, subTriangleProjection(2), 16, Q + 336, 56, 0.0, _tmp0, 16);
    for (int o = 0; o < 16; ++o) {
      for (int i = 0; i < 3; ++i) {
        subTriangleDofs(2)[0 + 1*o + 16*i] = 0.0;
      }
      subTriangleDofs(2)[0 + 1*o + 16*0] += 1.0 * _tmp0[0 + 1*o + 16*0] * selectVelocity[0 + 0];
      subTriangleDofs(2)[0 + 1*o + 16*1] += 1.0 * _tmp0[0 + 1*o + 16*1] * selectVelocity[0 + 1];
      subTriangleDofs(2)[0 + 1*o + 16*2] += 1.0 * _tmp0[0 + 1*o + 16*2] * selectVelocity[0 + 2];
    }
  }
  void kernel::subTriangleVelocity::execute3() {
    assert(Q != nullptr);
    assert(selectVelocity != nullptr);
    assert(subTriangleDofs(3) != nullptr);
    assert(subTriangleProjection(3) != nullptr);
    float *_tmp0;
    float _buffer0[168] __attribute__((aligned(32)));
    _tmp0 = _buffer0;
    for (int o = 0; o < 56; ++o) {
      for (int i = 0; i < 3; ++i) {
        _tmp0[0 + 1*o + 56*i] = 0.0;
      }
      _tmp0[0 + 1*o + 56*0] += 1.0 * Q[336 + 1*o + 56*0] * selectVelocity[0 + 0];
      _tmp0[0 + 1*o + 56*1] += 1.0 * Q[336 + 1*o + 56*1] * selectVelocity[0 + 1];
      _tmp0[0 + 1*o + 56*2] += 1.0 * Q[336 + 1*o + 56*2] * selectVelocity[0 + 2];
    }
    
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 64, 3, 56, 1.0, subTriangleProjection(3), 64, _tmp0, 56, 0.0, subTriangleDofs(3), 64);
  }
}
