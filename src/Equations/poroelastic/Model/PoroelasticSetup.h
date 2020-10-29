#include <armadillo>

#include <yateto/TensorView.h>

#include "Model/common.hpp"
#include "Kernels/common.hpp"
#include "Numerical_aux/Transformation.h"
#include "generated_code/init.h"
#include "PoroelasticJacobian.h"

namespace seissol {
  namespace model {
    template<typename T>
    inline void getTransposedSourceCoefficientTensor( PoroElasticMaterial const& material, T& ET) {
      ET.setZero();
      ET(10,10) = material.viscosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/(material.permeability*material.rho_fluid*(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity));
      ET(11,11) = material.viscosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/(material.permeability*material.rho_fluid*(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity));
      ET(12,12) = material.viscosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/(material.permeability*material.rho_fluid*(material.rho_fluid - material.tortuosity*(material.porosity*material.rho_fluid - material.rho*(material.porosity - 1))/material.porosity));

      ET(10,6) = material.porosity*material.viscosity/(material.permeability*material.tortuosity*(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1)));
      ET(11,7) = material.porosity*material.viscosity/(material.permeability*material.tortuosity*(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1)));
      ET(12,8) = material.porosity*material.viscosity/(material.permeability*material.tortuosity*(material.porosity*material.rho_fluid - material.porosity*material.rho_fluid/material.tortuosity - material.rho*(material.porosity - 1)));
    }

    template<>
    inline void getTransposedGodunovState( PoroElasticMaterial const&        local,
        PoroElasticMaterial const&        neighbor,
        FaceType                          faceType,
        init::QgodLocal::view::type&      QgodLocal,
        init::QgodNeighbor::view::type&   QgodNeighbor )
    {
      constexpr auto tolerance = 1.0e-4;

      using Matrix = typename arma::Mat<std::complex<double>>::template fixed<NUMBER_OF_QUANTITIES, NUMBER_OF_QUANTITIES>;
      using Vector = typename arma::Col<std::complex<double>>::template fixed<NUMBER_OF_QUANTITIES>;
      struct eigenDecomposition { Vector eigenvalues; Matrix eigenvectors; };
      auto getEigenDecomposition = [&tolerance](PoroElasticMaterial const& material) {
        Matrix coeff(arma::fill::zeros);
        getTransposedCoefficientMatrix(material, 0, coeff);

	Matrix arma_eigenvectors(arma::fill::zeros);
	Vector arma_eigenvalues(arma::fill::zeros);
	arma::eig_gen(arma_eigenvalues, arma_eigenvectors, coeff.t());

#ifndef NDEBUG
        int ev_neg = 0;
        int ev_pos = 0;
        for (int i = 0; i < NUMBER_OF_QUANTITIES; ++i) {
          if (evs(i).real() < -tolerance) {
            ++ev_neg;
          } else if (evs(i).real() > tolerance) {
            ++ev_pos;
          }
        }
        assert(ev_neg == 4);
        assert(ev_pos == 4);
#endif
        return eigenDecomposition({arma_eigenvalues, arma_eigenvectors});
      };
      auto eigen_local = getEigenDecomposition(local);
      auto eigen_neighbor = getEigenDecomposition(neighbor); 	

      Matrix chi_minus(arma::fill::zeros);
      Matrix chi_plus(arma::fill::zeros);
      for(int i = 0; i < 13; i++) {
        if(eigen_local.eigenvalues(i).real() < -tolerance) {
          chi_minus(i,i) = 1.0;
        }
        if(eigen_local.eigenvalues(i).real() > tolerance) {
          chi_plus(i,i) = 1.0;
        }
      }

      //Matrix R = eigen_local.eigenvectors() * chi_minus + eigen_neighbor.eigenvectors() * chi_plus;
      Matrix R = eigen_local.eigenvectors;
      if (faceType == FaceType::freeSurface) {
        logWarning() << "Poroelastic Free Surface is not tested yet.";
        //Matrix R_real = R.real().eval();
        //getTransposedFreeSurfaceGodunovState(false, QgodLocal, QgodNeighbor, R_real);
      } else {
	Matrix R_inv = inv(R);
        Matrix godunov_minus = R * chi_minus * R_inv;
        Matrix godunov_plus =  R * chi_plus * R_inv;

        for (unsigned i = 0; i < QgodLocal.shape(1); ++i) {
          for (unsigned j = 0; j < QgodLocal.shape(0); ++j) {
            QgodLocal(j,i) = godunov_plus(i,j).real();
            QgodNeighbor(j,i) = godunov_minus(i,j).real();
#ifndef NDEBUG
	    assert(std::abs(godunov_plus(j,i).imag()) < tolerance)
	    assert(std::abs(godunov_minus(j,i).imag()) < tolerance)
#endif
          }
        }
      }
    }

    inline void initializeSpecificLocalData( PoroElasticMaterial const& material,
        PoroelasticLocalData* localData )
    {
    }
  }
}
