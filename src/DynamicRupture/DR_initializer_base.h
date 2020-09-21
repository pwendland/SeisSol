//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_INITIALIZER_BASE_H
#define SEISSOL_DR_INITIALIZER_BASE_H

#include <c++/8.3.0/iostream>
#include <c++/8.3.0/unordered_map>
#include <Solver/Interoperability.h>
#include <yaml-cpp/yaml.h>
#include "Initializer/InputAux.hpp"
#include <DynamicRupture/DR_Parameters.h>

namespace seissol {
    namespace initializers {
      struct BaseDrInitializer;
      struct Init_FL_2;
      struct Init_FL_3; //aging law
      struct Init_FL_33;
      struct Init_FL_103;
      struct Init_FL_103_Thermal;
    }
}

class seissol::initializers::BaseDrInitializer {
protected:
  static constexpr int numberOfPoints = tensor::QInterpolated::Shape[0];
  static constexpr int numOfPointsPadded = init::QInterpolated::Stop[0];
  //YAML::Node m_InputParam;
  dr::DrParameterT m_Params;

public:
  virtual ~BaseDrInitializer() {}

  //set the parameters from .par file with yaml to this class attributes.
  void setInputParam(const YAML::Node& Params) {
    using namespace initializers;
    //TODO: maybe allocate dr::DrParameterT in MemoryManager and copy here only the reference
    m_Params.setAllInputParam(Params);
  }


  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
        initializers::LTSTree* dynRupTree,
        std::unordered_map<std::string,
        double*> faultParameters,
        unsigned* ltsFaceToMeshFace,
        seissol::Interoperability &e_interoperability) {
    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
      real  (*initialStressInFaultCS)[numOfPointsPadded][6] = it->var(dynRup->initialStressInFaultCS);  //get from fortran  EQN%InitialStressInFaultCS
      real  (*cohesion)[numOfPointsPadded]                  = it->var(dynRup->cohesion);                //get from faultParameters
      real  (*mu)[ numOfPointsPadded ]            = it->var(dynRup->mu);                                //get from fortran  EQN%IniMu(:,:)
      real  (*slip)[ numOfPointsPadded ]          = it->var(dynRup->slip);                              // = 0
      real  (*slip1)[numOfPointsPadded ]          = it->var(dynRup->slip1);                             // = 0
      real  (*slip2)[ numOfPointsPadded ]         = it->var(dynRup->slip2);                             // = 0
      real  (*slipRate1)[ numOfPointsPadded ]     = it->var(dynRup->slipRate1);                         //get from fortran  EQN%IniSlipRate1
      real  (*slipRate2)[numOfPointsPadded ]      = it->var(dynRup->slipRate2);                         //get from fortran  EQN%IniSlipRate2
      real  (*rupture_time)[ numOfPointsPadded ]  = it->var(dynRup->rupture_time);                      // = 0
      bool  (*RF)[ numOfPointsPadded ]            = it->var(dynRup->RF);                                //get from fortran
      real  (*peakSR)[ numOfPointsPadded ]        = it->var(dynRup->peakSR);                            // = 0
      real  (*tracXY)[ numOfPointsPadded ]        = it->var(dynRup->tracXY);                            // = 0
      real  (*tracXZ)[ numOfPointsPadded ]        = it->var(dynRup->tracXZ);                            // = 0

      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
        for (unsigned iBndGP = 0; iBndGP < init::QInterpolated::Stop[0]; ++iBndGP) {    //loop includes padded elements
          slip[ltsFace][iBndGP] = 0.0;
          slip1[ltsFace][iBndGP] = 0.0;
          slip2[ltsFace][iBndGP] = 0.0;
          rupture_time[ltsFace][iBndGP] = 0.0;
          peakSR[ltsFace][iBndGP] = 0.0;
          tracXY[ltsFace][iBndGP] = 0.0;
          tracXZ[ltsFace][iBndGP] = 0.0;
        }
        //get initial values from fortran
        for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
          if(faultParameters["cohesion"] != NULL ){
            cohesion[ltsFace][iBndGP] = static_cast<real>( faultParameters["cohesion"][meshFace * numberOfPoints] );
          }else{
            //TODO: maybe not log it = too much spam?
            //std::cout << "DR_initializer_base: cohesion set to 0, not found from faultParameters";
            cohesion[ltsFace][iBndGP] = 0;
          }
        }
        //initialize padded elements for vectorization
        for (unsigned iBndGP = numberOfPoints; iBndGP < numOfPointsPadded; ++iBndGP) {
          cohesion[ltsFace][iBndGP]               = 0.0;
        }
        e_interoperability.getDynRupParameters(ltsFace, meshFace, initialStressInFaultCS, mu, slipRate1, slipRate2, RF);

      }//lts-face loop
      layerLtsFaceToMeshFace += it->getNumberOfCells();
    }//leaf_iterator loop
    std::cout << "init DR for Base\n";
  }

};

class seissol::initializers::Init_FL_2 : public seissol::initializers::BaseDrInitializer {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
        initializers::LTSTree* dynRupTree,
        std::unordered_map<std::string,
        double*> faultParameters,
        unsigned* ltsFaceToMeshFace,
        seissol::Interoperability &e_interoperability) override {
    BaseDrInitializer::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters, ltsFaceToMeshFace, e_interoperability);
    seissol::initializers::DR_FL_2 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 *>(dynRup);

    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
      real (*d_c)[numOfPointsPadded]                       = it->var(ConcreteLts->d_c);                 //from faultParameters
      real (*mu_S)[numOfPointsPadded]                      = it->var(ConcreteLts->mu_S);                //from faultParameters
      real (*mu_D)[numOfPointsPadded]                      = it->var(ConcreteLts->mu_D);                //from faultParameters
      real (*forced_rupture_time)[numOfPointsPadded]       = it->var(ConcreteLts->forced_rupture_time); //from faultParameters
      bool (*DS)[numOfPointsPadded]                        = it->var(ConcreteLts->DS);                  //from parameter file
      real *averaged_Slip                                  = it->var(ConcreteLts->averaged_Slip);       // = 0
      real (*dynStress_time)[numOfPointsPadded]            = it->var(ConcreteLts->dynStress_time);      // = 0

      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
        for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
          d_c[ltsFace][iBndGP]                    = static_cast<real>( faultParameters["d_c"][meshFace * numberOfPoints] );
          mu_S[ltsFace][iBndGP]                   = static_cast<real>( faultParameters["mu_s"][meshFace * numberOfPoints] );
          mu_D[ltsFace][iBndGP]                   = static_cast<real>( faultParameters["mu_d"][meshFace * numberOfPoints] );
          if(faultParameters["forced_rupture_time"] != NULL ){
              forced_rupture_time[ltsFace][iBndGP]    = static_cast<real>( faultParameters["forced_rupture_time"][meshFace * numberOfPoints] );
          }else{
              forced_rupture_time[ltsFace][iBndGP]    = 0.0;
          }
        }
        //initialize padded elements for vectorization
        for (unsigned iBndGP = numberOfPoints; iBndGP < numOfPointsPadded; ++iBndGP) {
          d_c[ltsFace][iBndGP]                    = 0.0;
          mu_S[ltsFace][iBndGP]                   = 0.0;
          mu_D[ltsFace][iBndGP]                   = 0.0;
          forced_rupture_time[ltsFace][iBndGP]    = 0.0;
        }
        averaged_Slip[ltsFace]= 0.0;

        for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded; ++iBndGP) {    //loop includes padded elements
          dynStress_time[ltsFace][iBndGP] = 0.0;
          DS[ltsFace][iBndGP] = m_Params.IsDsOutputOn;
        }
      }//lts-face loop
      layerLtsFaceToMeshFace += it->getNumberOfCells();
    }//leaf_iterator loop
  }
};


class seissol::initializers::Init_FL_3 : public seissol::initializers::BaseDrInitializer {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                          initializers::LTSTree* dynRupTree,
                                          std::unordered_map<std::string,
                                              double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability &e_interoperability) override {
    BaseDrInitializer::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters, ltsFaceToMeshFace, e_interoperability);
    seissol::initializers::DR_FL_3 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_3 *>(dynRup);
    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
      real *RS_f0                                               = it->var(ConcreteLts->RS_f0);
      real *RS_a                                                = it->var(ConcreteLts->RS_a);
      real *RS_b                                                = it->var(ConcreteLts->RS_b);
      real *RS_sl0                                              = it->var(ConcreteLts->RS_sl0);
      real *RS_sr0                                              = it->var(ConcreteLts->RS_sr0);
      real (*StateVar)[numOfPointsPadded]                       = it->var(ConcreteLts->StateVar);

      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

        //get initial values from fortran
        //TODO: write this function_
        e_interoperability.getDynRupFL_3(ltsFace, meshFace, RS_f0, RS_a, RS_b, RS_sl0, RS_sr0, StateVar);

      }//lts-face loop
      layerLtsFaceToMeshFace += it->getNumberOfCells();
    }//leaf_iterator loop
  }
};



class seissol::initializers::Init_FL_33 : public seissol::initializers::BaseDrInitializer {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
          initializers::LTSTree* dynRupTree,
          std::unordered_map<std::string,
          double*> faultParameters,
          unsigned* ltsFaceToMeshFace,
          seissol::Interoperability &e_interoperability) override {
    BaseDrInitializer::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters, ltsFaceToMeshFace, e_interoperability);
    seissol::initializers::DR_FL_33 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 *>(dynRup);

    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {

      real  (*nucleationStressInFaultCS)[numOfPointsPadded][6]  = it->var(ConcreteLts->nucleationStressInFaultCS); //get from fortran
      real *averaged_Slip                                       = it->var(ConcreteLts->averaged_Slip);      // = 0

      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

        e_interoperability.getDynRupNucStress(ltsFace, meshFace, nucleationStressInFaultCS);
        averaged_Slip[ltsFace]= 0.0;

      }//lts-face loop
      layerLtsFaceToMeshFace += it->getNumberOfCells();
    }//leaf_iterator loop
  }
};

class seissol::initializers::Init_FL_103 : public seissol::initializers::BaseDrInitializer {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                          initializers::LTSTree* dynRupTree,
                                          std::unordered_map<std::string,
                                              double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability &e_interoperability) override {
    BaseDrInitializer::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters, ltsFaceToMeshFace, e_interoperability);
    seissol::initializers::DR_FL_103 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103 *>(dynRup);

    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {


      real  (*nucleationStressInFaultCS)[numOfPointsPadded][6]  = it->var(ConcreteLts->nucleationStressInFaultCS); //get from fortran
      real (*RS_sl0_array)[numOfPointsPadded]                   = it->var(ConcreteLts->RS_sl0_array);       //get from faultParameters
      real (*RS_a_array)[numOfPointsPadded]                     = it->var(ConcreteLts->RS_a_array);         //get from faultParameters
      real (*RS_srW_array)[numOfPointsPadded]                   = it->var(ConcreteLts->RS_srW_array);       //get from faultParameters
      bool (*DS)[numOfPointsPadded]                             = it->var(ConcreteLts->DS);                 //par file
      real *averaged_Slip                                       = it->var(ConcreteLts->averaged_Slip);      // = 0
      real (*stateVar)[numOfPointsPadded]                       = it->var(ConcreteLts->stateVar);           //get from Fortran = EQN%IniStateVar
      real (*dynStress_time)[numOfPointsPadded]                 = it->var(ConcreteLts->dynStress_time);     // = 0

      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

        e_interoperability.getDynRupStateVar(ltsFace, meshFace, stateVar);
        e_interoperability.getDynRupNucStress(ltsFace, meshFace, nucleationStressInFaultCS);

        for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded; ++iBndGP) {    //loop includes padded elements
          dynStress_time[ltsFace][iBndGP] = 0.0;
          DS[ltsFace][iBndGP] = m_Params.IsDsOutputOn;
        }
        averaged_Slip[ltsFace]= 0.0;

        for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
          RS_a_array[ltsFace][iBndGP] = static_cast<real>( faultParameters["rs_a"][meshFace * numberOfPoints] );
          RS_srW_array[ltsFace][iBndGP] = static_cast<real>( faultParameters["rs_srW"][meshFace * numberOfPoints] );
          RS_sl0_array[ltsFace][iBndGP] = static_cast<real>( faultParameters["RS_sl0"][meshFace * numberOfPoints] );
        }
        //initialize padded elements for vectorization
        for (unsigned iBndGP = numberOfPoints; iBndGP < numOfPointsPadded; ++iBndGP) {
          RS_a_array[ltsFace][iBndGP] = 0.0;
          RS_srW_array[ltsFace][iBndGP] = 0.0;
          RS_sl0_array[ltsFace][iBndGP] = 0.0;
        }

      }//lts-face loop
      layerLtsFaceToMeshFace += it->getNumberOfCells();
    }//leaf_iterator loop

  }
};

class seissol::initializers::Init_FL_103_Thermal : public seissol::initializers::Init_FL_103 {
public:
  static constexpr unsigned int TP_grid_nz = 60;  //todo: make this global?

  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                          initializers::LTSTree* dynRupTree,
                                          std::unordered_map<std::string,
                                              double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability &e_interoperability) override {
    BaseDrInitializer::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters, ltsFaceToMeshFace, e_interoperability);
    Init_FL_103::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters, ltsFaceToMeshFace, e_interoperability);

    seissol::initializers::DR_FL_103_Thermal *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103_Thermal *>(dynRup);

    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;



    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {


      real (*TP)[numOfPointsPadded][2]                          = it->var(ConcreteLts->TP);
      real (*TP_Theta)[numOfPointsPadded][TP_grid_nz]           = it->var(ConcreteLts->TP_Theta);
      real (*TP_sigma)[numOfPointsPadded][TP_grid_nz]           = it->var(ConcreteLts->TP_sigma);

      real (*TP_half_width_shear_zone)[numOfPointsPadded]       = it->var(ConcreteLts->TP_half_width_shear_zone);
      real (*alpha_hy)[numOfPointsPadded]                       = it->var(ConcreteLts->alpha_hy);


      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

        for (unsigned iBndGP = numberOfPoints; iBndGP < numOfPointsPadded; ++iBndGP) {
          TP[ltsFace][iBndGP][0] = m_Params.IniTemp;
          TP[ltsFace][iBndGP][1] = m_Params.IniPressure;
          TP_half_width_shear_zone[ltsFace][iBndGP] = static_cast<real>( faultParameters["TP_half_width_shear_zone"][meshFace * numberOfPoints] );
          alpha_hy[ltsFace][iBndGP] = static_cast<real>( faultParameters["alpha_hy"][meshFace * numberOfPoints] );
          for (unsigned iTP_grid_nz = TP_grid_nz; iTP_grid_nz < TP_grid_nz; ++iTP_grid_nz) {
            TP_Theta[ltsFace][iBndGP][iTP_grid_nz] = 0.0;
            TP_sigma[ltsFace][iBndGP][iTP_grid_nz] = 0.0;
          }
        }
        //TODO: initialize all TPs

      }//lts-face loop
      layerLtsFaceToMeshFace += it->getNumberOfCells();
    }//leaf_iterator loop

  }
};

#endif //SEISSOL_DR_INITIALIZER_BASE_H
