//
// Created by Ravil
//

#include "tree/Common.hpp"
#include "LTS.h"
#include "tree/LTSTree.hpp"

// TODO: Memkind
using namespace seissol::initializers;
  void LTS::addTo(LTSTree& tree) {
  #ifdef USE_PLASTICITY
      LayerMask plasticityMask = LayerMask(Ghost);
  #else
      LayerMask plasticityMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
  #endif
      tree.addVar(                    dofs, LayerMask(Ghost),     PAGESIZE_HEAP,      MEMKIND_DOFS );
  #if NUMBER_OF_RELAXATION_MECHANISMS > 0
      tree.addVar(                 dofsAne, LayerMask(Ghost),     PAGESIZE_HEAP,      MEMKIND_DOFS );
  #endif
      tree.addVar(                 buffers,      LayerMask(),                 1,      MEMKIND_TIMEDOFS );
      tree.addVar(             derivatives,      LayerMask(),                 1,      MEMKIND_TIMEDOFS );
      tree.addVar(         cellInformation,      LayerMask(),                 1,      MEMKIND_CONSTANT );
      tree.addVar(           faceNeighbors, LayerMask(Ghost),                 1,      MEMKIND_TIMEDOFS );
      tree.addVar(        localIntegration, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
      tree.addVar(  neighboringIntegration, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
      tree.addVar(                material, LayerMask(Ghost),                 1,      seissol::memory::Standard );
      tree.addVar(              plasticity,   plasticityMask,                 1,      seissol::memory::Standard );
      tree.addVar(               drMapping, LayerMask(Ghost),                 1,      MEMKIND_CONSTANT );
      tree.addVar(                  energy,   plasticityMask,     PAGESIZE_HEAP,      seissol::memory::Standard );
      tree.addVar(                 pstrain,   plasticityMask,     PAGESIZE_HEAP,      seissol::memory::Standard );
      tree.addVar(           displacements, LayerMask(Ghost),     PAGESIZE_HEAP,      seissol::memory::Standard );

      tree.addBucket(buffersDerivatives,                          PAGESIZE_HEAP,      MEMKIND_TIMEDOFS );
      tree.addBucket(displacementsBuffer,                         PAGESIZE_HEAP,      MEMKIND_TIMEDOFS );
  }