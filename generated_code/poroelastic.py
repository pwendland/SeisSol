#!/usr/bin/env python3
  
import numpy as np
from yateto import Tensor, Scalar, simpleParameterSpace
from yateto.input import parseXMLMatrixFile, parseJSONMatrixFile, memoryLayoutFromFile
from yateto.ast.node import Add
from yateto.ast.transformer import DeduceIndices, EquivalentSparsityPattern

from elastic import ElasticADERDG as ADERDGBase
from multSim import OptionalDimTensor

class PoroelasticADERDG(ADERDGBase):
  def __init__(self, order, multipleSimulations, matricesDir, memLayout, **kwargs):
    super().__init__(order, multipleSimulations, matricesDir, memLayout)
    clones = {
      'star': ['star(0)', 'star(1)', 'star(2)'],
    }
    self.db.update( parseXMLMatrixFile('{}/matrices_poroelastic.xml'.format(matricesDir), clones) )
    memoryLayoutFromFile(memLayout, self.db, clones)
    ET_spp = self.db['ET'].spp().as_ndarray()
    self.db.ET = Tensor('ET', ET_spp.shape, spp=ET_spp)

  def sourceMatrix(self):
    return self.db.ET

  def numberOfQuantities(self):
    return 13

  def godunov_spp(self):
    spp = np.ones((self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool)
    return spp

  def flux_solver_spp(self):
    return self.godunov_spp()

  def transformation_spp(self):
    spp = np.ones((self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool)
    return spp

  def transformation_inv_spp(self):
    return self.transformation_spp()

  def addInit(self, generator):
      super().addInit(generator)

  def add_include_tensors(self, include_tensors):
      super().add_include_tensors(include_tensors)
