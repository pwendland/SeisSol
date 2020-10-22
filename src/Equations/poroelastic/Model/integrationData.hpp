#ifndef MODEL_POROELASTIC_INTEGRATIONDATA_H_
#define MODEL_POROELASTIC_INTEGRATIONDATA_H_

#include <generated_code/tensor.h>

namespace seissol {
  namespace model {

    struct PoroelasticLocalData {
      real sourceMatrix[seissol::tensor::ET::size()];
    };
    struct PoroelasticNeighborData {};
  }
}

#endif
