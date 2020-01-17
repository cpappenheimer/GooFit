#pragma once

#include <vector>

#include "goofit/GlobalCudaDefines.h"

namespace GooFit {

class MathUtils final {
  public:
  static fptype doNeumaierSummation(const std::vector<fptype>& vals);
};

} // end namespace GooFit
