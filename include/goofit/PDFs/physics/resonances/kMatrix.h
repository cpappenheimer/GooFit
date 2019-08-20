#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

#include <array>

#define NPOLES 5
#define NCHANNELS 5

namespace GooFit {
namespace Resonances {

class kMatrix : public ResonancePdf {
  public:
    kMatrix(std::string name,
            unsigned int pterm, //< 0 or 1
            bool is_pole,       //< False for prod
	    Variable beta1_r,
	    Variable beta1_i,
            Variable sA0,
            Variable sA,
            Variable s0_prod,
            Variable s0_scatt,
            std::vector<Variable> fscat,
            std::vector<Variable> poles,
            Variable mass,
            Variable width,
            unsigned int L,
            unsigned int Mpair,
            FF FormFac    = FF::BL_Prime,
            fptype radius = 1.5);

    ~kMatrix() override = default;
};

} // namespace Resonances
} // namespace GooFit
