#include "goofit/PDFs/basic/BifurGaussPdf.h"

namespace GooFit {

__device__ fptype device_BifurGauss(fptype* evt, ParameterContainer &pc) {
    fptype x = evt[0]; // why does indices recall itself?
    fptype mean = pc.parameters[1];
    fptype sigmaLeft = pc.parameters[2];
    fptype sigmaRight = pc.parameters[3];

    // how to calculate the value of a bifurcated gaussian?
    fptype sigma = sigmaLeft;

    if(x > mean)
        sigma = sigmaRight;

    pc.incrementIndex(1, 3, 0, 0, 1);

    fptype ret = exp(-0.5*(x-mean)*(x-mean)/(sigma*sigma));
    return ret;
}

__device__ device_function_ptr ptr_to_BifurGauss = device_BifurGauss;

__host__ BifurGaussPdf::BifurGaussPdf(std::string n, Variable *_x, Variable *mean, Variable *sigmaL, Variable *sigmaR)
    : GooPdf(_x, n) {
    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(sigmaL));
    pindices.push_back(registerParameter(sigmaR));
    GET_FUNCTION_ADDR(ptr_to_BifurGauss);
    initialize(pindices);
}

__host__ void BifurGaussPdf::recursiveSetIndices () {
    GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_BifurGauss");
    GET_FUNCTION_ADDR(ptr_to_BifurGauss);

    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions ++;

    populateArrays ();
}

// q: how shall the normalization of a bifurcated gaussian be calculated?
// a: a "sum" of two half-gaussians?
__host__ fptype BifurGaussPdf::integrate(fptype lo, fptype hi) const {

    fptype sL = host_parameters[parametersIdx + 2];
    fptype sR = host_parameters[parametersIdx + 3];

    fptype normL = 1. / (sqrt(2 * M_PI) * sL);
    fptype normR = 1. / (sqrt(2 * M_PI) * sR);

    return .5 * normL + .5 * normR;
}
} // namespace GooFit
