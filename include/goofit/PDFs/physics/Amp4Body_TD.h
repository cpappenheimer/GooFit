/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#pragma once

#include <goofit/PDFs/physics/Amp4BodyBase.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/PDFs/physics/detail/NormSpinCalculator_TD.h>

#include <mcbooster/GContainers.h>

#include <thrust/remove.h>

#include <tuple>

namespace GooFit {

class LSCalculator_TD;
class AmpCalc_TD;
class SFCalculator_TD;
class NormIntegrator_TD;
class Lineshape;

struct MCNormBatchResult {
  const int _numAcc;
  const fptype _normValue;
  MCNormBatchResult(int numAcc, fptype normValue)
  : _numAcc(numAcc),
    _normValue(normValue) {}
};

class Amp4Body_TD : public Amp4BodyBase {
  public:
    Amp4Body_TD(std::string n,
                std::vector<Observable> observables,
                DecayInfo4t decay,
                MixingTimeResolution *r,
                GooPdf *eff,
                Observable *mistag,
		std::vector<int> normSeeds,
                unsigned int mcEventsNormPerBatch);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalization will get *really* confused and give wrong answers.

    __host__ fptype normalize() override;
    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 8);
    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }
    __host__ void setGenerationOffset(int off) { generation_offset = off; }
    __host__ void setMaxWeight(fptype wmax) { maxWeight = wmax; }

    __host__ std::
        tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h>
        GenerateSig(unsigned int numEvents, int seed = 0);

    __host__ void populateArrays() override;

    __host__ int getNumAccNormEvents() const { return _numAccNormEvents; }

    static void  printDeviceVecComplexVals(thrust::device_vector<fpcomplex>::const_iterator first, thrust::device_vector<fpcomplex>::const_iterator last,
				     thrust::iterator_difference<thrust::device_vector<fpcomplex>::const_iterator>::type stride,
				     const std::string &codeLoc);

  protected:
  private:
    __host__ MCNormBatchResult computeNormBatch(unsigned int batchNum);
    __host__ void computeValues();

    std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<unsigned int>>> AmpMap;
    std::vector<SpinFactor *> SpinFactors;
    std::vector<Lineshape *> LineShapes;
    std::vector<AmpCalc_TD *> AmpCalcs;
    NormIntegrator_TD *Integrator;
    std::vector<SFCalculator_TD *> sfcalculators;
    std::vector<LSCalculator_TD *> lscalculators;
    unsigned int efficiencyFunction;
    DecayInfo4t decayInfo;
    MixingTimeResolution *resolution;
    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    thrust::device_vector<fpcomplex> *cachedResSF{nullptr}; // Caches the BW values and Spins for each event.
    thrust::device_vector<fpcomplex> *cachedAMPs{nullptr};  // cache Amplitude values for each event.
    mutable bool generation_no_norm{false};
    mutable bool SpinsCalculated{false};
    bool *redoIntegral;
    mutable bool forceRedoIntegrals{true};
    fptype *cachedMasses;
    fptype *cachedWidths;
    int totalEventSize;
    int cacheToUse{0};
    unsigned int generation_offset{0};
    double maxWeight{0};
    int _numAccNormEvents{0};

    const std::vector<int> _NORM_SEEDS;
    const unsigned int _MC_EVENTS_NORM_PER_BATCH;

    static const int _PRINT_LIMIT;
};

} // namespace GooFit
