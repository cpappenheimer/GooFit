#include "goofit/PdfBase.h"
#include "goofit/fitting/FCN.h"
#include "goofit/PDFs/GooPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

FCN::FCN(Params& params) : params_(&params) {
    host_callnumber = 0;
    
}

double FCN::operator()(const std::vector<double>& pars) const {
    std::vector<double> gooPars; // Translates from Minuit indexing to GooFit indexing
    
    for(Variable* var : params_->vars_) {
        var->unchanged_ = var->value == pars.at(var->index);
        gooPars.push_back(pars.at(var->index));
    }

    params_->pdf_->copyParams(gooPars);
    double nll = params_->pdf_->calculateNLL();
    host_callnumber++;

    return nll;
}

// Get the number of variable parameters
Params* FCN::GetParams() {
    return params_;
}
    
}
