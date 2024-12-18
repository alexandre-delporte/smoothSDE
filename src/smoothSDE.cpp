
#include <TMB.hpp>
#include "nllk/nllk_sde.hpp"
#include "nllk/nllk_bm_ssm.hpp"
#include "nllk/nllk_ou_ssm.hpp"
#include "nllk/nllk_ctcrw.hpp"
#include "nllk/nllk_e_seal_ssm.hpp"
#include "nllk/nllk_racvm.hpp"
#include "nllk/nllk_crcvm_cubic.hpp"
#include "nllk/nllk_crcvm_test.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
    // SDE type
    DATA_STRING(type);

    if (type == "BM" || type == "BM_t" || type == "OU" || type == "CIR") {
        return nllk_sde(this);
    } else if (type == "BM_SSM") {
        return(nllk_bm_ssm(this));
    } else if (type == "OU_SSM") {
        return(nllk_ou_ssm(this));
    } else if (type == "CTCRW") {
        return nllk_ctcrw(this);
    } else if (type == "ESEAL_SSM") {
        return nllk_eseal_ssm(this);
    } else if (type=="RACVM") {
        return nllk_racvm(this);
    } else if (type=="CRCVM_cubic") {
        return nllk_crcvm_cubic(this);
    } else if (type=="CRCVM_test") {
        return nllk_crcvm_test(this);
    }else {
        error ("Unknown SDE type");
    }
    return 0;
}
