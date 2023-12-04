
#include <TMB.hpp>
#include "nllk/nllk_sde.hpp"
#include "nllk/nllk_bm_ssm.hpp"
#include "nllk/nllk_ou_ssm.hpp"
#include "nllk/nllk_ctcrw.hpp"
#include "nllk/nllk_e_seal_ssm.hpp"
#include "nllk/nllk_racvm1.hpp"
#include "nllk/nllk_racvm2.hpp"
#include "nllk/nllk_crcvm1.hpp"
#include "nllk/nllk_crcvm2.hpp"

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
    } else if (type=="RACVM1") {
        return nllk_racvm1(this);
    } else if (type=="RACVM2") {
        return nllk_racvm2(this);
    } else if (type=="CRCVM1") {
        return nllk_crcvm1(this);
    } else if (type=="CRCVM2") {
        return nllk_crcvm2(this);
    }  else {
        error ("Unknown SDE type");
    }
    return 0;
}
