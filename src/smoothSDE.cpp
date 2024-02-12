
#include <TMB.hpp>
#include "nllk/nllk_sde.hpp"
#include "nllk/nllk_bm_ssm.hpp"
#include "nllk/nllk_ou_ssm.hpp"
#include "nllk/nllk_ctcrw.hpp"
#include "nllk/nllk_e_seal_ssm.hpp"
#include "nllk/nllk_racvm1.hpp"
#include "nllk/nllk_racvm2.hpp"
#include "nllk/nllk_vm_crcvm.hpp"
#include "nllk/nllk_s_crcvm1.hpp"
#include "nllk/nllk_i_crcvm.hpp"
#include "nllk/nllk_s_crcvm2.hpp"
#include "nllk/nllk_si_crcvm.hpp"
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
    } else if (type=="RACVM3") {
        return nllk_racvm2(this);
    } else if (type=="VM_CRCVM") {
        return nllk_vm_crcvm(this);
    } else if (type=="S_CRCVM1") {
        return nllk_s_crcvm1(this);
    } else if (type=="I_CRCVM") {
        return nllk_i_crcvm(this);
    } else if (type=="S_CRCVM2") {
        return nllk_s_crcvm2(this); 
     } else if (type=="SI_CRCVM") {
        return nllk_si_crcvm(this); 
    } else {
        error ("Unknown SDE type");
    }
    return 0;
}
