// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/independence.h"
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// max_taustar
double max_taustar();
static SEXP _independence_max_taustar_try() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    rcpp_result_gen = Rcpp::wrap(max_taustar());
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _independence_max_taustar() {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_independence_max_taustar_try());
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// max_hoeffding
double max_hoeffding();
static SEXP _independence_max_hoeffding_try() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    rcpp_result_gen = Rcpp::wrap(max_hoeffding());
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _independence_max_hoeffding() {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_independence_max_hoeffding_try());
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// calc_taustar
double calc_taustar(const std::vector<unsigned long>& perm);
static SEXP _independence_calc_taustar_try(SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned long>& >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_taustar(perm));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _independence_calc_taustar(SEXP permSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_independence_calc_taustar_try(permSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// calc_hoeffding
double calc_hoeffding(const std::vector<unsigned long>& perm);
static SEXP _independence_calc_hoeffding_try(SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned long>& >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_hoeffding(perm));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _independence_calc_hoeffding(SEXP permSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_independence_calc_hoeffding_try(permSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// calc_refined
double calc_refined(const std::vector<unsigned long>& perm);
static SEXP _independence_calc_refined_try(SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::vector<unsigned long>& >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_refined(perm));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _independence_calc_refined(SEXP permSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_independence_calc_refined_try(permSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _independence_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("double(*max_taustar)()");
        signatures.insert("double(*max_hoeffding)()");
        signatures.insert("double(*.calc.taustar)(const std::vector<unsigned long>&)");
        signatures.insert("double(*.calc.hoeffding)(const std::vector<unsigned long>&)");
        signatures.insert("double(*.calc.refined)(const std::vector<unsigned long>&)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _independence_RcppExport_registerCCallable() { 
    R_RegisterCCallable("independence", "_independence_max_taustar", (DL_FUNC)_independence_max_taustar_try);
    R_RegisterCCallable("independence", "_independence_max_hoeffding", (DL_FUNC)_independence_max_hoeffding_try);
    R_RegisterCCallable("independence", "_independence_.calc.taustar", (DL_FUNC)_independence_calc_taustar_try);
    R_RegisterCCallable("independence", "_independence_.calc.hoeffding", (DL_FUNC)_independence_calc_hoeffding_try);
    R_RegisterCCallable("independence", "_independence_.calc.refined", (DL_FUNC)_independence_calc_refined_try);
    R_RegisterCCallable("independence", "_independence_RcppExport_validate", (DL_FUNC)_independence_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_independence_max_taustar", (DL_FUNC) &_independence_max_taustar, 0},
    {"_independence_max_hoeffding", (DL_FUNC) &_independence_max_hoeffding, 0},
    {"_independence_calc_taustar", (DL_FUNC) &_independence_calc_taustar, 1},
    {"_independence_calc_hoeffding", (DL_FUNC) &_independence_calc_hoeffding, 1},
    {"_independence_calc_refined", (DL_FUNC) &_independence_calc_refined, 1},
    {"_independence_RcppExport_registerCCallable", (DL_FUNC) &_independence_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_independence(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
