
#ifndef _CRCVM_
#define _CRCVM_

#include "make_RACVM_matrix.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla; 
using namespace density; 
using namespace Eigen; 


//' Penalised negative log-likelihood for CRCVM model
template <class Type>
Type nllk_crcvm(objective_function<Type>* obj) { 
    //======//
    // DATA //
    //======//
    DATA_VECTOR(ID); // Time series ID
    DATA_VECTOR(times); // Observation times
    DATA_VECTOR(theta); // Observed angles
    DATA_VECTOR(DistanceShore); // Observed distances to boundary
    DATA_MATRIX(obs); // Response variables
    DATA_SPARSE_MATRIX(X_fe); // Design matrix for fixed effects
    DATA_SPARSE_MATRIX(X_re); // Design matrix for random effects
    DATA_SPARSE_MATRIX(S); // Penalty matrix
    DATA_IMATRIX(ncol_re);

    // Number of observations
    int n = obs.rows();
    
    // Time intervals
    vector<Type> dtimes = diff(times);
    
    //============//
    // PARAMETERS //
    //============//
    PARAMETER_VECTOR(coeff_fe); // Fixed effect parameters
    PARAMETER_VECTOR(log_lambda); // Smoothness parameters
    PARAMETER_VECTOR(coeff_re); // Random effect parameters
    

    // Derived parameters (linear predictors)
    vector<Type> par_vec = X_fe * coeff_fe + X_re * coeff_re;
    matrix<Type> par_mat(n, par_vec.size()/n);
    for(int i = 0; i < par_mat.cols(); i++) {
        // Matrix with one row for each time step and
        // one column for each parameter
        par_mat.col(i) = par_vec.segment(i*n, n);
    }

    // Parameters of velocity process
    vector<Type> tau = exp(par_mat.col(0).array());
    vector<Type> nu = exp(par_mat.col(1).array());
    vector<Type> a = exp(par_mat.col(2).array());
    vector<Type> b = exp(par_mat.col(3).array());
    vector<Type> D0 = exp(par_mat.col(4).array());
    vector<Type> D1 = exp(par_mat.col(5).array());
    vector<Type> sigma_D = exp(par_mat.col(6).array());
    vector<Type> sigma_theta = exp(par_mat.col(7).array());
    vector<Type> beta = 1/tau;
    vector<Type> sigma = 2 * nu / sqrt(M_PI * tau);
    vector<Type> omega= a*(theta-M_PI/2)*(theta+M_PI/2)*exp(-DistanceShore*DistanceShore/D0/D0)+
                        b*(exp(-1/2*((theta+M_PI/2/sqrt(3))*(theta+M_PI/2/sqrt(3))/sigma_theta/sigma_theta
                        +(DistanceShore-D1)*(DistanceShore-D1)/sigma_D/sigma_D))-
           exp(-1/2*((theta-M_PI/2/sqrt(3))*(theta-M_PI/2/sqrt(3))/sigma_theta/sigma_theta+(DistanceShore-D1)*(DistanceShore-D1)/sigma_D/sigma_D)));
    
    //============//
    // Likelihood //
    //============//
    // Initialise log-likelihood
    Type llk = 0;
    vector<Type> llk_all(n);
    llk_all.setZero();
    // Loop over observations
    for(int i = 1; i < n; i ++) {
        // No contribution if first observation of the track
        if(ID(i-1) == ID(i)) {

            matrix<Type> T = makeT_racvm(beta(i-1),omega(i-1),dtimes(i-1));
            matrix<Type> Q = makeQ_racvm(beta(i-1), sigma(i-1), omega(i-1), dtimes(i-1));

            vector<Type> u =  obs.row(i-1).transpose();     // Last observation 
            vector<Type> v =  obs.row(i).transpose();     // Current observation

            vector<Type> diff = v - T*u;

            // Compute the Mahalanobis term: diff' * Q^-1 * diff
            matrix<Type> Qinv = Q.inverse();
            vector<Type> Qinvdiff = Qinv*diff;
            Type mahalanobis = (diff* Qinvdiff).sum();

            // Compute log(det(Q))
            Type detQ = det(Q);

            // Multivariate normal log-likelihood
            llk_all(i) = -0.5 * (mahalanobis + log(detQ) + obs.cols() * log(2 * M_PI));
            
            llk = llk + llk_all(i);
        }
    }
    
 //===================//
    // Smoothing penalty //
    // ===================//
    Type nllk = -llk;
    // Are there random effects?
       if(ncol_re(0,0) > -1) {
        // Index in matrix S
        int S_start = 0;
        
        // Loop over smooths
        for(int i = 0; i < ncol_re.cols(); i++) {
            // Size of penalty matrix for this smooth
            int Sn = ncol_re(1,i) - ncol_re(0,i) + 1;
            
            // Penalty matrix for this smooth
            Eigen::SparseMatrix<Type> this_S = S.block(S_start, S_start, Sn, Sn);
            
            // Coefficients for this smooth
            vector<Type> this_coeff_re = coeff_re.segment(ncol_re(0,i) - 1, Sn);
            
            // Add penalty
            nllk = nllk -
                Type(0.5) * Sn * log_lambda(i) +
                Type(0.5) * exp(log_lambda(i)) * 
                density::GMRF(this_S).Quadform(this_coeff_re);
            
            // Increase index
            S_start = S_start + Sn;
        }
    }  
    
    return nllk;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif