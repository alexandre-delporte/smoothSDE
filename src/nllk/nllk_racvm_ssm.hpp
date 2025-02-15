#ifndef _RACVM_SSM_
#define _RACVM_SSM_

#include "utility.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla;
using namespace density;
using namespace Eigen;


//' Penalised negative log-likelihood for CTCRW
//'
//' This function was inspired by the source code of the package
//' crawl, authored by Devin Johnson and Josh London
//'
//' All derivations, including for the matrices T and Q defined above, are detailed
//' in Section 6.2.2 of Michelot (2019), Stochastic models of animal movement and
//' habitat selection. PhD thesis, University of Sheffield.
//' (etheses.whiterose.ac.uk/23688/1/TheoMichelot_PhD_thesis_April2019.pdf)
template <class Type>
Type nllk_racvm_ssm(objective_function<Type>* obj) {
    //======//
    // DATA //
    //======//
    DATA_VECTOR(ID); // Time series ID
    DATA_VECTOR(times); // Observation times
    DATA_MATRIX(obs); // Response variables
    DATA_SPARSE_MATRIX(X_fe); // Design matrix for fixed effects
    DATA_SPARSE_MATRIX(X_re); // Design matrix for random effects
    DATA_SPARSE_MATRIX(S); // Penalty matrix
    DATA_IMATRIX(ncol_re); // Start and end indexes of S and X_re for each random effect
    DATA_MATRIX(a0); // Initial state estimate for Kalman filter
    DATA_MATRIX(P0); // Initial state covariance for Kalman filter

    DATA_ARRAY(H_array); // Covariance matrices for observation error

    // Number of observations
    int n = obs.rows();

    // Number of dimensions
    int n_dim = obs.cols();

    // Time intervals (needs to be of length n)
    vector<Type> dtimes(n);
    for(int i = 0; i < n-1; i++)
        dtimes(i) = times(i+1) - times(i);
    dtimes(n-1) = 1;

    //============//
    // PARAMETERS //
    //============//
    // SD of measurement error
    PARAMETER(log_sigma_obs);
    Type sigma_obs = exp(log_sigma_obs);

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
    matrix<Type> mu = par_mat.block(0, 0, n, 2).array();
    vector<Type> tau = exp(par_mat.col(2).array());
    vector<Type> nu = exp(par_mat.col(3).array());
    vector<Type> omega = par_mat.col(4).array();
    vector<Type> beta = 1/tau;
    vector<Type> sigma = 2 * nu / sqrt(M_PI * tau);

    //================================//
    // Likelihood using Kalman filter //
    //================================//
    // Define all matrices and vectors needed below
    matrix<Type> Z(2, 4);
    Z.setZero();
    Z(0,0)=1;
    Z(1,1)=1;
    matrix<Type> H = makeH(sigma_obs,n_dim);
    matrix<Type> T(4, 4);
    matrix<Type> Q(4, 4);
    matrix<Type> F(2, 2);
    F.setZero();
    matrix<Type> K(4, 2);
    K.setZero();
    matrix<Type> L(4, 4);
    L.setZero();
    vector<Type> u(2);
    u.setZero();
    Type detF;

    // Initial state mean
    vector<Type> aest(4);
    aest = a0.row(0);

    // Initial state covariance matrix
    matrix<Type> Pest(4, 4);
    Pest = P0;

    // Counter for ID (to initialise a0)
    int k = 1;

    // Kalman filter iterations
    Type llk = 0;
    matrix<Type> aest_all(n, 4);
    aest_all.setZero();
    aest_all.row(0) = aest;
    for(int i = 1; i < n; i++) {
        if(ID(i) != ID(i-1)) {
            // If first location of track, re-initialise state vector

            aest = a0.row(k);
            k = k + 1;
            Pest = P0;
        } else {
            // Compute Kalman filter matrices
            if(H_array.size() > 1) {
                H = H_array.col(i).matrix();
            }
            matrix<Type> T = makeT_racvm(beta(i),omega(i),dtimes(i));
            matrix<Type> Q = makeQ_racvm(beta(i), sigma(i), omega(i), dtimes(i));
            matrix<Type> B = makeB_racvm(beta(i),omega(i),dtimes(i));

            // Mean velocity component of state update
            vector<Type> mu_i = mu.row(i).transpose();
            vector<Type> B_times_mu = B * mu_i;

            if(R_IsNA(asDouble(obs(i,0)))) {
                // If missing observation
                aest = T * aest + B_times_mu;
                Pest = T * Pest * T.transpose() + Q;
            } else {
                // Measurement residual
                vector<Type> obsrow =  obs.row(i).transpose();
                u = obsrow - Z * aest;
                // Residual covariance
                F = Z * Pest * Z.transpose() + H;
                detF = det(F);

                if(detF <= 0) {
                    aest = T * aest;
                    Pest = T * Pest * T.transpose() + Q;
                } else {
                    // Update log-likelihood
                    matrix<Type> FinvT = F.inverse().transpose();
                    vector<Type> FinvTu = FinvT * u;
                    Type uFu = (u * FinvTu).sum();
                    llk = llk - (log(detF) + uFu)/2;
                    // Kalman gain
                    K = T * Pest * Z.transpose() * F.inverse();
                    // Update state estimate
                    aest = T * aest + K * u + B_times_mu;
                    // Update estimate covariance
                    L = T - K * Z;
                    Pest = T * Pest * L.transpose() + Q;
                }
            }
        }

        aest_all.row(i) = aest;
    }

    REPORT(aest_all)

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
