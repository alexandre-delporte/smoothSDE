
#ifndef _CRCVM_SSM_
#define _CRCVM_SSM_

#include "utility.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla;
using namespace density;
using namespace Eigen;


//' Penalised negative log-likelihood for CRCVM


template <class Type>
Type nllk_crcvm_ssm(objective_function<Type>* obj) {
    //======//
    // DATA //
    //======//
    DATA_VECTOR(ID); // Time series ID
    DATA_VECTOR(times); // // Observation times
    DATA_MATRIX(obs); // Response variables
    DATA_MATRIX(interpolated_BoundaryDistance); // Interpolated BoundaryDistance values
    DATA_MATRIX(interpolated_BoundaryAngle); // Interpolated BoundaryAngle values
    DATA_SPARSE_MATRIX(X_fe); // Design matrix for fixed effects
    DATA_SPARSE_MATRIX(X_re); // Design matrix for random effects
    DATA_SPARSE_MATRIX(S); // Penalty matrix
    DATA_IMATRIX(ncol_re); // Start and end indexes of S and X_re for each random effect
    DATA_MATRIX(a0); // Initial state estimate for Kalman filter
    DATA_MATRIX(P0); // Initial state covariance for Kalman filter

    DATA_ARRAY(H_array); // Covariance matrices for observation error

    // Number of observations
    int n = obs.rows();

    // Number of interpolations
    int m = interpolated_BoundaryDistance.cols();

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
         
    //================================//
    // Likelihood using Kalman filter //
    //================================//
    // Define all matrices and vectors needed below
    matrix<Type> Z(2, 4);
    Z.setZero();
    Z(0,0)=1;
    Z(1,1)=1;
    matrix<Type> H = makeH(sigma_obs,n_dim);
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

    matrix<Type> omega_all(n,m);
    omega_all.setZero();

    array<Type> T_matrices(4, 4, n);
    array<Type> Q_matrices(4, 4, n);


    matrix<Type> dtimes_all(n,m);
    dtimes_all.setZero();


    Type omega_size = 0;

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


            vector<Type> angle = interpolated_BoundaryAngle.row(i);
            vector<Type> distance = interpolated_BoundaryDistance.row(i);

            vector<Type> omega= a(i)*angle*(angle-0.5*M_PI)*(angle+0.5*M_PI)*exp(-distance/D0(i))/distance+
                b(i)*(exp(-0.5*((angle+0.5*M_PI/sqrt(3))*(angle+0.5*M_PI/sqrt(3))/sigma_theta(i)/sigma_theta(i)+(distance-D1(i))*(distance-D1(i))/sigma_D(i)/sigma_D(i)))-
                    exp(-0.5*((angle-0.5*M_PI/sqrt(3))*(angle-0.5*M_PI/sqrt(3))/sigma_theta(i)/sigma_theta(i)+(distance-D1(i))*(distance-D1(i))/sigma_D(i)/sigma_D(i))));

            omega_size=omega.size();
            omega_all.row(i) = omega;

            vector<Type> dtimes_i (m);

            // Fill the vector with values from 0 to dtimes(i)
            for (int j = 0 ; j < m; j++) {
                dtimes_i(j) = dtimes(i)/Type(m);
            }

            dtimes_all.row(i)=dtimes_i;

            matrix<Type> T = makeT_crcvm(beta(i),omega,dtimes_i);
            for (int r = 0; r < 4; r++) {
                for (int c = 0; c < 4; c++) {
                 T_matrices(r, c, i) = T(r,c);
                }
            }

            matrix<Type> Q = makeQ_crcvm(beta(i), sigma(i), omega, dtimes_i);

            for (int r = 0; r < 4; r++) {
                for (int c = 0; c < 4; c++) {
                 Q_matrices(r, c, i) = Q(r,c);
                }
            }

            
            if(R_IsNA(asDouble(obs(i,0)))) {
                // If missing observation
                aest = T * aest ;
                Pest = T * Pest * T.transpose() + Q;
            } else {
                // Measurement residual
                vector<Type> obsrow =  obs.row(i).transpose();
                u = obsrow - Z * aest;
                // Residual covariance
                F = Z * Pest * Z.transpose() + H;
                detF = det(F);

                if(detF <= 0) {
                    aest = T * aest ;
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
                    aest = T * aest + K * u ;
                    // Update estimate covariance
                    L = T - K * Z;
                    Pest = T * Pest * L.transpose() + Q;
                }
            }
        }

        aest_all.row(i) = aest;
    }

    REPORT(aest_all);
    REPORT(interpolated_BoundaryDistance); 
    REPORT(interpolated_BoundaryAngle);
    REPORT(omega_all);   
    REPORT(par_mat);
    REPORT(T_matrices);
    REPORT(Q_matrices);
    REPORT(dtimes_all);
    REPORT(omega_size);
   


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

    REPORT(nllk);  
    
    return nllk;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
