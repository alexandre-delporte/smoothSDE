#ifndef _OU_SSM_
#define _OU_SSM_

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla;
using namespace density;
using namespace Eigen;

//' Make H matrix for Kalman filter
//'
//' @param sigma_obs SD of measurement error
//' @param n_dim Number of dimensions
template<class Type>
matrix<Type> makeH_ou_ssm(Type sigma_obs, int n_dim) {
    matrix<Type> H(n_dim, n_dim);
    H.setZero();
    for(int i = 0; i < n_dim; i ++) {
        H(i, i) = sigma_obs * sigma_obs;
    }
    return H;
}

//' Make T matrix for Kalman filter
//'
//' @param tau Time scale parameter
//' @param dt Length of time interval
//' @param n_dim Number of dimensions
template<class Type>
matrix<Type> makeT_ou_ssm(Type tau, Type dt, int n_dim) {
    matrix<Type> T(n_dim, n_dim);
    T.setZero();
    for(int i = 0; i < n_dim; i ++) {
        T(i, i) = exp(-dt/tau);
    }
    return T;
}

//' Make B matrix for Kalman filter
//'
//' @param tau Time scale parameter
//' @param dt Length of time interval
//' @param n_dim Number of dimensions
template<class Type>
matrix<Type> makeB_ou_ssm(Type tau, Type dt, int n_dim) {
    matrix<Type> B(n_dim, n_dim);
    B.setZero();
    for(int i = 0; i < n_dim; i ++) {
        B(i, i) = 1 - exp(-dt/tau);
    }
    return B;
}

//' Make Q matrix for Kalman filter
//'
//' @param tau Time scale parameter
//' @param kappa Spatial scale parameter
//' @param dt Length of time interval
//' @param n_dim Number of dimensions
template<class Type>
matrix<Type> makeQ_ou_ssm(Type tau, Type kappa, Type dt, int n_dim) {
    matrix<Type> Q(n_dim, n_dim);
    Q.setZero();
    for(int i = 0; i < n_dim; i ++) {
        Q(i, i) = kappa * (1 - exp(-2*dt/tau));
    }
    return Q;
}

//' Penalised negative log-likelihood for Ornstein-Uhlenbeck
//' state-space model (with measurement error)
template <class Type>
Type nllk_ou_ssm(objective_function<Type>* obj) {
    //======//
    // DATA //
    //======//
    DATA_VECTOR(ID); // Time series ID
    DATA_VECTOR(times); // Observation times
    DATA_MATRIX(obs); // Response variables
    DATA_SPARSE_MATRIX(X_fe); // Design matrix for fixed effects
    DATA_SPARSE_MATRIX(X_re); // Design matrix for random effects
    DATA_SPARSE_MATRIX(S); // Penalty matrix
    DATA_IVECTOR(start_ncol_re); // Number of columns of S and X_re for each random effect
    DATA_IVECTOR(end_ncol_re);
    DATA_MATRIX(a0); // Initial state estimate for Kalman filter
    DATA_MATRIX(P0); // Initial state covariance for Kalman filter

    DATA_ARRAY(H_array); // Covariance matrices for observation error
    
    // Number of dimensions
    int n_dim = obs.cols();
    // Number of observations
    int n = obs.rows();

    // Time intervals
    vector<Type> dtimes(n);
    for(int i = 0; i < n-1; i++)
        dtimes(i) = times(i+1) - times(i);
    dtimes(n - 1) = 1;
    
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

    // Parameters of latent process
    matrix<Type> mu = par_mat.block(0, 0, n, n_dim).array();
    vector<Type> tau = exp(par_mat.col(n_dim).array());
    vector<Type> kappa = exp(par_mat.col(n_dim + 1).array());
    
    //================================//
    // Likelihood using Kalman filter //
    //================================//
    // Define all matrices and vectors needed below
    matrix<Type> Z(n_dim, n_dim);
    Z.setIdentity();
    matrix<Type> H = makeH_ou_ssm(sigma_obs, n_dim);
    matrix<Type> T(n_dim, n_dim);
    T.setIdentity();
    matrix<Type> B(n_dim, n_dim);
    T.setIdentity();
    matrix<Type> Q(n_dim, n_dim);
    matrix<Type> F(n_dim, n_dim);
    F.setZero();
    matrix<Type> K(n_dim, n_dim);
    K.setZero();
    matrix<Type> L(n_dim, n_dim);
    L.setZero();
    vector<Type> u(n_dim);
    u.setZero();
    Type detF = 1;

    // Initial state mean
    vector<Type> aest(n_dim);
    aest = a0.row(0);
    // Initial state covariance matrix
    matrix<Type> Pest(n_dim, n_dim);
    Pest = P0;

    // Counter for ID (to initialise a0)
    int k = 1;
    
    // Kalman filter iterations
    Type llk = 0;
    matrix<Type> aest_all(n, n_dim);
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
            matrix<Type> T = makeT_ou_ssm(tau(i), dtimes(i), n_dim);
            matrix<Type> B = makeB_ou_ssm(tau(i), dtimes(i), n_dim);
            matrix<Type> Q = makeQ_ou_ssm(tau(i), kappa(i), dtimes(i), n_dim);
            vector<Type> B_mu = B * mu.row(i).transpose();

            if(R_IsNA(asDouble(obs(i,0)))) {
                // If missing observation
                aest = T * aest + B_mu;
                Pest = T * Pest * T.transpose() + Q;
            } else {
                // Measurement residual
                vector<Type> obsrow =  obs.row(i).transpose();
                u = obsrow - Z * aest;
                
                // Residual covariance
                F = Z * Pest * Z.transpose() + H;
                detF = exp(atomic::logdet(F));

                if(detF<=0) {
                    aest = T * aest + B_mu;
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
                    aest = T * aest + K * u + B_mu;
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
       if(start_ncol_re(0) > 0) {
        // Index in matrix S
        int S_start = 0;
        
        // Loop over smooths
        for(int i = 0; i < start_ncol_re.size(); i++) {
            // Size of penalty matrix for this smooth
            int Sn = end_ncol_re(i) - start_ncol_re(i) + 1;
            
            // Penalty matrix for this smooth
            Eigen::SparseMatrix<Type> this_S = S.block(S_start, S_start, Sn, Sn);
            
            // Coefficients for this smooth
            vector<Type> this_coeff_re = coeff_re.segment(start_ncol_re(i) - 1, Sn);
            
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
