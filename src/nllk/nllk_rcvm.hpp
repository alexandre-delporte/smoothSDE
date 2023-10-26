#ifndef _RCVM_
#define _RCVM_

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace R_inla;
using namespace density;
using namespace Eigen;

//' Matrix determinant
template<class Type>
Type det(matrix<Type> M) {
    int n_dim = M.cols();
    Type det = 0;
    if(n_dim == 1) {
        det = M(0, 0);
    } else if(n_dim == 2) {
        det = M(0,0) * M(1,1) - M(1,0) * M(0,1);
    } else {
        det = exp(atomic::logdet(M));
    }
    return det;
}

//' Make H matrix for Kalman filter
//'
//' @param sigma_obs SD of measurement error
//' @param n_dim Number of dimensions
template<class Type>
matrix<Type> makeH_ctcrw(Type sigma_obs, int n_dim) {
    matrix<Type> H(n_dim, n_dim);
    H.setZero();
    for(int i = 0; i < n_dim; i ++) {
        H(i, i) = sigma_obs * sigma_obs;
    }
    return H;
}

//' Make T matrix for Kalman filter
//'
//' @param beta Parameter beta of RCVM
//' @param omega Parameter omega of RCVM
//' @param dt Length of time interval
template<class Type>
matrix<Type> makeT_rcvm(Type beta, Type omega, Type delta) {

    matrix<Type> T(4, 4);
    T.setZero();

    matrix<Type> A(2,2);
    A << beta,-omega,omega,beta;

    double C = beta^2+omega^2;
    matrix<Type> invA(2,2);
    invA <<beta/C,omega/C,-omega/C,beta/C;

    matrix<Type> I(2,2);
    I << 1,0,0,1;

    matrix<Type> O(2,2);
    O.setZero();

    matrix<Type> R(2,2);
    R << cos(omega*delta),sin(omega*delta),-sin(omega*delta),cos(omega*delta);
    matrix<Type> expAdelta(2,2);
    expAdelta<<exp(-beta*delta)*R;

    // Combine the matrices into T
    T.block<2, 2>(0, 0) = I; // Top-left block
    T.block<2, 2>(0, 2) = invA*(I-expAdelta); // Top-right block
    T.block<2, 2>(2, 0) = O; // Bottom-left block
    T.block<2, 2>(2, 2) = expAdelta; // Bottom-right block

    return T;
}

//' Make Q matrix for Kalman filter
//'
//' @param beta Parameter beta of RCVM
//' @param sigma Parameter sigma of RCVM
//' @param omega Parameter omega od RCVM
//' @param delta Length of time interval
template<class Type>
matrix<Type> makeQ_rcvm(Type beta, Type sigma,Type omega, Type delta) {

    //initialize matrix with zeros
    matrix<Type> Q(4, 4);
    Q.setZero();

    //constants
    double tau=1/beta;
    double C=beta^2+omega^2;

    // variances and covariances values
    double var_xi=sigma^2/C*(delta+(omega^2-3/tau^2)/(2/tau*C)-exp(-2*delta/tau)/(2/tau)+
                        2*exp(-delta/tau)*(1/tau*cos(omega*delta)-omega*sin(omega*delta))/C);
    double var_zeta=sigma^2*tau/2*(1-exp(-2*delta/tau));

    double cov1=sigma^2/(2*C)*(1+exp(-2*delta/tau)-2*exp(-delta/tau)*cos(omega*delta));

    double cov2=sigma^2/C*(exp(-delta/tau)*sin(omega*delta)-omega/(2/tau)*(1-exp(-2*delta/tau)));

    // diagonal elements
    Q(0,0)=var_xi;
    Q(1,1)=var_xi;
    Q(2,2)=var_zeta;
    Q(3,3)=var_zeta;

    //zero off diagonal elements
    Q(0,1)=0;
    Q(1,0)=0;
    Q(2,3)=0;
    Q(3,2)=0;

    //non-zero off diagonal elements
    Q(0,2)=cov1;
    Q(2,0)=cov1;
    Q(0,3)=cov2;
    Q(3,0)=cov2;
    Q(1,2)=cov2;
    Q(2,1)=cov2;
    Q(1,3)=cov1;
    Q(3,1)=cov1;

    return Q;
}


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
Type nllk_rcvm(objective_function<Type>* obj) {
    //======//
    // DATA //
    //======//
    DATA_VECTOR(ID); // Time series ID
    DATA_VECTOR(times); // Observation times
    DATA_MATRIX(obs); // Response variables
    DATA_SPARSE_MATRIX(X_fe); // Design matrix for fixed effects
    DATA_SPARSE_MATRIX(X_re); // Design matrix for random effects
    DATA_SPARSE_MATRIX(S); // Penalty matrix
    DATA_IVECTOR(ncol_re); // Number of columns of S and X_re for each random effect
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
    vector<Type> tau = exp(par_mat.col(0).array());
    vector<Type> nu = exp(par_mat.col(1).array());
    vector<Type> omega = par_mat.col(2).array();
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
    matrix<Type> H = makeH_ctcrw(sigma_obs, n_dim);
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
            matrix<Type> T = makeT_rcvm(beta(i),omega(i),dtimes(i));
            matrix<Type> Q = makeQ_rcvm(beta(i), sigma(i), omega(i), dtimes(i));


            if(R_IsNA(asDouble(obs(i,0)))) {
                // If missing observation
                aest = T * aest;
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
                    aest = T * aest + K * u;
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
    if(ncol_re(0) > 0) {
        // Index in matrix S
        int S_start = 0;

        // Loop over smooths
        for(int i = 0; i < ncol_re.size(); i++) {
            // Size of penalty matrix for this smooth
            int Sn = ncol_re(i);

            // Penalty matrix for this smooth
            Eigen::SparseMatrix<Type> this_S = S.block(S_start, S_start, Sn, Sn);

            // Coefficients for this smooth
            vector<Type> this_coeff_re = coeff_re.segment(S_start, Sn);

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
