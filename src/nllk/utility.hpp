#ifndef _utility_
#define _utility_

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

//' Make size 2 rotation matrix
//'
//' @param angle 
template <class Type>
matrix<Type> makeRotationMatrix(Type angle) {
    matrix<Type> R(2, 2);
    R << cos(angle), -sin(angle), sin(angle), cos(angle);
    return R;
}



//' Make H matrix for Kalman filter
//'
//' @param sigma_obs SD of measurement error
//' @param n_dim Number of dimensions
template<class Type>
matrix<Type> makeH(Type sigma_obs, int n_dim) {
    matrix<Type> H(n_dim, n_dim);
    H.setZero();
    for(int i = 0; i < n_dim; i ++) {
        H(i, i) = sigma_obs * sigma_obs;
    }
    return H;
}

//' Make T matrix for Kalman filter in CTCRW
//' 
//' @param beta Parameter beta of OU velocity process
//' @param dt Length of time interval
//' @param n_dim Number of dimensions of CTCRW process
template<class Type>
matrix<Type> makeT_ctcrw(Type beta, Type dt, int n_dim) {
    matrix<Type> T(2*n_dim, 2*n_dim);
    T.setZero();
    for(int i = 0; i < n_dim; i++) {
        T(2*i, 2*i) = 1;
        T(2*i, 2*i + 1) = (1-exp(-beta*dt))/beta;
        T(2*i + 1, 2*i + 1) = exp(-beta*dt);
    }
    return T;
}

//' Make Q matrix for Kalman filter in CTCRW
//' 
//' @param beta Parameter beta of OU velocity process
//' @param sigma Parameter sigma of OU velocity process
//' @param dt Length of time interval
//' @param n_dim Number of dimensions of CTCRW process
template<class Type>
matrix<Type> makeQ_ctcrw(Type beta, Type sigma, Type dt, int n_dim) {
    matrix<Type> Q(2*n_dim, 2*n_dim);
    Q.setZero();
    for(int i = 0; i < n_dim; i++) {
        Q(2*i, 2*i) = (sigma/beta)*(sigma/beta)*(dt - 2/beta*(1-exp(-beta*dt)) + 
            1/(2*beta)*(1-exp(-2*beta*dt)));
        Q(2*i, 2*i + 1) = sigma*sigma/(2*beta*beta) * (1 - 2*exp(-beta*dt) + exp(-2*beta*dt));
        Q(2*i + 1, 2*i) = Q(2*i, 2*i + 1);
        Q(2*i + 1, 2*i + 1) = sigma*sigma/(2*beta) * (1-exp(-2*beta*dt));
    }
    return Q;
}

//' Make B matrix for Kalman filter in CTCRW
//' 
//' @param beta Parameter beta of OU velocity process
//' @param dt Length of time interval
//' @param n_dim Number of dimensions of CTCRW process
template<class Type>
matrix<Type> makeB_ctcrw(Type beta, Type dt, int n_dim) {
    matrix<Type> B(2*n_dim, n_dim);
    B.setZero();
    for(int i = 0; i < n_dim; i++) {
        B(2*i, i) = dt - (1 - exp(-beta*dt))/beta;
        B(2*i + 1, i) = 1 - exp(-beta*dt);
    }
    return B;
}



//' Make T matrix in transition density in RACVM
//'
//' @param beta Parameter beta of RACVM
//' @param omega Parameter omega of RACVM
//' @param dt Length of time interval
template<class Type>
matrix<Type> makeT_racvm(Type beta, Type omega, Type dt) {

    matrix<Type> A(2,2);
    A << beta,-omega,omega,beta;

    Type C = beta*beta+omega*omega;
    matrix<Type> invA(2,2);
    invA << beta/C,omega/C,-omega/C,beta/C;

    matrix<Type> I(2,2);
    I << 1,0,0,1;

    matrix<Type> R(2,2);
    R << cos(omega*dt),sin(omega*dt),-sin(omega*dt),cos(omega*dt);
    matrix<Type> expAdt(2,2);
    expAdt<<exp(-beta*dt)*R;

    matrix<Type> IntexpAdt(2,2);
    IntexpAdt << invA*(I-expAdt);

    matrix<Type> T(4, 4);
    T.setZero();

    // Combine the matrices into T

    // Top-left block
    T(0, 0) = 1;
    T(1, 1) = 1;

    // Top-right block
    T(0, 2) = IntexpAdt(0, 0);
    T(0, 3) = IntexpAdt(0, 1);
    T(1, 2) = IntexpAdt(1, 0);
    T(1, 3) = IntexpAdt(1, 1);

    // Bottom-right block
    T(2, 2) = expAdt(0, 0);
    T(2, 3) = expAdt(0, 1);
    T(3, 2) = expAdt(1, 0);
    T(3, 3) = expAdt(1, 1);

   
    return T;
}

//' Make Q matrix in transition density in RACVM
//'
//' @param beta Parameter beta of RACVM
//' @param sigma Parameter sigma of RACVM
//' @param omega Parameter omega od RACVM
//' @param dt Length of time interval
template<class Type>
matrix<Type> makeQ_racvm(Type beta, Type sigma,Type omega, Type dt) {


    //constants
    Type tau = 1/beta;
    Type C = beta*beta+omega*omega;

    // variances and covariances values
    Type var_xi = sigma*sigma/C*(dt+(omega*omega-3/(tau*tau))/(2/tau*C)-exp(-2*dt/tau)/(2/tau)+
                        2*exp(-dt/tau)*(1/tau*cos(omega*dt)-omega*sin(omega*dt))/C);
    Type var_zeta = sigma*sigma*tau/2*(1-exp(-2*dt/tau));

    Type cov1 = sigma*sigma/(2*C)*(1+exp(-2*dt/tau)-2*exp(-dt/tau)*cos(omega*dt));

    Type cov2 = sigma*sigma/C*(exp(-dt/tau)*sin(omega*dt)-omega/(2/tau)*(1-exp(-2*dt/tau)));

    //initialize matrix with zeros
    matrix<Type> Q(4, 4);
    Q.setZero();

    // diagonal elements
    Q(0, 0) = Q(1, 1) = var_xi;
    Q(2, 2) = Q(3, 3) = var_zeta;
    
    //non-zero off diagonal elements
    Q(0, 2) = Q(2, 0) = cov1;
    Q(1, 3) = Q(3, 1) = cov1;
    Q(0, 3) = Q(3, 0) = cov2;
    Q(1, 2) = Q(2, 1) = -cov2;

    return Q;
}



//' Make B matrix in transition density in RACVM
//'
//' @param beta Parameter beta of RACVM
//' @param omega Parameter omega of RACVM
//' @param dt Length of time interval
template<class Type>
matrix<Type> makeB_racvm(Type beta,Type omega, Type dt) {
    
    matrix<Type> A(2,2);
    A << beta,-omega,omega,beta;

    Type C = beta*beta+omega*omega;
    matrix<Type> invA(2,2);
    invA <<beta/C,omega/C,-omega/C,beta/C;

    matrix<Type> I(2,2);
    I << 1,0,0,1;

    matrix<Type> R(2,2);
    R << cos(omega*dt),sin(omega*dt),-sin(omega*dt),cos(omega*dt);
    matrix<Type> expAdt(2,2);
    expAdt<<exp(-beta*dt)*R;
    matrix<Type> IntexpAdt(2,2);
    IntexpAdt << invA*(I-expAdt);

    matrix<Type> B(4,2);
    B.setZero();

    // Top block
    B(0, 0) = dt * I(0, 0) - IntexpAdt(0, 0);
    B(0, 1) = dt * I(0, 1) - IntexpAdt(0, 1);
    B(1, 0) = dt * I(1, 0) - IntexpAdt(1, 0);
    B(1, 1) = dt * I(1, 1) - IntexpAdt(1, 1);

    // Bottom block
    B(2, 0) = I(0, 0) - expAdt(0, 0);
    B(2, 1) = I(0, 1) - expAdt(0, 1);
    B(3, 0) = I(1, 0) - expAdt(1, 0);
    B(3, 1) = I(1, 1) - expAdt(1, 1);

    return B;
}


//' Make T matrix in transition density in CRCVM with interpolated distances and angles
//'
//' @param beta Parameter beta of CRCVM
//' @param omega Vector of parameters omega of CRCVM
//' @param dt Vector of length of time intervals
template<class Type>
matrix<Type> makeT_crcvm(Type beta, vector<Type> omega, vector<Type> dt) {

    // // constants
    int m = omega.size();

    matrix<Type> T(4, 4);
    T << 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1;

    for(int i= 0; i < m; i++) {
        
        matrix<Type> Ti = makeT_racvm(beta,omega(i),dt(i));
        T=Ti*T;
    }

return T;
}



//' Make Q matrix in transition density in CRCVM with interpolated distances and angles
//'
//' @param beta Parameter beta of CRCVM
//' @param sigma Parameter sigma of CRCVM
//' @param omega Vector of parameters omega of CRCVM
//' @param dt Vector of length of time intervals
template<class Type>
matrix<Type> makeQ_crcvm(Type beta, Type sigma,vector<Type> omega, vector<Type> dt) {


    int n = omega.size();
    Type tau = 1 / beta;
   

    // Block matrices initialization
    matrix<Type> Q(4, 4);
    Q.setZero();

    // // Loop through time steps
     for(int k = 0; k < n; k++) {

        matrix<Type> Pk(4, 4);
        Pk << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;

        for (int j = k+1; j <n; j++) {

            matrix<Type> Tj=makeT_racvm(beta,omega(j),dt(j));
            Pk =  Pk * Tj;
        }


        matrix<Type> Qk= makeQ_racvm(beta,sigma,omega(k),dt(k));

        Q=Q+Pk*Qk*Pk.transpose();
    }
        
    return Q;

}

#endif
