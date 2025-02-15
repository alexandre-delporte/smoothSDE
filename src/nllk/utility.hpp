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


//' Slice a vector
//'
//' @param v vector
//' @param start_idx integer
//' @param end_idx integer
template<class Type>
vector<Type> slice(vector<Type> v, int start_idx, int end_idx) {
    

    // Return empty vector if indices are invalid
    if (start_idx > end_idx) {
        return vector<Type>(0);
    }

    vector<Type> subvector(end_idx - start_idx + 1);
    
    for (int i = start_idx; i <= end_idx; ++i) {
        subvector(i - start_idx) = v(i); 
    }

    return subvector;
}

//' Cumulative sum of a vector
//'
//' @param v vector
template <class Type>
vector<Type> cumsum(vector<Type> v) {
    int n = v.size();
    vector<Type> result(n);
    
    if (n == 0) return result;
    
    result[0] = v[0];
    
    for (int i = 1; i < n; ++i) {
        result[i] = result[i - 1] + v[i];
    }

    return result;
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
    matrix<Type> I(2, 2);
    I << 1, 0, 0, 1;

    // Block matrices initialization
    matrix<Type> top_left_block(2, 2), top_right_block(2, 2);
    top_left_block.setZero(); top_right_block.setZero();

    matrix<Type> bottom_right_block = sigma * sigma * tau / 2 * (1 - exp(-2 * dt.sum() / tau))*I;

     // Make last rotation matrix
     matrix<Type> R_n=makeRotationMatrix(-omega(n-1)*dt(n-1));

    // // Loop through time steps
     for(int j = 1; j < n+1; j++) {

        Type C_j = beta * beta + omega(j-1) * omega(j-1);

        // Compute q1_j and q2_j terms
        Type q1_j = sigma * sigma / C_j * (
            dt(j-1) + (omega(j-1) * omega(j-1) - 3 / (tau * tau)) / (2 / tau * C_j) -
            exp(-2 * dt(j-1) / tau) / (2 / tau) +
            2 * exp(-dt(j-1) / tau) * (1 / tau * cos(omega(j-1)*dt(j-1)) - omega(j-1) * sin(omega(j-1)*dt(j-1))) / C_j
        );

        Type q2_j = sigma * sigma * tau / 2 * (1 - exp(-2 * dt(j-1) / tau));

        // Compute covariance matrix Gamma_j
        matrix<Type> Gamma_j(2, 2);
        Gamma_j.setZero();
        Gamma_j(0, 0) = Gamma_j(1, 1) = sigma * sigma / (2 * C_j) * (1 + exp(-2 * dt(j-1) / tau) - 2 * exp(-dt(j-1) / tau) * cos(omega(j-1)*dt(j-1)));
        Gamma_j(0, 1) = sigma * sigma / C_j * (exp(-dt(j-1) / tau) * sin(omega(j-1)*dt(j-1)) - omega(j-1)/ (2 / tau) * (1 - exp(-2 * dt(j-1)/ tau)));
        Gamma_j(1, 0) = -Gamma_j(0, 1);
        
        // Compute M_j
        matrix<Type> M_j(2,2);
        M_j.setZero();

        for(int k = j; k < n; k++) {

            // Compute Pjk
            matrix<Type> P_jk = I;

            for(int i = j; i < k; i++) {

                matrix <Type> R =makeRotationMatrix(-omega(i)*dt(i));

                P_jk=P_jk * exp(-beta*dt(i))*R;

            }

            // Compute S_k
            matrix<Type> A_k(2,2);
            A_k << beta,-omega(k),omega(k),beta;
            Type C_k = beta*beta+omega(k)*omega(k);
            matrix<Type> invA_k(2,2);
            invA_k << beta/C_k,omega(k)/C_k,-omega(k)/C_k,beta/C_k;

            matrix<Type> R_k = makeRotationMatrix(-omega(k)*dt(k));
            matrix<Type> S_k = invA_k*(I-exp(-beta*dt(k))*R_k);

            M_j = M_j + P_jk * S_k;
        }

        top_left_block = top_left_block + q1_j*I+M_j.transpose()*Gamma_j+Gamma_j.transpose()*M_j+
                        M_j*M_j.transpose()*q2_j*I;

       
        
        //Compute P_nj and factor in velocity variance
        matrix<Type> P_jn = I;

        for(int i = j; i < n; i++) {

            matrix <Type> R =makeRotationMatrix(-omega(i)*dt(i));

            P_jn=P_jn * exp(-beta*dt(i))*R;
        }

        top_right_block = top_right_block+P_jn.transpose()*Gamma_j
                        + M_j*P_jn.transpose()*q2_j*I;
        
     }


    matrix<Type> bottom_left_block=top_right_block.transpose();

    //initialize matrix with zeros
    matrix<Type> Q(4, 4);
    Q.setZero();

    // diagonal elements
    Q(0, 0) = Q(1, 1) = top_left_block(0,0);
    Q(2, 2) = Q(3, 3) = bottom_right_block(0,0);
    
    //non-zero off diagonal elements
    Q(0, 2) = Q(2, 0) = top_right_block(0,0);
    Q(1, 3) = Q(3, 1) = top_right_block(1,1);
    Q(0, 3) = Q(3, 0) = top_right_block(0,1);
    Q(1, 2) = Q(2, 1) = top_right_block(1,0);


    return Q;

}

#endif
