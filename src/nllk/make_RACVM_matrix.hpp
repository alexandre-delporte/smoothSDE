#ifndef _make_RACVM_matrix_
#define _make_RACVM_matrix_

using namespace R_inla;
using namespace density;
using namespace Eigen;


//' Make T matrix in transition density
//'
//' @param beta Parameter beta of RACVM
//' @param omega Parameter omega of RACVM
//' @param dt Length of time interval
template<class Type>
matrix<Type> makeT_racvm(Type beta, Type omega, Type delta) {

    matrix<Type> A(2,2);
    A << beta,-omega,omega,beta;

    Type C = beta*beta+omega*omega;
    matrix<Type> invA(2,2);
    invA << beta/C,omega/C,-omega/C,beta/C;

    matrix<Type> I(2,2);
    I << 1,0,0,1;

    matrix<Type> R(2,2);
    R << cos(omega*delta),sin(omega*delta),-sin(omega*delta),cos(omega*delta);
    matrix<Type> expAdelta(2,2);
    expAdelta<<exp(-beta*delta)*R;
    matrix<Type> IntexpAdelta(2,2);
    IntexpAdelta << invA*(I-expAdelta);

    matrix<Type> T(4, 4);
    T.setZero();

    // Combine the matrices into T

    // Top-left block
    T(0, 0) = 1;
    T(1, 1) = 1;

    // Top-right block
    T(0, 2) = IntexpAdelta(0, 0);
    T(0, 3) = IntexpAdelta(0, 1);
    T(1, 2) = IntexpAdelta(1, 0);
    T(1, 3) = IntexpAdelta(1, 1);

    // Bottom-right block
    T(2, 2) = expAdelta(0, 0);
    T(2, 3) = expAdelta(0, 1);
    T(3, 2) = expAdelta(1, 0);
    T(3, 3) = expAdelta(1, 1);

   
    return T;
}

//' Make Q matrix in transition density
//'
//' @param beta Parameter beta of RACVM
//' @param sigma Parameter sigma of RACVM
//' @param omega Parameter omega od RACVM
//' @param delta Length of time interval
template<class Type>
matrix<Type> makeQ_racvm(Type beta, Type sigma,Type omega, Type delta) {


    //constants
    Type tau = 1/beta;
    Type C = beta*beta+omega*omega;

    // variances and covariances values
    Type var_xi = sigma*sigma/C*(delta+(omega*omega-3/(tau*tau))/(2/tau*C)-exp(-2*delta/tau)/(2/tau)+
                        2*exp(-delta/tau)*(1/tau*cos(omega*delta)-omega*sin(omega*delta))/C);
    Type var_zeta = sigma*sigma*tau/2*(1-exp(-2*delta/tau));

    Type cov1 = sigma*sigma/(2*C)*(1+exp(-2*delta/tau)-2*exp(-delta/tau)*cos(omega*delta));

    Type cov2 = sigma*sigma/C*(exp(-delta/tau)*sin(omega*delta)-omega/(2/tau)*(1-exp(-2*delta/tau)));

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



//' Make B matrix in transition density
//'
//' @param beta Parameter beta of RACVM
//' @param omega Parameter omega of RACVM
//' @param dt Length of time interval
template<class Type>
matrix<Type> makeB_racvm(Type beta,Type omega, Type delta) {
    
    matrix<Type> A(2,2);
    A << beta,-omega,omega,beta;

    Type C = beta*beta+omega*omega;
    matrix<Type> invA(2,2);
    invA <<beta/C,omega/C,-omega/C,beta/C;

    matrix<Type> I(2,2);
    I << 1,0,0,1;

    matrix<Type> R(2,2);
    R << cos(omega*delta),sin(omega*delta),-sin(omega*delta),cos(omega*delta);
    matrix<Type> expAdelta(2,2);
    expAdelta<<exp(-beta*delta)*R;
    matrix<Type> IntexpAdelta(2,2);
    IntexpAdelta << invA*(I-expAdelta);

    matrix<Type> B(4,2);
    B.setZero();

    // Top block
    B(0, 0) = delta * I(0, 0) - IntexpAdelta(0, 0);
    B(0, 1) = delta * I(0, 1) - IntexpAdelta(0, 1);
    B(1, 0) = delta * I(1, 0) - IntexpAdelta(1, 0);
    B(1, 1) = delta * I(1, 1) - IntexpAdelta(1, 1);

    // Bottom block
    B(2, 0) = I(0, 0) - expAdelta(0, 0);
    B(2, 1) = I(0, 1) - expAdelta(0, 1);
    B(3, 0) = I(1, 0) - expAdelta(1, 0);
    B(3, 1) = I(1, 1) - expAdelta(1, 1);

    return B;
}

#endif
