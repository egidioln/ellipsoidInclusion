#ifndef ELLINCHECK_LIB_HPP
#define ELLINCHECK_LIB_HPP
#include <armadillo>



#define MAX_ITER 10000
#define CONTAINED 1
#define EPS_STOP 1E-6
#define EPS_CONST 1E-6

using namespace arma;

class l_cp
{
private:
    unsigned int _n;
    vec _csqr;
    vec _lb;
    vec _lbsqr;
    vec _c;
    double _l;
    double _dl;
    double _ddl;
    

public:
    double betaMin;
    double betaMax;
    l_cp(vec c, vec lb);
    void update(const double beta);
    double f();
    double df();
    double ddf();
    double newtownIterate();
    unsigned int max(double& beta);
};



#endif // ELLINCHECK_LIB_HPP
