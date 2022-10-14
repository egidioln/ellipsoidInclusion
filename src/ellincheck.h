#ifndef ELLINCHECK_LIB_HPP
#define ELLINCHECK_LIB_HPP
#include <armadillo>



#define MAX_ITER 10000
#define EPS_STOP 1E-6
#define EPS_CONST 1E-6
#define CONTAINED true
#define NOT_CONTAINED false

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
    double _max;
    double _beta_ast;
    unsigned int _i_ast;
    

public:
    double betaMin;
    double betaMax;
    l_cp();
    l_cp(vec c, vec lb);
    void update(const double);
    double f();
    double df();
    double ddf();
    double newtownIterate();
    double max();
    double max(double&);
    double max(double&, unsigned int &);
    double max(double&, unsigned int &, unsigned int);
};


class ellipsoid
{
private:
    unsigned int _n;
    double _volume;

public:
    mat P;
    vec c;
    ellipsoid(mat , vec);
    double volume();

    bool included_in(ellipsoid);
    bool included_in(ellipsoid, l_cp&);
};

bool ellincheck(double*, double*, double*, double*, unsigned int);

#endif // ELLINCHECK_LIB_HPP
