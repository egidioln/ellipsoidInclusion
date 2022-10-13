#include <armadillo>
#include "ellincheck.h"
#include <cstdlib>
using namespace arma;


l_cp::l_cp(vec c, vec lb)
{
    this->_c = c;
    this->_csqr = arma::square(c);
    this->_lb = lb;
    this->_lbsqr = lb % lb;
    this->_n = lb.size();
    this->_l = 0;
    this->_dl = 0;
    this->_ddl = 0;

    this->betaMax = 1-arma::dot(c,c); 
    this->betaMin = 1/arma::min(lb); 
    
}
void l_cp::update(const double beta){        
    this->_l = 1 - beta;
    this->_dl  = -1.0;
    this->_ddl =  0.0;
    for (size_t i = 0; i < this->_n; i++)
    {
        const double lbi_beta = this->_lb[i]*beta;
        const double one_mlbibt2 = (1-lbi_beta)*(1-lbi_beta);
        const double one_mlbibt3 = (1-lbi_beta)*one_mlbibt2;  
        // cout << "x: " << lbi_beta    << endl;
        // cout << "y: " << one_mlbibt2 << endl;
        // cout << "z: " << one_mlbibt3 << endl;
        this->_l += this->_csqr[i] * lbi_beta / (1 - lbi_beta); 
        this->_dl += this->_csqr[i] * this->_lb[i] / one_mlbibt2; 
        this->_ddl += 2 * this->_csqr[i] * this->_lbsqr[i] / one_mlbibt3; 
    }     
}
double l_cp::f(){
    return this->_l;
}
double l_cp::df(){
    return this->_dl;
}
double l_cp::ddf(){
    return this->_ddl;
}


double l_cp::newtownIterate(){
    return 1*this->_dl/this->_ddl;
}

unsigned int l_cp::max(double& beta){
    unsigned int i = 0;
    cout << "beta in [" << this->betaMin << ", " << this->betaMax << "]" << endl;
    if (this->betaMax<this->betaMin)
    {
        cout << ((CONTAINED == 0)? "CORRECT" : "**INCORRECT**") << endl;
        return 0;
    }
    
    this->update(beta);
    if(this->df()<=0)
    {
        for ( ; i < MAX_ITER; i++)
        {
            if (std::abs(this->df())<EPS_STOP)
                break;
            beta -= this->newtownIterate();
            if (beta>this->betaMax)
                beta = this->betaMax; 
            else if (beta<this->betaMin)
                beta = this->betaMin+EPS_CONST;
            this->update(beta);
        }
    }
    else
        return -1;
    return i;
}

class ellipsoid
{
private:
    /* data */
public:
    ellipsoid(float** P, float* c);
    ~ellipsoid();
};

// ellipsoid::ellipsoid(/* args */)
// {
// }

ellipsoid::~ellipsoid()
{
}
