#include <armadillo>
#include "ellincheck.h"
#include <cstdlib>
#include <chrono>

using namespace arma;

static unsigned int _i = 0;

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
    this->_max = 1;
    this->_beta_ast = 1;

    this->betaMax = 1-arma::dot(c, c); 
    this->betaMin = 1/arma::min(lb); 
    
}
void l_cp::update(const double beta){        
    this->_l = 1 - beta;
    this->_dl  = -1.0;
    this->_ddl =  0.0;
    for (size_t i = 0; i < this->_n; i++)
    {
        const double lbi_beta = this->_lb[i]*beta;
        const double one_mlbibt = (1 - lbi_beta);
        const double cisqr_over_one_mlbibt = this->_csqr[i]/one_mlbibt;
        const double one_mlbibt2 = one_mlbibt*one_mlbibt;
        // const double one_mlbibt3 = one_mlbibt*one_mlbibt2;  
        // cout << "x: " << lbi_beta    << endl;
        // cout << "y: " << one_mlbibt2 << endl;
        // cout << "z: " << one_mlbibt3 << endl;
        this->_l += cisqr_over_one_mlbibt * lbi_beta; 
        this->_dl += cisqr_over_one_mlbibt * this->_lb[i] / one_mlbibt; 
        this->_ddl += 2 * cisqr_over_one_mlbibt * this->_lbsqr[i] / one_mlbibt2; 
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

double l_cp::max(){
    double beta;
    unsigned int i;
    unsigned int max_iter = MAX_ITER;
    return l_cp::max(beta, i, max_iter);

}
double l_cp::max(double& beta){
    unsigned int i;
    unsigned int max_iter = MAX_ITER;
    return l_cp::max(beta, i, max_iter);

}
double l_cp::max(double& beta, unsigned int& i){
    unsigned int max_iter = MAX_ITER;
    return l_cp::max(beta, i, max_iter);

}

double l_cp::max(double& beta, unsigned int &i = _i, unsigned int max_iter = MAX_ITER){
    i = 0;
    beta = this->_beta_ast;
    if (this->_max < 1){
        cout << "Warning: cached" << endl; // TODO no cache
        return this->_max;
    }
    // cout << "beta in [" << this->betaMin << ", " << this->betaMax << "]" << endl;
    if (this->betaMax<this->betaMin)
    {
        cout << "Warning: empty interval" << endl;
        return this->_max = -INFINITY; // empty interval
    }
    
    this->update(beta);
    if(this->df()<=0)
    {
        for ( ; i < max_iter; i++)
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
    this->_beta_ast = beta;
    return this->_max = this->f();
}

ellipsoid::ellipsoid(mat P, vec c){
    this->P = P;
    this->c = c;
    this->_n = c.size();
}
double ellipsoid::volume(){ // TODO
    cout << "TODO: volume not yet implemented" << endl;
    if (this->_volume == 0){
        this->_volume = 0;
        return 0; 
    }
    return this->_volume;
}


bool ellipsoid::in(ellipsoid el0) {
    
using namespace std::chrono;    
    auto t_start = high_resolution_clock::now();
    // INIT
    vec c;
    vec lb(this->_n);
    double max;
    
    auto t_preconfig = high_resolution_clock::now();
    // CONFIG
    const mat L0 = arma::chol(el0.P);
    const mat L0_inv = arma::inv(L0);

    mat Pt = L0 * this->P * L0.t();
    const mat ct = L0.t() * (this->c - el0.c);
    
    mat V(this->_n, this->_n);

    eig_sym(lb, V, Pt);
    c = V.t() * ct;
    auto t_preldef = high_resolution_clock::now();
    // LDEF
    l_cp l(c,lb);
    auto t_preopt = high_resolution_clock::now();
    // OPT
    max = l.max()>=0;
    auto t_final = high_resolution_clock::now();

    auto duration_init = duration_cast<microseconds>(t_preconfig - t_start);
    auto duration_config = duration_cast<microseconds>(t_preldef - t_preconfig);
    auto duration_ldef = duration_cast<microseconds>(t_preopt - t_preldef);
    auto duration_opt = duration_cast<microseconds>(t_final - t_preopt);
    cout << "===Profilling:" << endl;
    cout << "duration_init:\t" << duration_init.count() << endl;
    cout << "duration_config:\t" << duration_config.count() << endl;
    cout << "duration_ldef:\t" << duration_ldef.count() << endl;
    cout << "duration_opt:\t" << duration_opt.count() << endl;
    return max; 
}


