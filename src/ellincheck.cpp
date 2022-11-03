#include <armadillo>
#include "ellincheck.h"
#include <cstdlib>
#include <chrono>

using namespace arma;

static unsigned int _i = 0;

l_cp::l_cp(){

}

l_cp::l_cp(vec c, vec lb)
{
    using namespace arma;
    this->_c = c;
    this->_csqr = square(c);
    this->_lb = lb;
    this->_lbsqr = lb % lb;
    this->_n = lb.size();
    this->_l = 0;
    this->_dl = 0;
    this->_ddl = 0;
    this->_max = 1;
    this->_beta_ast = 1;
    this->_i_ast = 1;
    
    unsigned int idx = index_min(lb);
    const double lmin = lb[idx];
    const double cmin = c[idx];
    const double b = cmin*cmin*lmin-1-lmin;

    this-> _delta_ub = b*b-4*lmin;

    if (_delta_ub>0)
    {
        this->betaMax = std::min(1-dot(c, c), (-b-sqrt(_delta_ub))/2/lmin);
        this->betaMin = std::max(1/lmin, (-b+sqrt(_delta_ub))/2/lmin);
    }
    else
    {
        this->betaMax = 1-dot(c, c); 
        this->betaMin = 1-dot(c, c);//1/lmin; 
    }

    // cout << "c:\t" << c << endl;
    // cout << "lb:\t" << lb << endl;

    
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

double l_cp::max(double& beta, unsigned int &i, unsigned int max_iter){
    i = 0;
    beta = this->_beta_ast;
    if (this->_max < 1){
        cout << "Warning: cached" << endl; // TODO no cache
        i = this->_i_ast;
        return this->_max;
    }
    cout << "beta in [" << this->betaMin << ", " << this->betaMax << "]" << endl;
    if (this->betaMax<this->betaMin)
    {
        cout << "Warning: empty interval" << endl;
        return this->_max = -INFINITY; // empty interval
    }
    beta = this->betaMax;
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
    this->_i_ast = i;
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

bool ellipsoid::included_in(ellipsoid el0){
    l_cp l;
    return ellipsoid::included_in(el0, l);
}

bool ellipsoid::included_in(ellipsoid el0, l_cp& l) {
        
    using namespace std::chrono;   
    using namespace arma; 
    auto t_start = high_resolution_clock::now();
    // INIT
    vec c;
    vec lb(this->_n);
    mat V(this->_n, this->_n);
    double max;
    
    auto t_preconfig = high_resolution_clock::now();
    // CONFIG
    mat L0 = chol(el0.P);
    mat L0_inv = inv(L0);

    mat Pt = L0_inv * this->P * L0_inv.t();
    // mat Pt = solve(L0, solve(L0, this->P).t()).t();
    mat ct = L0.t() * (this->c - el0.c);
    
    // cout << "Pt:\t" << Pt << endl;
    // cout << "ct:\t" << ct << endl;
    auto t_preldef = high_resolution_clock::now();
    // LDEF

    eig_sym(lb, V, Pt);
    c = V.t() * ct;

    l = l_cp(c,lb);
    auto t_preopt = high_resolution_clock::now();
    // OPT
    max = l.max();
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
    return max>=0; 
}

bool ellincheck(double *_c, double *_P, double *_c0, double *_P0, unsigned int n){
    using namespace std::chrono;
    vec c(_c, n);
    mat P(_P, n, n);
    vec c0(_c0, n);
    mat P0(_P0, n, n);
    cout << "c:\t" << c << endl;
    cout << "P:\t" << P << endl;
    cout << "c0:\t" << c0 << endl;
    cout << "P0:\t" << P0 << endl;

    ellipsoid el(P, c);
    ellipsoid el0(P0, c0);
    bool resp;
    l_cp l;
    auto start = high_resolution_clock::now();
    resp = el.included_in(el0, l);
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    cout <<  "Test el in el0" << endl;
    cout << "t:\t" << duration.count() << endl;
    cout << "In:\t" << resp << endl;
    double beta;
    cout << "l:\t" << l.max(beta) << endl;
    cout << "beta:\t" << beta << endl;
    return resp;
}

