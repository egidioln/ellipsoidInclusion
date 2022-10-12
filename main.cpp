#include <iostream>
#include <armadillo>
#include <chrono>
#define MAX_ITER 10000
#define CONTAINED 0
#define EPS_STOP 1E-6
#define EPS_CONST 1E-6

using namespace std;
using namespace arma;
using namespace std::chrono;

class l_cp
{
private:
    unsigned int _n;
    vec _csqr;
    vec _lb;
    vec _lbsqr;
    double _l;
    double _dl;
    double _ddl;
    

public:
    double lbmin;
    l_cp(vec c, vec lb)
    {
        this->_csqr = arma::square(c);
        this->_lb = lb;
        this->_lbsqr = lb % lb;
        this->_n = lb.size();
        this->_l = 0;
        this->_dl = 0;
        this->_ddl = 0;
        this->lbmin = min(lb);
    }
    void update(const double beta){        
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
    double f(){
        return this->_l;
    }
    double df(){
        return this->_dl;
    }
    double ddf(){
        return this->_ddl;
    }
    

    double newtownIterate(){
        return 1*this->_dl/this->_ddl;
    }
};





int main()
   {
    arma_rng::set_seed_random();
    const unsigned int n = 200;
    unsigned int i = 0;
    vec c(n, arma::fill::randu);
    vec lb(n, arma::fill::randu);

    lb = 1 / arma::abs(lb); //ensure positive
    if (CONTAINED)
        c = c * 0.25*(1-1/sqrt(min(lb))) / norm(c); // ensure inside
    else 
        c = c * 0.51*(1-1/sqrt(max(lb))) / norm(c); // ensure inside
    
    // define function l_cp
    l_cp l(c,lb);

    const double betaMax = 1-arma::dot(c,c); 
    const double betaMin = 1/l.lbmin; 
    cout << "beta in [" << betaMin << ", " << betaMax << "]" << endl;
    if (betaMax<betaMin)
    {
        cout << ((CONTAINED == 0)? "CORRECT" : "**INCORRECT**") << endl;
        return 0;
    }
    // iterate
    double beta = betaMax;

    l.update(beta);
    auto start = high_resolution_clock::now();
    if(l.df()<=0)
        for ( ; i < MAX_ITER; i++)
        {
            if (norm(l.df())<EPS_STOP)
                break;
            beta -= l.newtownIterate();
            if (beta>betaMax)
                beta = betaMax; 
            else if (beta<betaMin)
                beta = betaMin+EPS_CONST;
            l.update(beta);
        }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << duration.count() << endl;
    cout << "iter:\t" << i << endl;
    cout << "beta:\t" << beta << endl;
    cout << "gk:\t" << l.newtownIterate() << endl;
    cout << "l:\t" << l.f() << endl;
    cout << "dl:\t" << l.df() << endl;
    cout << "ddl:\t" << l.ddf() << endl << endl;
    cout << ((CONTAINED == (l.f()>=0))? "CORRECT" : "**INCORRECT**") << endl;
    // double res;
    // const int N = 10000;
    // // f1 timming
    // {
    //     auto start = high_resolution_clock::now();
    //     for (size_t i = 0; i < N; i++)
    //     {
    //         res = l.f2(0.2);
    //     }
    //     auto stop = high_resolution_clock::now();
    //     auto duration = duration_cast<microseconds>(stop - start);
    //     cout << duration.count() << endl;
    // }
    // // f2 timming
    // {
    //     auto start = high_resolution_clock::now();
    //     for (size_t i = 0; i < N; i++)
    //     {
    //         res = l.f(0.2);
    //     }
    //     auto stop = high_resolution_clock::now();
    //     auto duration = duration_cast<microseconds>(stop - start);
    //     cout << duration.count() << endl;
    // }

    
    return 0; 
}