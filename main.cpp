#include <iostream>
#include <armadillo>
#include <chrono>
#include "ellincheck/ellincheck.h"

using namespace std;
using namespace arma;
using namespace std::chrono;



bool test_lb(bool contained){
    arma_rng::set_seed_random();
    const unsigned int n = 3;
    unsigned int i = 0;
    
    // test lb;
    vec c(n, arma::fill::randu);
    vec lb(n, arma::fill::randu);
    
    lb = 100 / arma::abs(lb); //ensure positive
    if (contained)
        c = c * 0.9*(1-1/sqrt(min(lb))) / norm(c); // ensure inside
    else 
        c = c * (1-0.5/sqrt(max(lb))) / norm(c); // ensure outside
    
    // define function l_cp
    l_cp l(c,lb);

    // iterate
    double beta = l.betaMax;

    auto start = high_resolution_clock::now();
    l.max(beta, i);
    auto stop = high_resolution_clock::now();


    auto duration = duration_cast<microseconds>(stop - start);
    cout << duration.count() << endl;
    cout << "iter:\t" << i << endl;
    cout << "beta:\t" << beta << endl;
    cout << "gk:\t" << l.newtownIterate() << endl;
    cout << "l:\t" << l.f() << endl;
    cout << "dl:\t" << l.df() << endl;
    cout << "ddl:\t" << l.ddf() << endl << endl;
    return contained == (l.f()>=0);

}



bool test_ells(bool contained){
    arma_rng::set_seed_random();
    const unsigned int n = 2;
    
    // test lb;
    // vec c(n, arma::fill::randu);
    // mat P(n, n, arma::fill::randu);
    // vec c0(n, arma::fill::randu);
    // mat P0(n, n, arma::fill::randu);
    

    // P *= P.t()*100;
    // P0 *= P0.t()*0.0001;
    // c = c0 + c*0.0001;
    if (contained){
        double c[n] = {1.5, 1.5};
        double P[n][n] = {{4.0, 0.5},       
                    {0.5, 6.0}};
        double c0[n] = {1.6, 1.4};
        double P0[n][n] = {{0.4, -0.1},
                    {-0.1, 0.5}};
        return contained == ellincheck(c, *P, c0, *P0, n);
    }
    else{
        double c[n] = {1.6, 1.4};
        double P[n][n] = {{0.4, -0.1},
                    {-0.1, 0.5}};
        double c0[n] = {1.5, 1.5};
        double P0[n][n] = {{4.0, 0.5},       
                    {0.5, 6.0}};
        return contained == ellincheck(c, *P, c0, *P0, n);
    }
            



}


int main()
   {
    bool allOk = true;
    allOk &= test_lb(CONTAINED);

    allOk &= test_lb(NOT_CONTAINED);
    allOk &= test_ells(CONTAINED); 
    allOk &= test_ells(NOT_CONTAINED); 
    cout << "Test res: " << allOk << endl;

    return !allOk; 
}