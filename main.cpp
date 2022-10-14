#include <iostream>
#include <armadillo>
#include <chrono>
#include "ellincheck/ellincheck.h"

using namespace std;
using namespace arma;
using namespace std::chrono;



void test_lb(){
    arma_rng::set_seed_random();
    const unsigned int n = 1000;
    unsigned int i = 0;
    
    // test lb;
    vec c(n, arma::fill::randu);
    vec lb(n, arma::fill::randu);
    
    lb = 100 / arma::abs(lb); //ensure positive
    if (CONTAINED)
        c = c * 0.9*(1-1/sqrt(min(lb))) / norm(c); // ensure inside
    else 
        c = c * 0.51*(1-1/sqrt(max(lb))) / norm(c); // ensure inside
    
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
    cout << ((CONTAINED == (l.f()>=0))? "CORRECT" : "**INCORRECT**") << endl;

}



void test_ells(){
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

    double c[n] = {1.5, 1.5};
    double P[n][n] = {{4.0, 0.5},       
                 {0.5, 6.0}};
    double c0[n] = {1.6, 1.4};
    double P0[n][n] = {{0.4, -0.1},
                   {-0.1, 0.5}};

            



    ellincheck(c, *P, c0, *P0, n);
}


int main()
   {
    test_lb();
    test_ells(); 
    return 0; 
}