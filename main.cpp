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
    
    lb = 1 / arma::abs(lb); //ensure positive
    if (CONTAINED)
        c = c * 0.25*(1-1/sqrt(min(lb))) / norm(c); // ensure inside
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
    const unsigned int n = 100;
    
    // test lb;
    vec c(n, arma::fill::randu);
    mat P(n, n, arma::fill::randu);
    vec c0(n, arma::fill::randu);
    mat P0(n, n, arma::fill::randu);
    
    P *= P.t();
    P0 *= P0.t()*0.0001;
    c = c0 + c*0.0001;

    ellipsoid el(P, c);
    ellipsoid el0(P0, c0);
    bool resp;
    auto start = high_resolution_clock::now();
    resp = el.in(el0);
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    cout <<  "Test el in el0" << endl;
    cout << "t:\t" << duration.count() << endl;
    cout << "In:\t" << resp << endl;
    
}


int main()
   {
    test_lb();
    test_ells(); 
    return 0; 
}