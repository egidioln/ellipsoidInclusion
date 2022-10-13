#include <iostream>
#include <armadillo>
#include <chrono>
#include "ellincheck.hpp"

using namespace std;
using namespace arma;
using namespace std::chrono;






int main()
   {
    arma_rng::set_seed_random();
    const unsigned int n = 10000;
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


    // iterate
    double beta = l.betaMax;

    auto start = high_resolution_clock::now();
    i = l.max(beta);
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