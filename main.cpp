#include <iostream>
#include <armadillo>
#include <chrono>

using namespace std;
using namespace arma;
using namespace std::chrono;

class l_cp
{
private:
    int n;
    vec csqr;
    vec lb;

public:
    l_cp(vec c,  vec lb)
    {
        this->csqr = arma::square(c);
        this->lb = lb;
        this->n = lb.size();
    }

    double f(const double beta){
        const vec lb_beta = this->lb * beta;
        double out = 1 - beta;
        for (size_t i = 0; i < this->n; i++)
        {
            out += this->csqr[i] *lb_beta[i]/(1-lb_beta[i]); 
        }
        return out;
    }
    double f2(const double beta){
        const vec aux = this->csqr % (this->lb * beta) / (1-(this->lb * beta));
        return (1 - beta) + sum(aux);
    }

    double df(const double beta){
        double out = 1 - beta;
        return out;
    }

    double ddf(const double beta){
        double out = 1;
        return out;
    }
};





int main()
   {
    const int n = 600;
    vec c(n, fill::randu);
    vec lb(n, fill::randu);
    //const vec* csqr;
    //const vec* lb;

    c = c *0.001; // ensure inside
    lb = arma::square(lb); //ensure positive

    // algorithm starts here
    l_cp l(c,lb);
    double beta = 1;
    for (size_t i = 0; i < 50; i++)
    {
        if (beta>1)
            break;
        beta += l.df(beta) * 0.01;
    }
    double res;
    const int N = 10000;
    // f1 timming
    {
        auto start = high_resolution_clock::now();
        for (size_t i = 0; i < N; i++)
        {
            res = l.f2(0.2);
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << duration.count() << endl;
    }
    // f2 timming
    {
        auto start = high_resolution_clock::now();
        for (size_t i = 0; i < N; i++)
        {
            res = l.f(0.2);
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << duration.count() << endl;
    }
    return 0; 
}