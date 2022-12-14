# ellipsoidInclusion

| **Build Status** |
|:----------------:|
| [![Build Status][build-img]][build-url] |
| [![codecov][codecov-img]][codecov-url] |
<!-- |  [![Codecov branch][codecov-img]][codecov-url] | -->
[build-img]: https://github.com/egidioln/ellipsoidInclusion/workflows/CMake/badge.svg?branch=main
[build-url]: https://github.com/egidioln/ellipsoidInclusion/actions?query=workflow%3ACMake
[codecov-img]: https://codecov.io/gh/egidioln/ellipsoidInclusion/branch/main/graph/badge.svg?token=8DUhQe22qD
[codecov-url]: https://codecov.io/gh/egidioln/ellipsoidInclusion
This is a c++  implementation of a function that checks the inclusion of one n-ellipsoid in another. For a positive definite matrix $P\succ0\in\mathbb{R}^{n\times n}$ and a vector $c\in\mathbb{R}^{n}$, an *ellipsoid shaped by* $P$ *and cetered at* $c$ is defined as $E(P,c) := \\{x\in\mathbb{R}^{n}:(x-c)^\top P(x-c)\leq 1\\}$.

The following code fragment tests the inclusion $E(P,c) \subseteq E(P_0,c_0) $

```cpp
    double c[n] = {1.5, 1.5};
    double P[n][n] = {{4.0, 0.5},       
                 {0.5, 6.0}};
    double c0[n] = {1.6, 1.4};
    double P0[n][n] = {{0.4, -0.1},
                   {-0.1, 0.5}};

    if ellincheck(c, *P, c0, *P0, n){
     // do stuff
    }
```

Classes can also be defined as follows:

```cpp
  ellipsoid el(P, c);
  ellipsoid el0(P0, c0);

  el.included_in(el0)

```
where `P` and `P0` must be `arma::mat` and `c` and `c0` must be `arma::vec` ([check the Armadillo c++ library](https://arma.sourceforge.net/docs.html) for more details on these types).

The method implemented in this library is available in [this paper](https://arxiv.org/abs/2211.06237). Please, cite it as:
```
@misc{calbert2022efficient,
  doi = {10.48550/ARXIV.2211.06237},
  url = {https://arxiv.org/abs/2211.06237},
  author = {Calbert, Julien and Egidio, Lucas N. and Jungers, Raphaël M.},
  title = {An Efficient Method to Verify the Inclusion of Ellipsoids},
  publisher = {preprint (arXiv)},
  year = {2022},
}

```

