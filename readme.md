## Madsen-Reid Polynomial Root Finder in C++

This is a black-box polynomial root finder written in modern C++ (header only, no dependencies except the standard library).

It is based on the 1975 paper: K. Madsen and J. Reid, "Fortran Subroutines for Finding Polynomial Zeros."

To find and print the real and complex roots of a polynomial with real coefficients, such as 5x^3 + 4x^2 + 3x + 2:

```
std::vector<double> polynomial{ 5, 4, 3, 2 };

bitlush::madsen_root_finder<double> madsen;

madsen.find_roots(polynomial.data(), polynomial.size() - 1);

for (std::size_t i = 0; i < polynomial.size() - 1; i++)
{
	std::cout << madsen[i] << std::endl;
}
```

Similarly, polynomials with complex coefficients can be handled using the class `bitlush::madsen_root_finder<std::complex<double>>`.

### Why another root finder?
The goal of this code is to present a clean and clear implementation of the algorithm described in the Madsen and Reid paper. It is not a direct translation of the original Fortran (which is difficult to read). The majority of root-finding implementations are not easy to understand and are hard to build upon.

Examples of the Madsen-Reid algorithm are scarce, yet it performs well and is a good choice despite being over 50 years old. Hopefully, this implementation will help increase its popularity.

### Performance and Accuracy
The Madsen-Reid algorithm performs favorably compared to alternative polynomial root-finding algorithms in terms of speed and robustness, while having the advantage of being easy to understand.

### Future Improvements
One improvement would be to add forward and backward polynomial deflation, as the coefficients of the remainder polynomial after division for roots with magnitudes greater than one increase unbounded with the polynomial degree.

Another improvement would be to integrate tighter stopping criteria based on research post-dating the Madsen-Reid paper.

### The State of the Art
The state of the art in polynomial root finding is to find the eigenvalues of the companion matrix formed from the polynomial coefficients. However, this method is still ill-conditioned for pathological cases, such as the Wilkinson polynomials.

For best numerical stability with extremely large polynomials, the current state of the art is to find the eigenvalues of the colleague matrix for the Chebyshev approximation of the polynomial. See the Chebfun project. This technique works for any piecewise-smooth univariate function, not just polynomials.

Bini and Pan have written extensively on the subject of polynomial root finding, but it's unclear what their work offers over the above techniques, given the complexity of the proposed algorithms and the lack of example implementations.

### Ill-Conditioned Polynomials
Ill-conditioning results from finite precision floating-point operations. The only sure solution for dealing with ill-conditioned polynomials is to use higher-precision arithmetic.

### Jenkins-Traub
This is the most cited black-box polynomial root finder, although it is not superior to the alternative approaches. It has been shown that this algorithm is equivalent to finding eigenvalues of the companion matrix formed from the polynomial coefficients. However, finding the eigenvalues of a matrix has been much more researched, and implementations are prolific and well-tested.