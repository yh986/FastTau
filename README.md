# FastTau
R project for computing Kendall Tau quickly (O nlogn) while handling incomplete data appropriately

---
**Introduction**
Correlations are important in science! The standard Pearson correlation calculation does not work well in non-normally-distributed data; Kendall's tau is a measure of correlation which is resistant to non-normality.

However, the R base function cor(X, Y, method = c("kendall")) is very difficult to use on large data sets because it runs at speed O(n^2). Furthermore, cor() is notorious for handling missing values poorly when correlating sets of vectors, giving all sorts of unexpected values depending on the option used.

There exists an algorithm by Christensen, D. (2005) to compute Kendall tau at O(nlogn), but existing implementations in R packages handle missing values or ties poorly (from the ones I've found, at least). So, here is an implementation to do it quickly and with support for missing values, simply, it'll return NaN if a pair of vectors cannot have tau computed.

Christensen, D. Fast algorithms for the calculation of Kendall’s τ. Computational Statistics 20, 51–62 (2005). https://doi.org/10.1007/BF02736122

---
**The files**
Tau.R has the main function Tau(A, B) to correlate between vectors or matrices A and B.
Rcpp's sourcecpp function compiles the C++ file, IntegerTau.cpp, which implements Christensen (2005)'s fast Tau algorithm.
