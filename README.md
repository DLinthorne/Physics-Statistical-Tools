# Physics Statistical and Analysis Tools

Analysis tools for physics computations and statistics

## `ROOT`Statistical/ML Methods

The following folders focus using CERN's `ROOT C++ toolkit to conduct statistical tests, simulations, and analysis:

- [bayesian_routines](bayesian_routines/) shows example macros and methods on bayesian comparisons of parameter,
estimation and confidence intervals and limits. A local P-value search is done on ATLAS Higgs Boson data.
- [central_limit_theorem](central_limit_theorem/) shows how the central limit theorem can be realized using simple
parametric distributions.
- [chi2_contours](chi2_contours/) Calculates the two dimensional contour intevals using the Neyman-Pearson convention.
- [chi2_methods](chi2_methods/) Calculates the $chi^2$ distribution through exact methods and monte-carlo. An examples of $chi^2$ fitting methods is shown.
- [estimator_statistics](estimator_statistics/) Calculates simple estimators/moments (mean, variances, bias, ...)
- [particle_classifier](particle_classifier/) A full analytic investigation of hypothesis classifying using a fully implemented fisher linear discriminator. Real particle data from CERN's ATLAS experiment is used as an example.
- [resonances](resonances/) Shows the effects of convolving detector noise on bump-hunt resonances (Breit-Wigner).
- [simple_monte_carlo](simple_monte_carlo/) Showing a Monte Carlo geometry acceptance calculation (dart on a sphere example)

