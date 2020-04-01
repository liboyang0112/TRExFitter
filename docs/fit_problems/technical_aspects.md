# Technical aspects of the likelihood creation and minimisation procedure

This is, strictly speaking, not the main point of the tutorial, however, many problems that occur can be understood with the knowledge of some of the technicalities.


## Likelihood creation

The first thing that needs to be understood is that for each NP we provide only two histograms (in many cases only one, but we symmetrise the effect).
With the nominal histogram this leaves only three histograms that are used as an input for each NP.
However, the fit requires a continuous parameter describing systematic uncertainties.
Thus, some interpolation is required.


### Splitting shape and normalisation component

Each source of NP can, in principle, have three different impacts on the nominal distribution, it can affect only the normalisation of the nominal distribution or affect only the shape of the distribution or both shape and normalisation (most common).
If one source affects both the normalisation and shape of a distribution, these two effects are automatically split into pure shape component (no normalisation effect) and pure normalisation effect (no shape effect).
The reason for this split is that the each components undergoes different interpolation strategy.
The shape component uses bin-by-bin linear interpolation, while the normalisation component uses exponential interpolation.

Shape interpolation:

$$
    \sigma_{b,p}(\alpha) = \sigma^{0}_{b,p} + I_{b,p,\text{lin}}, \text{where} I_{b,p,\text{lin}}(b,p,\alpha, I^0_{b,p}, I^+_{b,p},I^-_{b,p}) = \begin{cases} \alpha(I^+_{b,p}-I^0_{b,p}) \quad \alpha \geq 0 \\ \alpha(I^0_{b,p}-I^-_{b,p}) \quad \alpha < 0 \end{cases},
$$

where $I^{+}_{b,p}, I^{-}_{b,p}$ and $I^{0}_{b,p}$ terms represent the expected yields for the systematic up variation, down variation and the nominal prediction for a process $p$ in a bin $b$, respectively.

Normalisation interpolation:

$$
    \eta_p(\alpha) = I_{\text{exp}}, \text{where} I_{\text{exp}}(p,\alpha, I^0_p, I^+_p,I^-_p) = \begin{cases} \left(I^+_p/I^0_p\right)^\alpha \quad &\alpha \geq 0 \\ \left(I^-_p/I^0_p\right)^{-\alpha} \quad &\alpha < 0 \end{cases}.
$$

The exponential interpolation is motivated by preventing the normalisation to go to nonphysical negative values.
Although the normalisation and shape component of each NP are split, they share the same Gaussian (by default) prior, hence only one prior per NP is available.
The exponential interpolation with a Gaussian prior leads to log-normal distribution for the normalisation component.
It should be noted that the exponential interpolation is approximately linear as long as there are no strong pulls, however, with large pulls it deviates from a linear function and leads to non-symmetric effects.
Indeed, this is usually responsible for non-symmetric behaviour of the fit.


## Minimisation procedure

The whole procedure of PL relies on finding the global minimum of the negative logarithm of the likelihood formula, which in many cases present a problem of finding a global minimum in few-hundred dimensional space.
Additionally, the minimisation procedure is not only responsible for the minimum, but it also estimates the correlation of the input parameters (normalisation factors and systematics NPs)!
TRExFitter uses MIGRAD as a minimiser which implements "Davidon–Fletcher–Powell" (DFP) approach.

??? info "DFP References"
    - W. C. Davidon, *Variable Metric Method for Minimization* (1959), [doi:10.2172/4252678](https://doi.org/10.2172%2F4252678),
    - M. J. D. Powell, *Variable Metric Methods for Constrained Optimization* (1983), [doi:10.1007/978-3-642-68874-4_12](https://doi.org/10.1007/978-3-642-68874-4_12),
    - R. Fletcher, *A new approach to variable metric algorithms* (1970), [doi:10.1093/comjnl/13.3.317](https://doi.org/10.1093/comjnl/13.3.317).

The minimisation procedure can be summarised in the following steps

- Start from a given value of parameters, $X$.
- Calculate the gradient (first derivatives) $G$ in the given point, assume the Hessian matrix (matrix of second derivatives) is unity - this assumes the NPs are uncorrelated, which a reasonable starting point physics-wise.
- Perform a linear search, along the direction of the gradient: find $\alpha$ which minimises $F (x -\alpha V \times G)$, where $F$ is the negative logarithm of the likelihood function.
  $V$ is a covariance matrix of the parameters of the fit which is equal to the inverse of the Hessian matrix.
- Correct the covariance matrix V using formulae from DFP.
- Repeat until the estimated distance to minimum (EDM), EDM = $G^T \times V \times G$ is sufficiently small (EDM < 0.001).
  EDM represents vertical distance to the minimum in case of a quadratic function.

The outlined approach allows the minimisation procedure to "climb hills" - to overcome local maxima around local minima.
Additionally, the correlation matrix of the NPs, that can be trivially obtained from the covariance matrix, is available almost as a by-product of the minimisation procedure.
The diagonal elements of the correlation matrix are used to derive the post-fit uncertainties of the NPs and the POI.
The correlation matrix is a symmetric matrix, thus the post-fit uncertainties are also symmetric.
However, it provides a good description of the region around the minimum only if the region is represented well by a quadratic function of the parameters.
To obtain a more accurate estimation of the uncertainty of the POI, the MINOS technique, which takes into account correlations of the parameters and does not rely on the quadratic shape of the logarithm of the likelihood function, is used which may lead to non-symmetric uncertainties on POI or NPs.
You can set the MINOS calculation for parameters in TRExFitter with `UseMinos: paramNames`.
