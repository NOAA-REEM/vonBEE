github repo: <https://github.com/kholsman/vonBEE>

<img src="Figs/vonBEE.png" style="width:25.0%" />

Last updated: Oct 04, 2024

# Install the package vonBEE

``` r
  missing <- setdiff(c("devtools","usethis"),   installed.packages()[, 1])
  if (length(missing) > 0) 
    install.packages(missing)

  devtools::install_github("NOAA-REEM/vonBEE")
```

# Introduction

-   [ ] Add intro paragraph about climate effects on growth

Water temperature is known to directly impact growth through influencing
metabolic and digestion rates, which often scale exponentially with body
weight and temperature (see Hanson et al., 1997 for an overview).

(based on <https://seananderson.ca/2014/10/17/tmb/>)

# Methods

## Temperature specific weight at age

We modified the generalized formulation of the von Bertalanffy growth
function (VBGF; von Bertalanffy 1938; Pauly 1981; Temming 1994) to
predict temperature-dependent growth by allowing the allometric scaling
parameter *d* to increase with temperature. Essington et al. (2010) and
Holsman and Aydin (2015), and Holsman et al. (2016) describe the
derivation and application of the VBGF towards bioenergetics modeling in
great detail, so we do not repeat it here. Essentially, in this
formulation *d* represents the realized allometric slope of consumption,
which integrates both the direct effect of temperature on consumption
and indirect ecological interactions that scale with temperature and
influence relative foraging rates (see Essington et al., 2010; Holsman
and Aydin, 2015). We fit the VBGF to otolith-based length- and
weight-at-age data (*n* = 21,388, 14,362, and 772, for pollock, Pacific
cod, and arrowtooth flounder, respectively) collected during AFSC Bering
Sea surveys and analyzed at the AFSC (REF).

For each species, we used the vonBT() model to fit a state-space
environmental growth model for fish weight and age data using the base
state-space model (based on Gompertz et al. 20XX):

Eq. 1
$$ ln(\hat{W}{i}) = W^\infty{i} + \frac{1}{(1-d\_{i})}log(1-e^{-K(1-d\_{i})(A_i-t_0)})+\varepsilon_i\sim N(0,\sigma^2{obs})$$
, where

Eq. 2
$$W^\infty\_{i} =(\frac{H}{K})^{\frac{1}{(1 - d\_{i})}} $$

where *t*<sub>0, *i*</sub> is the age at which *W*<sub>*i*</sub> = 0,
*W*<sub>*i*</sub><sup>∞</sup> is the asymptotic mass which can vary by
individual *i* cohort year effects, *H* is the assimilation constant *K*
is the energy loss constant (Essington et al., 2010), and
*ε*<sub>*i*</sub> is a normally and independently distributed random
variable with mean 0 and variance *σ*<sub>*o**b**s*</sub><sup>2</sup>.
Essington et al. (2010) and Holsman and Aydin, (2015) statistically
estimated the *d*, *K* and *H* parameters for various species to
estimate consumption rates. In particular, Holsman and Aydin (2015)
found that the *d* parameter varied between species and regions in
Alaska (USA). We further modified this approach to estimate the
environment impacts growth on growth annually for each year *y* through
a logistic model that includes a vector of annual covariates
(*X*<sub>*c*, *y*</sub>) effects on *d* (which ranges between 0 and 1).
where *α*0<sub>*d*, *i*</sub> and *α*<sub>*d*, *i*, *y*</sub> represent
the mean the *d* consumption parameter intercept and *β*<sub>*c*</sub>
is the coefficient for the residual effect of an environmental variable
on the *d* consumption parameter, such that:

Eq. 3
*d*<sub>*i*</sub> = 1/(1+*e*<sup>−(*U*<sub>*y*</sub> + *β*<sub>0, *y*<sub>*y* − *A*<sub>*i*</sub></sub></sub> + ∑<sub>*c*</sub>(*β*<sub>*c*</sub>\**X*<sub>*c*, *y*</sub>)</sup>)

We chose this formulation based on the empirical relationship between
temperature and consumption, assuming that *d* would capture the
differential effects of temperature on growth, and that waste rates
scale proportionally with weight but do not vary over time with diet or
temperature (i.e. *K* is constant but *d* can vary with temperature).
This formulation allows both the slope and asymptotic limit of growth to
vary with temperature. Similar approaches, with slightly different
modifications to the VBGF, including temperature and prey specific terms
for *d* and *K*, respectively, have been used elsewhere to evaluate
climate impacts on fish growth (e.g., Cheung et al., 2015; Hamre, 2003).

We further modeled the *d* consumption parameter as a state-space model
that estimates random effects on *U*<sub>*y*</sub> as the unobserved
state vector:

Eq. 4
*U*<sub>*y*</sub> = *μ* + *β*<sub>1</sub>*U*<sub>*y* − 1</sub> + *ε*<sub>*y*</sub>
Process error is then modeled as a normal distribution with mean of 0
and standard deviation of *σ*<sub>*p**r**o**c*</sub>:

Eq. 5
*ε*<sub>*y*</sub> ∼ *N*(0,*σ*<sub>*p**r**o**c*</sub><sup>2</sup>)

Optionally, there can be auto regressive (state-space) terms included as
*β*<sub>1</sub> (e.g., if growth the previous year impacts growth in
observation year). Cohort year (*Y*<sub>*y* − *A*<sub>*i*</sub></sub> )
effects can also be included as random or fixed effects via
*β*<sub>0, *y*<sub>*y* − *A*<sub>*i*</sub></sub></sub>. If these terms
are not included the model simplifies to a random (process error) annual
intercept ( *μ*) plus covariate effects. If covariate terms are further
not included and the random process error is not estimated the model
simplifies to the generalized VonB.

Weight at age is then modeled as:

Eq. 6
*l**n*(*W*<sub>*i*</sub>) ∼ *N*(*l**n*(*Ŵ*<sub>*i*</sub>),*σ*<sub>*o**b**s*</sub><sup>2</sup>)
where *σ*<sub>*o**b**s*</sub> is the standard deviation of observation
error (log scale).

<!-- ###------- -->
<!-- Eq. 2 $~~~~~~W_{ij,y}=W_{\infty,iy} (1-e^{(-K_i (1-d_{i,y} )(j-t_{0,i} )) } )^{1/(1-d_{i,y} )} e^\varepsilon$, where $\varepsilon~N(0,\sigma_{d,i}^2 )$ -->
<!-- Eq. 3 $~~~~~~d_{i,y}=e^{(\alpha_{d,i,y}+\alpha0_{d,i}+\beta_{d,i}T_y) }$ -->

# Results

Comparative model fits using average observed bottom temperature are
included below.

-   [ ] Next add in G potential from energetic indices as a covariate.

<figure>
<img src="Figs/model_plots.jpg"
alt="“Fig.1. Comparative model fits to observed weight at age data for Walleye pollock in the EBS”" />
<figcaption aria-hidden="true">“Fig.1. Comparative model fits to
observed weight at age data for Walleye pollock in the EBS”</figcaption>
</figure>

<!-- <!-- !["Fig.2. Temperature model fits to observed weight at age data for Walleye pollock in the EBS"](data/out/BS/WALLEYE%20POLLOCK/Figs/model_Temp.jpg){width="400"} -->

–\>

<figure>
<img src="Figs/model_Temp_byage2.jpg"
alt="“Fig.2. Temperature effects on weight at age" />
<figcaption aria-hidden="true">“Fig.2. Temperature effects on weight at
age</figcaption>
</figure>
