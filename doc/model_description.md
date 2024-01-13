# Model of forest dynamics / forest transition

## Symbols

State variables:
- `P`: Primary forests area
- `S`: Secondary forests area
- `A`: Agriculture area
- `V`: Timber volumes in secondary forest

Control variables:
- `d`: deforested area (primary forest)
- `r`: regeneration area (secondary forest)
- `h`: harvested area (secondary forest) 

Parameters:
- `δ`: discount rate
- `γ`: growth rate of the logistic function in terms of density
- `K`: max density of secondary forest (i.e. carrying capacity)
- `D`: density of the primary forest (constant)
- `L`: total land area 

Abbreviations:
- `PF`: primary forests
- `SF`: secondary forests
- `AG`: agricultural area

## Assumptions
- h1: Primary forests are in stationary state, i.e. their volume variation is only given by clear cutting and its density `D` is fixed
- h2: Harvested area of $PF$ can be allocated to either $SF$ or $AG$
Question: why not area transfer possible from SF to A ? Aside adding complexity to the model, one reason is that in real world land is heterogeneous. If there has been an initial land allocation from PF to SF, it is likely that this choice woudn't change later. Modelling no land tranfer from SF to A is hence a way to implicitly consider heterogeneous land in our simple model.
- h3: Harvesting of $SF$ doesn't change its area, i.e. if clear-cut, the area remains $SF$
- h4: At each moment in time, benefits for the society depend from area of $PF$ (environmental benefits), sum of harvesting volumes of $PF$ and $SF$ (i.e. indifferentiated wood), area of $AG$
- h5: Costs of harvesting $PF$ are inversely proportional to the area of $PF$
- h6: The society can "control" the harvesting of primary forest $d$, the harvesting of secondary forest $h$ and the allocation of harvested primary forest


## Model constraints
- $\dot P_t = - d_t$
- $\dot S_t = r_t$
- $\dot A_t = d_t - r_t$
- $P_t + S_t + A_t = L$ (normalisation)
- $d_t <= P_t$
- $h_t <= S_t$
- $r_t <= d_t$
- $\dot V_{t} = (V_t / S_t) * \gamma * (1 - V_t/ (S_t * K))*S_t -  h_t *(V_t/S_t)$ The first component is the volumes coming from growth, with the logistic equation expressed in terms of density, the second term is the harvested volumes. Note that here we ignore the instantateous growth (of volumes) from the new added area, but that's ok, as it tends to be zero when V=0 in the logistic equation. The expression on the RHS can be simplified to $V_t * \gamma  - V_t^2 * \gamma / (S_t * K) -  h_t *(V_t/S_t)$

## Welfare to maximise (on each time moment)
- $W_{(t)} = b_e\{P\} + b_w\{d,h \frac{V}{S}\} + b_a\{A\}  - c_{p}\{d,P\}  - c_r\{r\} - c_s(h\frac{V}{S})$

## Notes
- Secondary forest is considered an highly aggregated level of forest stands at different ages

## Analytical resolution 
### Current time Hamiltonian

(1) H_C = 
      $b_e\{P_t\} + b_w\{d_t,h_t*V_t/S_t\} + b_a\{A_t\}$
     - $c_{p}\{d_t,P_t\}  - c_r\{r_t\} - c_s\{h_t*V_t/S_t\}$
     + $\lambda_t * (-d_t)$
     + $\phi_t * (r_t)$
     + $\mu_t * (d_t - r_t)$
     + $\rho_t * (V_t  * \gamma * (1 - V_t/ (S_t * K)) -  h_t *(V_t/S_t))$
     + $\sigma_t * (L - P_t - S_t - A_t)$

where:
- $\lambda_t$: co-state variable of the state variable $P_t$
- $\phi_t$: co-state variable of the state variable $S_t$
- $\mu_t$: co-state variable of the state variable $A_t$
- $\rho_t$: co-state variable of the state variable $V_t$
- $\sigma_t$: shadow price of the total land area $L$

### First order conditions

Control variables: $d_t, h_t, r_t$ 

- (2) $\frac{\partial H_C}{\partial d_t} = 0$ ⇒ $\frac{\partial b_w}{\partial d_t} - \frac{\partial c_p}{\partial d_t} - \lambda_t + \mu_t = 0$
- (3) $\frac{\partial H_C}{\partial h_t} = 0$ ⇒ $\frac{\partial b_w}{\partial h_t} - \frac{\partial c_s}{\partial h_t} - \rho_t * \frac{V_t}{S_t} = 0$
- (4) $\frac{\partial H_C}{\partial r_t} = 0$ ⇒ $-\frac{\partial c_r}{\partial r_t} + \phi_t - \mu_t = 0$
- (5) $\frac{\partial H_C}{\partial \sigma_{t}} = 0$ ⇒ $(L - P_t - S_t - A_t) = 0$

#### Interpretations
Equations (2) and (4) can be interpreted as static efficiency rule regarding land transfer decisions.
Eq (2) can be rewritten as $\lambda_t = \frac{\partial b_w}{\partial d_t} - \frac{\partial c_p}{\partial d_t} + \mu_t$, that is the marginal value of the primary forest land, i.e. the cumulative sum of the net benefits arising from the primary forest, must be the same as converting the land to agricultural land, valued at its marginal value $\mu_t$ (in turn the cumulative sum of the agricultural benefits). In the process of converting the land however we do also perform clear-cutting (deforestation), obtaining the net marginal benefits of harvesting primary forest timber $\frac{\partial b_w}{\partial d_t} - \frac{\partial c_p}{\partial d_t}$. Of course, we can also decide to convert the primary forest to secondary forest instead, the marginal value of which is $\phi_t$. As eq. (4) tells us that $\mu_t = -\frac{\partial c_r}{\partial r_t} + \phi_t$ we obtain that $\lambda_t$ is also equal to $\frac{\partial b_w}{\partial d_t} - \frac{\partial c_p}{\partial d_t} -\frac{\partial c_r}{\partial r_t} + \phi_t$. In this case, on top of the net benefits of the wood use from the deforestation we need to account for the regeneration costs.
Finally, equation (3) is a static efficiency condition in terms of harvest intensity of secondary forest. $\rho_t$ is the marginal value of secondary forest volumes. We can use the chain rule to rewrite eq (3) as $\frac{\partial b_w}{\partial hV_t} * \frac{V_t}{S_t} - \frac{\partial c_s}{\partial hV_t} * \frac{V_t}{S_t} = \rho_t * \frac{V_t}{S_t}$, where $hV_t = h_t * \frac{V_t}{S_t}$ are the harvesting _volumes_, to show that the marginal value of the secondary forest volumes are always equal to the net benefits from the harvested timber.


### Equations of motion:

- $\dot \lambda_t = \delta \lambda - \frac{\partial H_C}{\partial P_t} = \delta \lambda - \frac{\partial b_e}{\partial P_t } + \frac{\partial c_p}{\partial P_t } + \sigma_t$ 
- $\dot \phi_t = \delta \phi - \frac{\partial H_C}{\partial S_t} = \delta \phi - \frac{\partial b_w}{\partial S_t} + \frac{\partial c_s}{\partial S_t } - \frac{\rho_t V_t^2 \gamma}{S_t^2 K} - \frac{\rho_t V_t h_t}{S_t^2} + \sigma_t$
- $\dot \mu_t = \delta \mu - \frac{\partial H_C}{\partial A_t} = \delta \mu - \frac{\partial b_a}{\partial A_t } + \sigma_t $
- $\dot \rho_t = \delta \rho - \frac{\partial H_C}{\partial V_t} = \delta \rho - \frac{\partial b_w}{\partial V_t} + \frac{\partial c_s}{\partial V_t} - \rho_t \gamma + \frac{2 \rho_t V_t \gamma}{S_t K} + \frac{\rho_t h_t}{S_t} $




# Numerical implementation


While we did want to apply our model to realistic data, we don't have at this time reliable estimate of the parameters of the benefits and cost functions.
In order to avoid any misinterpretation of our paper as a policy analysis and our results as normative, we therefore decided to apply it to a fictional region, which we call 'Lisarb'.

Lisarb is a large region with significant forest resources. Specifically, it has a land surface of 851.5 Mha, of which 57.0% are natural, primary forests, 1.3% are secondary, plantation forests and 28.0% is agricultural land.

## Functional forms

All costs and benefits are in 10^6 US$. Parameters on the left of the semi-column are endogenous variables computed at each time mpont by the optimization routine, whereas parameters on the right are exogenous and fixed for the whole simulation period.

| Name                                   | Function      |
| -------------------------------------  | ------------- | 
| Environmental benefits                 | `ben_env(F;c1,c2) = c1*F^c2` |
| Agricultural use benefits              | `ben_agr(A;c3,c4) = c3*A^c4` |
| Timber use benefits                    | `ben_timber(S,V,d,h;c5,c6,D) = c5*(d * D + h * V/S)^c6` |
| Carbon benefits                        | `ben_carbon(S,V,d,h,t; D,γ,K,c7,c8) = c7*exp(c8 * t) * (((V)*γ*(1-(V / (S * K) )) - h * (V/S)) - (d * D )  )` |
| Harvesting primary forest costs        | `cost_pfharv(F,d;c9,c10,c11) = c9 * (d*D)^c10 * F^c11` |
| Harvesting secondary forest costs      | `cost_sfharv(S,V,h;c12,c13) = c12 * (h * V/S)^c13` |
| Regeneration of secondary forest costs | `cost_sfreg(r,h;c14,c15) = c14 * (r+h) ^ c15`|

The welfare (the objective of the maximisation) is then the sum of all benefits less the costs at each time point `t`, discounted by `exp(-t*σ)`.

## Parametrization

The `base` scenario assumes the following exogenous parameters. Each individual scenario then overrides one or more of them as specified in the table `Scenarios`.

### Initial values

| Variable | Interpretation | Unit | Value
| ------ | ------- | ----- | ----- |
| F₀ | Init Primary Forest area | M ha | 485.4 |
| S₀ | Init Secondary Forest area | M ha | 11.2 |
| A₀ | Init Agricultural area     | M ha | 238.7 |
| V₀ | Init timber volumes of secondary forests | M m^3 | 3057 |
| d₀ | Init prim forest harvesting area | M ha | 0.49 |
| h₀ | Init sec forest harvesting area | M ha | 0.11 |
| r₀ | Init sec forest regeneration area | M ha | 0.24 |

### Parameters

| Parameter | Interpretation | Unit | Value
| ------ | ------- | ----- | ----- |
| D | Density of the primary forest (constant) | m^3/ha | 244.0 |
| K | Maximum density of the secondary forest (i.e. carrying capacity) | m^3/ha | 600 |
| σ | Discount rate | rate | 0.04 | 
| γ | Growth rate of the logistic function in terms of density | rate | 0.1 |

### Costs and benefits coefficients

The following table provides the values of all the coefficients appearing in the benefit and cost functions. While power coefficients (influencing the degree of curvature of the functions) are simply assumed, the multipliers are calibrated in such a way that the initial values of the relevant costs or benefits assume realistic values. 

| Coef | Interpretation | Value | Origin |
| ---- | -------------- | ----- | ------ |
| c1 | Multiplier of the environmental benefits | 17846 | `ben_env₀ = (0.3*2700) * F₀` |
| c2 | Power of the environmental benefit | 1/2 | assumed |
| c3 | Multiplier of the agricultural benefits | 9330 | `ben_agr₀ = 603 * A₀` |
| c4 | Power of the agricultural benefit | 1/2 | assumed |
| c5 | Multiplier of the timber benefits | 1845 | `ben_timber₀ = 300 * hV₀` |
| c6 | Power of the timber benefits | 1/2 |  assumed |
| c7 | Multiplier of the carbon benefits | 0 | `ben_carbon₀ = 0 * ΔV₀` |
| c8 | Growth rate of the carbon benefits | 0 |  assumed |
| c9 | Multiplier of the PF harvesting costs | 5.55e11 | `cost_pfharv₀ = 10 * hV₀_pf` |
| c10 | Power of the harvesting costs of PF (on harvested area) | 1 |  assumed |
| c11 | Power of the harvesting costs of PF (on area) | -4 |  assumed |
| c12 | Multiplier of the SF harvesting costs | 5 | `cost_sfharv₀ = 5 * hV₀_sf` |
| c13 | Power of the harvesting costs of SF (on harvested area) | 1 |  assumed |
| c14 | Multiplier of the SF regeneration costs | 17.8 | `cost_sfreg₀ = 50 * (r₀+h₀)` |
| c15 | Power of the harvesting costs of SF (on harvested area) | 0 |  assumed |


# Optimization options
| Option | Description | Value |
| ------ | ----------- | ----- |
| optimizer  | Desired optimizer (solver engine) | Ipopt.Optimizer |
| opt_options | Solver options | Dict("max_cpu_time" => 20.0, "print_level" => 3) |
| T | Time horizon (years) | 2000 | 
| ns | Number of points in the time grid | 401 |




## Scenarios