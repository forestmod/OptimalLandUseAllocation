
module OLUA # OptimalLandUseAllocation

using Ipopt
using InfiniteOpt # To automatically "discretize" a continuous variable, like time in this case
using Markdown
using Revise

# Load exoxes parameters and options
include("default_data.jl")

export luc_model


md"""

Symbols used within the LUC model:

State variables:
- `F`: Primary forests area
- `S`: Secondary forests area
- `A`: Agricultural area
- `V`: Timber volumes in secondary forest

Control variables:
- `d`: deforested area (primary forest)
- `r`: regeneration area (secondary forest)
- `h`: harvested area (secondary forest)
"""

function luc_model(;
    σ           = σ,        # Discount rate 
    K           = K,        # Maximum density of the secondary forest (i.e. carrying capacity), 
    γ           = γ,        # Growth rate of the logistic function in terms of density
    benv_c1     = benv_c1,  # Multiplier of the environmental benefits
    benv_c2     = benv_c2,  # Power of the environmantal benefit
    bagr_c1     = bagr_c1,  # Multiplier of the agricultural benefits
    bagr_c2     = bagr_c2,  # Power of the agricultural benefits
    bwood_c1    = bwood_c1, # Multiplier of the wood-use benefits
    bwood_c2    = bwood_c2, # Power of the wood-use benefits
    chpf_c1     = chpf_c1,  # Multiplier of the harvesting costs of primary forest
    chpf_c2     = chpf_c2,  # Power of the harvesting costs of primary forest (harvested area)
    chpf_c3     = chpf_c3,  # Power of the harvesting costs of primary forest (primary forest area)
    chsf_c1     = chsf_c1,  # Multiplier of the harvesting costs of secondary forest
    chsf_c2     = chsf_c2,  # Power of the harvesting costs of secondary forest
    crsf_c1     = crsf_c1,  # Multiplier of the regeneration costs of secondary forest
    crsf_c2     = crsf_c2,  # Power of the regeneration costs of secondary forest
    D           = D,        # Density of the primary forest  (constant)
    bc_c1       = bc_c1 ,   # Carbon initial price
    bc_c2       = bc_c2,    # Carbon price growth rate

    # Init values...
    F₀          = F₀,       # Initial primary-forest area
    S₀          = S₀,       # Initial secondary forest area
    A₀          = A₀,       # Initial agricultural area
    V₀          = V₀,       # Initial secondary forest volumes
    d₀          = d₀,       # First prim for harvesting
    # Options
    optimizer   = optimizer,  # Desired optimizer (solver)
    opt_options = opt_options, # Optimizer options
    T           = T,     # Time horizont
    ns          = ns       # nNmber of supports on which to divide the time horizon
  ) 

  # Functions...
  discount(t; σ=σ)                                                 = exp(-σ*t) # discount function
  ben_env(F;benv_c1=benv_c1,benv_c2=benv_c2)                       = (benv_c1*F^benv_c2) # Environmental benefits [M$]
  ben_agr(A;bagr_c1=bagr_c1,bagr_c2=bagr_c2)                       = (bagr_c1*A^bagr_c2) # Agricultural use benefits [M$]
  ben_wood(S,V,d,h;bwood_c1=bwood_c1,bwood_c2=bwood_c2,D=D)        = (bwood_c1*(d * D + h * V/S)^bwood_c2) # Wood use benefits [M$]
  ben_carbon(S,V,d,h,t; D=D,γ=γ,K=K,bc_c1=bc_c1,bc_c2=bc_c2)       = bc_c1*exp(bc_c2 * t) * (var_vol_sf(S,V,h,γ=γ,K=K) - (d * D )  ) # Carbon benefits [M$]
  cost_pfharv(F,d;chpf_c1=chpf_c1,chpf_c2=chpf_c2,chpf_c3=chpf_c3) = (chpf_c1 * (d*D)^chpf_c2 * F^chpf_c3) # Harvesting primary forest costs [€]
  cost_sfharv(S,V,h;chsf_c1=chsf_c1,chsf_c2=chsf_c2)               = (chsf_c1 * (h * V/S)^chsf_c2)  # Harvesting secondary forest costs [€]
  cost_sfreg(r,h;crsf_c1=crsf_c1,crsf_c2=crsf_c2)                  = (crsf_c1 * (r+h) ^ crsf_c2)  # Regeneration of secondary forest costs [€]
  

  function welfare(F,S,A,V,d,h,r,t;)
      return (
        ben_env(F;)
      + ben_agr(A;)
      + ben_wood(S,V,d,h;) 
      + ben_carbon(S,V,d,h,t)
      - cost_pfharv(F,d;)
      - cost_sfharv(S,V,h;)
      - cost_sfreg(r,h)
      )
  end

  # Definition of the equations of motion of the state variables
  var_area_pf(d)          = -d
  var_area_sf(r)          = r
  var_area_ag(d,r)        = d - r
  var_vol_sf(S,V,h;γ=γ,K=K) = (V)*γ*(1-(V / (S * K) )) - h * (V/S) # Mm³ See https://en.wikipedia.org/wiki/Logistic_function#In_ecology:_modeling_population_growth (the logistic growth is in terms of density: (V_sf/A_sf)*γ*(1-((V_sf/A_sf) /maxD ))*A_sf -hV_sf )


  # To check only...
  #step = 1
  #ts = 1:step:800
  #acost = 2
  #v0 = 450 * acost
  #hVcost = 3 * acost
  #Vtest = copy(v0)
  #v_by_step = zeros(length(ts))
  #for (it,t) in enumerate(ts)
  #  #println(it)
  #  v_by_step[it] = Vtest
  #  hA_step = hVcost * acost/Vtest 
  #  Vtest += var_vol_sf(acost,Vtest,hA_step*step)*step
  #end
  #Vtest
  #plot(ts, collect(v_by_step[i] for i in 1:length(ts))  )



  solver = optimizer_with_attributes(optimizer, opt_options...)
  m      = InfiniteModel(solver)
  @infinite_parameter(m, t in [0, T], num_supports = ns)

  # Variables declaration...
  @variable(m, F >= 0, Infinite(t), start = F₀ )  # prim forest area
  @variable(m, S >= 0, Infinite(t), start = S₀ )  # sec forest area
  @variable(m, A >= 0, Infinite(t), start = A₀ )  # agr forest area
  @variable(m, V >= 0, Infinite(t), start = V₀ )  # sec vor vol
  @variable(m, d >= 0, Infinite(t), start = d₀ )  # prim for harv area
  @variable(m, r >= 0, Infinite(t), start = r₀ )  # sec for reg area
  @variable(m, h >= 0, Infinite(t), start = h₀  ) # sec for harv area

  # Initial conditions....
  @constraint(m, F(0)  == F₀)
  @constraint(m, S(0)  == S₀)
  @constraint(m, A(0)  == A₀)
  @constraint(m, V(0)  == V₀)
  @constraint(m, d(0)  == d₀)
  #@constraint(m, h(0)  == 0)
  #@constraint(m, r(0)  == 0)


  # Other conditions...
  @constraint(m, (F+S+A)  == (F₀+S₀+A₀) )
  @constraint(m, d <= F)
  @constraint(m, h <= S)
  @constraint(m, r <= d )

  # Final conditions...
  # not specified.
  
  # Addition of the equations of motion of state variable to the problem...
  @constraint(m, dA_pf, deriv(F, t) == var_area_pf(d))
  @constraint(m, dA_sf, deriv(S, t) == var_area_sf(r))
  @constraint(m, dA_ag, deriv(A, t) == var_area_ag(d,r))
  @constraint(m, dV_sf, deriv(V, t) == var_vol_sf(S,V,h))


  @objective(m, Max, integral(welfare(F,S,A,V,d,h,r,t), t, weight_func = discount))

  # Optimisation and retrival of optimal data/time path
  optimize!(m)
  status = termination_status(m)

  F_opt     = value.(F)
  S_opt     = value.(S)
  A_opt     = value.(A)
  V_opt     = value.(V)
  d_opt    = value.(d)
  h_opt    = value.(h)
  r_opt    = value.(r)
  ts       = supports(t)
  opt_obj  = objective_value(m) 
  ben_env_opt = ben_env.(F_opt)
  ben_agr_opt     = ben_agr.(A_opt)
  ben_wood_opt    = ben_wood.(S_opt,V_opt,d_opt,h_opt)
  ben_carbon_opt  = ben_carbon.(S_opt,V_opt,d_opt,h_opt,ts)
  cost_pfharv_opt = cost_pfharv.(F_opt,d_opt,)
  cost_sfharv_opt = cost_sfharv.(S_opt,V_opt,h_opt)
  cost_sfreg_opt  = cost_sfreg.(r_opt,h_opt)
  welfare_opt     = welfare.(F_opt,S_opt,A_opt,V_opt,d_opt,h_opt,r_opt,ts)

  return (F=F_opt, S=S_opt, A=A_opt, V=V_opt, d=d_opt, h=h_opt, r=r_opt,
          obj=opt_obj, support= ts, status=status,
          ben_env = ben_env_opt, ben_agr = ben_agr_opt, ben_wood= ben_wood_opt, ben_carbon = ben_carbon_opt,
          cost_pfharv = cost_pfharv_opt, cost_sfharv = cost_sfharv_opt, cost_sfreg = cost_sfreg_opt,
          welfare = welfare_opt
          )
end


end # end module