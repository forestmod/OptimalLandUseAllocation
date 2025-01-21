
module OLUA # OptimalLandUseAllocation

using Ipopt
using InfiniteOpt # To automatically "discretize" a continuous variable, like time in this case
using Markdown
using Revise

# Load exoxes parameters and options
include("default_data.jl")

export luc_model, welfare,
       ben_env, ben_agr, ben_wood, ben_carbon,
       cost_pfharv, cost_sfharv, cost_sfreg

ben_env(F;benv_c1=benv_c1,benv_c2=benv_c2)                       = (benv_c1*F^benv_c2) # Environmental benefits [M$]
ben_agr(A;bagr_c1=bagr_c1,bagr_c2=bagr_c2)                       = (bagr_c1*A^bagr_c2) # Agricultural use benefits [M$]
ben_wood(S,V,d,h;bwood_c1=bwood_c1,bwood_c2=bwood_c2,D=D)        = (bwood_c1*(d * D + h * V/S)^bwood_c2) # Wood use benefits [M$] - Here there is an important simplification that I harvest the forest homogeneously (instead of selectively the mature one)
ben_carbon_seq(S,V,d,h,t; D=D,γ=γ,K=K,bc_seq_c1=bc_seq_c1,bc_seq_c2=bc_seq_c2,co2seq=co2seq)  = bc_seq_c1*exp(bc_seq_c2 * t) * (((V)*γ*(1-(V / (S * K) )) - h * (V/S)) - (d * D )  ) * co2seq # Carbon sequestration benefits [M$] - here we don't yeat know var_vol_sf, so rewriting it explicitly
ben_carbon_sub(S,V,d,h,t; D=D,bc_sub_c1=bc_sub_c1,bc_sub_c2=bc_sub_c2,co2sub=co2sub)  = bc_sub_c1*exp(bc_sub_c2 * t) * (d * D + h * V/S)  * co2sub # Carbon substitution benefits [M$]
cost_pfharv(F,d;chpf_c1=chpf_c1,chpf_c2=chpf_c2,chpf_c3=chpf_c3) = (chpf_c1 * (d*D)^chpf_c2 * F^chpf_c3) # Harvesting primary forest costs [€]
cost_sfharv(S,V,h;chsf_c1=chsf_c1,chsf_c2=chsf_c2)               = (chsf_c1 * (h * V/S)^chsf_c2)  # Harvesting secondary forest costs [€]
cost_sfreg(r,h;crsf_c1=crsf_c1,crsf_c2=crsf_c2)                  = (crsf_c1 * (r+h) ^ crsf_c2)  # Regeneration of secondary forest costs [€]

co2_seq(S,V,d,h; D=D,γ=γ,K=K,co2seq=co2seq) = (((V)*γ*(1-(V / (S * K) )) - h * (V/S)) - (d * D )  ) * co2seq # CO2eq sequestered in forest resourses that year [tons of CO2 eq]
co2_sub(S,V,d,h; D=D,co2sub=co2sub) = (d * D + h * V/S)  * co2sub # CO2eq substituted by harvested timber that year [tons of CO2 eq]

function welfare(F,S,A,V,d,h,r,t;D=D,γ=γ,K=K,
                co2seq=co2seq,co2sub=co2sub,
                benv_c1=benv_c1,benv_c2=benv_c2,
                bagr_c1=bagr_c1,bagr_c2=bagr_c2,
                bwood_c1=bwood_c1,bwood_c2=bwood_c2,
                bc_seq_c1=bc_seq_c1,bc_seq_c2=bc_seq_c2,
                bc_sub_c1=bc_sub_c1,bc_sub_c2=bc_sub_c2,
                chpf_c1=chpf_c1,chpf_c2=chpf_c2,chpf_c3=chpf_c3,
                chsf_c1=chsf_c1,chsf_c2=chsf_c2,
                crsf_c1=crsf_c1,crsf_c2=crsf_c2
  )
  return (
    ben_env(F;benv_c1=benv_c1,benv_c2=benv_c2)
  + ben_agr(A;bagr_c1=bagr_c1,bagr_c2=bagr_c2)
  + ben_wood(S,V,d,h;bwood_c1=bwood_c1,bwood_c2=bwood_c2,D=D) 
  + ben_carbon_seq(S,V,d,h,t;D=D,γ=γ,K=K,bc_seq_c1=bc_seq_c1,bc_seq_c2=bc_seq_c2,co2seq=co2seq)
  + ben_carbon_sub(S,V,d,h,t;D=D,bc_sub_c1=bc_sub_c1,bc_sub_c2=bc_sub_c2,co2sub=co2sub)
  - cost_pfharv(F,d;chpf_c1=chpf_c1,chpf_c2=chpf_c2,chpf_c3=chpf_c3)
  - cost_sfharv(S,V,h;chsf_c1=chsf_c1,chsf_c2=chsf_c2)
  - cost_sfreg(r,h;crsf_c1=crsf_c1,crsf_c2=crsf_c2)
  )
end 

md"""

Compute land use optimizing welfare.

A note on the time dimension:
state_var[t] = state_var[t-1] + flow_var[t]
state_var[t=0] is the current state
The welfare optimization relates hence to times [t=1,t=T], i.e. [1,T], but the variables in the model relate to time [t=0,t=T] as we need the state variables at time t=0 (and these are fixed, as well for the flow faviables t=0 - that don't influence the model)
"""

function luc_model(;
    σ           = σ,        # Discount rate 
    K           = K,        # Maximum density of the secondary forest (i.e. carrying capacity), 
    γ           = γ,        # Growth rate of the logistic function in terms of density
    co2seq      =  co2seq,  # Coefficient m³ wood resourses -> ton CO2eq sequestered
    co2sub      =  co2sub,  # Coefficient m³ harvested timber -> ton CO2eq substituted
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
    bc_seq_c1       = bc_seq_c1 ,   # Carbon (seq) initial price
    bc_seq_c2       = bc_seq_c2,    # Carbon (seq) price growth rate
    bc_sub_c1       = bc_sub_c1 ,   # Carbon (sub) initial price
    bc_sub_c2       = bc_sub_c2,    # Carbon (sub) price growth rate

    # Init values...
    F₀          = F₀,       # Initial primary-forest area
    S₀          = S₀,       # Initial secondary forest area
    A₀          = A₀,       # Initial agricultural area
    V₀          = V₀,       # Initial secondary forest volumes
    d₀          = d₀,       # Initial prim for harvesting
    h₀          = h₀,       # Initial sec for harvesting
    r₀          = r₀,       # Initial sec for regeneration
    # Options
    optimizer   = optimizer,   # Desired optimizer (solver)
    opt_options = opt_options, # Optimizer options
    T           = T,     # Time horizont
    ns          = ns     # nNmber of supports on which to divide the time horizon
  ) 

  # Functions...
  discount(t; σ=σ)                                                 = t == 0 ? 0.0 : exp(-σ*t) # We don't consider welfare from t=0 where flow variables are not influential

  # Definition of the equations of motion of the state variables
  var_area_pf(d)          = -d
  var_area_sf(r)          = r
  var_area_ag(d,r)        = d - r
  var_vol_sf(S,V,h;γ=γ,K=K) = (V)*γ*(1-(V / (S * K) )) - h * (V/S) # Mm³ See https://en.wikipedia.org/wiki/Logistic_function#In_ecology:_modeling_population_growth (the logistic growth is in terms of density: (V_sf/A_sf)*γ*(1-((V_sf/A_sf) /maxD ))*A_sf -hV_sf )
  # Here there is an important simplification that I harvest the forest homogeneously (instead of selectively the mature one)

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
  @constraint(m, h(0)  == h₀)
  @constraint(m, r(0)  == r₀)
  #@constraint(m, r(0)  == 0)
  >

  # Other conditions...
  #@constraint(m, tot_land, (F+S+A)  == (F₀+S₀+A₀) )
  @constraint(m, d <= F)
  @constraint(m, h <= S)
  #@constraint(m, h == 0)
  #@constraint(m, r == 0)
  @constraint(m, r <= d )

  # Final conditions...
  # not specified.
  
  # Addition of the equations of motion of state variable to the problem...
  @constraint(m, dA_pf, deriv(F, t) == var_area_pf(d))
  @constraint(m, dA_sf, deriv(S, t) == var_area_sf(r))
  @constraint(m, dA_ag, deriv(A, t) == var_area_ag(d,r))
  @constraint(m, dV_sf, deriv(V, t) == var_vol_sf(S,V,h))


  @objective(m, Max, integral(
    welfare(F,S,A,V,d,h,r,t;
            D=D,γ=γ,K=K,
            co2seq=co2seq,co2sub=co2sub,
            benv_c1=benv_c1,benv_c2=benv_c2,
            bagr_c1=bagr_c1,bagr_c2=bagr_c2,
            bwood_c1=bwood_c1,bwood_c2=bwood_c2,
            bc_seq_c1=bc_seq_c1,bc_seq_c2=bc_seq_c2,
            bc_sub_c1=bc_sub_c1,bc_sub_c2=bc_sub_c2,
            chpf_c1=chpf_c1,chpf_c2=chpf_c2,chpf_c3=chpf_c3,
            chsf_c1=chsf_c1,chsf_c2=chsf_c2,
            crsf_c1=crsf_c1,crsf_c2=crsf_c2),
    t, weight_func = discount)
  )

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
  
  pF       = .- dual.(dA_pf)
  pS       = .- dual.(dA_sf)
  pA       = .- dual.(dA_ag)
  pV       = .- dual.(dV_sf)
  #pTL      = .- dual.(tot_land)

  ts       = supports(t)
  opt_obj  = objective_value(m) 
  ben_env_opt = ben_env.(F_opt)
  ben_agr_opt     = ben_agr.(A_opt)
  ben_wood_opt    = ben_wood.(S_opt,V_opt,d_opt,h_opt)
  ben_carbon_seq_opt  = ben_carbon_seq.(S_opt,V_opt,d_opt,h_opt,ts)
  ben_carbon_sub_opt  = ben_carbon_sub.(S_opt,V_opt,d_opt,h_opt,ts)
  cost_pfharv_opt = cost_pfharv.(F_opt,d_opt,)
  cost_sfharv_opt = cost_sfharv.(S_opt,V_opt,h_opt)
  cost_sfreg_opt  = cost_sfreg.(r_opt,h_opt)
  welfare_opt     = welfare.(F_opt,S_opt,A_opt,V_opt,d_opt,h_opt,r_opt,ts)
  co2_seq_opt     = co2_seq.(S_opt,V_opt,d_opt,h_opt)
  co2_sub_opt     = co2_sub.(S_opt,V_opt,d_opt,h_opt)

  return (F=F_opt, S=S_opt, A=A_opt, V=V_opt, d=d_opt, h=h_opt, r=r_opt,
          obj=opt_obj, support= ts, status=status,
          ben_env = ben_env_opt, ben_agr = ben_agr_opt, ben_wood= ben_wood_opt, ben_carbon_seq = ben_carbon_seq_opt, ben_carbon_sub = ben_carbon_sub_opt,
          cost_pfharv = cost_pfharv_opt, cost_sfharv = cost_sfharv_opt, cost_sfreg = cost_sfreg_opt,
          welfare = welfare_opt, co2_seq = co2_seq_opt, co2_sub = co2_sub_opt, 
          pF=pF, pS=pS, pA=pA, pV=pV,
          # pTL=pTL
          )
end


end # end module