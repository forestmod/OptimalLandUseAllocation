# # Forest dynamics numerical model

# Literate.markdown("harvest_dynamics.jl", "."; flavor=Literate.CommonMarkFlavor(), execute=true) #src
# p{break-inside: avoid} #src

# ### Setting up the environment / packages...
cd(@__DIR__)
using Pkg
Pkg.activate(".")
#Pkg.resolve()
#Pkg.instantiate()

using Plots
using Plots.PlotMeasures
using Markdown
using Revise
using Test

# ### Load the model structure (function "luc_model")
includet("model.jl")
using .OLUA


# ------------------------------------------------------------------------------
# ### Compute the "base" optimisation...
base     = luc_model() # compute the optimization for the "base" scenario
ix       = 2:(length(base.support)-320) # index for plotting (we discard the distant future and the first year that is not part of the optimization)
times    = base.support[ix] 
ntpoints = length(times)
step     = times[end] / (ntpoints)  # We choosen to discretize every 5 years

# Computing the indexes up to year 200
ix200    = findfirst(t -> t >= 200, times)
ixs200   = ix[1:ix200]
times200 = base.support[ixs200]
ntpoints200 = length(times200)

# Computing the indexes up to year 100
ix100    = findfirst(t -> t >= 100, times)
ixs100   = ix[1:ix100]
times100 = base.support[ixs100]
ntpoints100 = length(times100)

# ------------------------------------------------------------------------------
# ### Model validation...
# As the overall model is analytically too complex, we opted by a validation "by parts", that is, we set the parameters of the models such that some components are disables (for example setting wood benefits to zero disable harvesting) so that the remaining parts have predictable characteristics.
# We performed 4 validation exercises: in the first one we test the effect of the specific time discretization emploied; in the second test we verify that with 0 discount rate, the SF is harvested such to reach in the long run the MSY; in the third test we verify that when we have benefits only from harvesting the PF we end up to a "eating the cake" model, and in particular the marginal penefit of wood (net price) follows the Hotelling's rule. Finally in the last validation test we check that without harvesting the SF follows the Verhulst (logistic) model.
# Each validation is performed with a in-code test that the results follow the assumptions and by plotting a chart of the relevant outputs.  

# ------------------------------------------------------------------------------
# #### Test 1: discretization choice (using less or more time points) doesn't influence much the results
out_dense    = luc_model(ns=1001,opt_options = Dict("max_cpu_time" => 60.0))
ix_dense     = 2:(length(out_dense.support)-800) # index for plotting
times_dense  = out_dense.support[ix_dense] # every 2 years
out_sparse   = luc_model(ns=201)
ix_sparse    = 2:(length(out_sparse.support)-160) # index for plotting
times_sparse = out_sparse.support[ix_sparse] # every 10 years
@testset "Number of discretization points" begin
    @test isapprox(base.r[findfirst(t-> t==80,times)],out_dense.r[findfirst(t-> t==80,times_dense)],atol=0.05)
    @test isapprox(base.r[findfirst(t-> t==80,times)],out_sparse.r[findfirst(t-> t==80,times_sparse)],atol=0.05)
    @test isapprox(base.F[findfirst(t-> t==80,times)],out_dense.F[findfirst(t-> t==80,times_dense)],rtol=0.05)
    @test isapprox(base.F[findfirst(t-> t==80,times)],out_sparse.F[findfirst(t-> t==80,times_sparse)],rtol=0.05)
end

# Graphically...
plot(times, base.r[ix], xlabel="years", label="Base time point density (5 y)", title="SF reg area (flow) under different time densities");
plot!(times_dense, out_dense.r[ix_dense], label="Dense time point density (2 y)");
plot!(times_sparse, out_sparse.r[ix_sparse], label="Sparse time point density (10 y)")
#-
plot(times, base.F[ix], xlabel="years", label="Base time point density (5 y)", title="Forest prim area (stock) under different time densities");
plot!(times_dense, out_dense.F[ix_dense], label="Dense time point density (2 y)");
plot!(times_sparse, out_sparse.F[ix_sparse], label="Sparse time point density (10 y)")

# **Take home**
# Altought for stock variables being cumulative,we have some differences, discretization doesn't significantly influence the nature of the results. All our comparisions of the scenario with base are made using the same time discretization.

# ------------------------------------------------------------------------------
# #### Test 2: setting no interest rate leads to MSY in secondary forest (density = half the maximum density)

out = luc_model(σ=0)

Dsf = out.V ./ out.S # Secondary forest density
@testset "No discount leads to MSY in SF" begin
  @test isapprox(Dsf[findfirst(t-> t==100,times)], OLUA.K/2, rtol=0.01)
  @test isapprox(Dsf[findfirst(t-> t==200,times)], OLUA.K/2, rtol=0.01)
end

# Graphically...
plot(times, out.V[ix] ./ out.S[ix], lab = "Secondary forest density", xlabel="years", linecolor="darkseagreen3", title= "Secondary forest density under no discount");
plot!(times, fill(OLUA.K,ntpoints), lab="Max density")

# **Take home**
# Under no discount, intuitivelly the long term management of the forest is such that it allows the forest volumes to reache the stocks associated to the MSY (that, in the employed logistic model are at half the carrying capacity) in order to harvest the maximum possible timber.

# ------------------------------------------------------------------------------
# #### Test 3: setting cost and benefits parameters in order to have no harvesting in SF and no benefits other than timber leads to the "eating the cake" model / Hotelling rule

# The first version is with only benefits, the second one with some harvesting costs
out1 = luc_model(benv_c1=0.0,bagr_c1=0,chpf_c3=0,chpf_c1=0,crsf_c1=1000,chsf_c1=10000,γ=0) # with only timber benefits
out2 = luc_model(benv_c1=0.0,bagr_c1=0,chpf_c3=0,chpf_c1=10,chpf_c2=1.2,crsf_c1=1000,chsf_c1=10000,γ=0)# d₀=1.13) # with timber benefits and costs

dbw_dV(V,c1=OLUA.bwood_c1,c2=OLUA.bwood_c2) = c1*c2*V^(c2-1) # marginal benetits
dcw_dV(V,c1=10,c2=1.2) = c1*c2*V^(c2-1) # marginal costs

hv_1 = @. out1.d * OLUA.D + out1.h * out1.V / out1.S
hv_2 = @. out2.d * OLUA.D + out2.h * out2.V / out2.S

np1 = dbw_dV.(hv_1)                  # (net, as no costs) price
np2 = dbw_dV.(hv_2) .- dcw_dV.(hv_2) # net price

r1a = log(np1[28]/np1[27])/step     # growth rate of the price for 1 y
r1b = log(np1[25]/np1[20])/(step*5) # growth rate of the price for 5 y
r2a = log(np2[28]/np2[27])/step 
r2b = log(np2[25]/np2[20])/(step*5) 

@testset "Checking net price growth follows Hotelling rule" begin
  @test isapprox(r1a,OLUA.σ,rtol=0.01)
  @test isapprox(r1b,OLUA.σ,rtol=0.01)
  @test isapprox(r2a,OLUA.σ,rtol=0.05)
  @test isapprox(r2b,OLUA.σ,rtol=0.05)
end

# Graphically...
plot(times, np1[ix], lab = "Considering only benefits", xlabel="years", title= "Net price of the timber resource")
plot!(times, np2[ix], lab = "Considering harvesting costs")

# **Take home**
# The growth rate is the same, just it looks lower
# The 10k€ choke price that can be seen in the figure is only an artifact derived by the fact that we didn't change the structure of the model for this validation exercise, but just set the secondary forest harvesting costs at 10k€ to effectively inibit any sf harvesting in that price range. 

# ------------------------------------------------------------------------------
# ### Test 4: setting no Sf area change and no harvesting leads SF volumes to the Verhulst model

out = luc_model(bwood_c1=0.01,h₀=0,bagr_c1=0.0001)
Dsf = out.V ./ out.S 

@testset "No wood benefits lead no harvesting, no area change and SF following the Verhulst model" begin
  @test isapprox(Dsf[findfirst(t-> t==100,times)],OLUA.K,rtol=0.01)
  @test isapprox(Dsf[findfirst(t-> t==200,times)],OLUA.K,rtol=0.01)
  @test isapprox(out.V[findfirst(t-> t==200,times)],OLUA.K*OLUA.S₀,rtol=0.01)
end

# Graphically...
plot(times, Dsf[ix],   lab = "D: secondary forest density", linecolor="darkseagreen3", xlabel="years", title= "Sec. forest density under no harv and no reg");
plot!(times, fill(OLUA.K,ntpoints), lab="Max density")

#-

plot(times,out.V[ix], label="V: model output logistic", xlabel="years", title= "Sec. forest volumes under no harv and no reg");

# Manual computation of the V increase (by loop)
var_vol_test(S,V;γ=OLUA.γ,K=OLUA.K) = (V)*γ*(1-(V / (S * K) ))
ts2 = 0:step:times[end]
Vtest = copy(OLUA.V₀)
v_by_step = zeros(length(ts2))
for (it,t) in enumerate(ts2)
  v_by_step[it] = Vtest
  global Vtest += var_vol_test(OLUA.S₀,Vtest)*step
end
Vtest
plot!(ts2, collect(v_by_step[i] for i in 1:length(ts2)), label="V: discrete steps logistic" );

# Computation of the volumes using the continuous logistic function
logistic(x;k,r,x0) = k/(1+((k-x0)/x0)*exp(-r*x))
plot!(x->logistic(x,k=OLUA.K*OLUA.S₀,r=OLUA.γ,x0=OLUA.V₀),0,times[end], label = "V: true logistic funcion")

# ------------------------------------------------------------------------------
# ### Test 5: verification of the steady state

# Verification of h equilibrium.. in the steady state (that seems to be reached in the model) we must have constant
# volumes, where the growth is balanced by the harvesting

rhs = @. OLUA.γ*(base.S * OLUA.K - base.V) / OLUA.K
lhs = base.h

plot(times, rhs[ix]);
plot!(times,base.h[ix])

@test isapprox(base.h[findfirst(t-> t==200,times)], rhs[findfirst(t-> t==200,times)], rtol=0.001)


# ------------------------------------------------------------------------------
# ### Base Case Analysis

# Areas..
plot(times,  base.F[ix], lab = "F: primary forest area", linecolor="darkgreen", xlabel="years", title="Land Areas (base scen)");
plot!(times, base.S[ix], lab = "S: secondary forest area", linecolor="darkseagreen3");
plot!(times, base.A[ix], lab = "A: agricultural area",linecolor="sienna")

# Area transfers..
plot(times,  base.d[ix] .- base.r[ix], lab = "F -> A", linecolor="darkgreen", xlabel="years", title="Land Areas transfers (base scen)");
plot!(times, base.r[ix], lab = "F -> S", linecolor="darkseagreen3")
# plot!(times, base.a[ix], lab = "S -> A",linecolor="sienna") # This only with model exteded version ! #src

# Shadow prices..
plot(times[2:end],  base.pF[ix[2:end]], lab = "pF: primary forest area price", linecolor="darkgreen", xlabel="years", title="Land shadow prices (base scen)")
plot!(times[2:end], base.pS[ix[2:end]], lab = "pS: secondary forest area price", linecolor="darkseagreen3")
plot!(times[2:end], base.pA[ix[2:end]], lab = "pA: agricultural area price",linecolor="sienna")
plot!(times[2:end], base.pV[ix[2:end]], lab = "pV: secondary forest timber volumes price",linecolor="darkgrey")
#plot!(times[2:end], base.pTL[ix[2:end]], lab = "pTL: shadow price of the total land",linecolor="darkgrey")




# Note that the s.p. of primary forest is negative, because we are constraining our model to a fixed area, so 1 ha more of primary prodict means 1 ha less of secondary forest or agricultural area

# Volumes...
plot(times, base.F[ix] .* OLUA.D, lab = "V (pf): primary forest volumes",linecolor="darkgreen", title="Growing Volumes (base scen)", xlabel="title");
plot!(times, base.V[ix], lab = "V (sf): secondary forest volumes", linecolor="darkseagreen3");
plot!(times, base.F[ix] .* OLUA.D .+ base.V[ix],  lab = "TOT V: total forest volumes", linecolor="brown")

# Welfare analysis...
plot(times,  base.ben_env[ix], lab = "Environmental benefits", linecolor="darkgreen", title="Welfare balance (base scen)"); 
plot!(times, base.ben_agr[ix], lab = "Agr benefits",linecolor="sienna") ;
plot!(times, base.ben_wood[ix], lab = "Wood use benefits", linecolor="darkseagreen3");
plot!(times, base.ben_carbon_seq[ix], lab = "Carbon benefits (seq)", linecolor="grey");
plot!(times, base.ben_carbon_sub[ix], lab = "Carbon benefits (sub)", linecolor="darkgrey");
plot!(times, .- base.cost_pfharv[ix], lab = "PF harvesting costs", linecolor="darkgreen");
plot!(times, .- base.cost_sfharv[ix], lab = "SF harvesting costs", linecolor="darkseagreen3");
plot!(times, .- base.cost_sfreg[ix], lab = "SF regeneration costs",linecolor="sienna")
#-
plot(times, base.welfare[ix], lab = "Total welfare (base scen)", ylims=(0,1E6) )

# ------------------------------------------------------------------------------
# ### Scenario analysis
#
# TODO: Divide in 3 sets of scenarios
# A: Environmental analysis : `no_env_ben`, with_carb_ben_1, with_carb_ben_3, with_carb_ben_3, with_carb_ben_grp, restricted_pf_harv
# B: Price / market analysis: incr_timber_demand, lower_disc_rate, 
# C: CC efffect: cc_effect_pf, cc_effect_sf, cc_effect_ag
# We decouple CC effect with different anaysis on the specific ecosystem that we cosider mostly impacted by cc, at the time pf, sf or ag


# ------------------------------------------------------------------------------
# #### Scen 1: `no_env_ben`: PF environmental benefits not considered

out = luc_model(benv_c1=0.0)

plot(times, base.F[ix], lab = "F - base", linecolor="darkgreen", linewidth = 2, xlabel="years", title="Land allocation under no env benefits");
plot!(times, out.F[ix], lab = "F - no_env_ben", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - no_env_ben", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - no_env_ben", linecolor="sienna", ls=:dot)

harv_base    = copy(base.h)
harv_no_benv = copy(out.h)

# **Take home**
# Non considering environmental benefits would, as expect, largelly faster the deforestation of primary forests. Interesting, the increased freed areas would largelly be allocated to agriculture, as the deforestation would imply larger supply of timber that would not benefits the secondary forests (our timber benefits function is concave, and the timber from primary and secondary forest is a completelly homogeneous product)

# ------------------------------------------------------------------------------
# #### Scen 2: Carbon benefits (storage) accounted for
# ##### `with_carb_ben_1`:  maxD PF << maxD sf
out = luc_model(bc_seq_c1=100.0)

plot(times, base.F[ix], lab = "F - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith carb ben (maxD PF << maxD SF)");
plot!(times, out.F[ix], lab = "F - with_carb_ben_1", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - with_carb_ben_1", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - with_carb_ben_1", linecolor="sienna", ls=:dot)
#-
plot(times, base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times, out.V[ix] ./ out.S[ix],   lab = "with_carb_ben_1", linecolor="darkseagreen3", ls=:dot)

harv_with_carb_benv = copy(out.h)
plot(times, harv_base[ix], lab = "base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Harvesing in SF (h) under different scenarios")
plot!(times, harv_no_benv[ix], lab = "no_env_ben", linewidth = 2, linecolor="darkred")
plot!(times, harv_with_carb_benv[ix], lab = "considering_carbon_seq", linewidth = 2, linecolor="darkblue")
savefig("h_by_scenario.png")

# #####  `with_carb_ben_2`:  maxD PF < maxD sf

out1 = luc_model(bc_seq_c1=000.0, D=OLUA.K*0.8)
out2 = luc_model(bc_seq_c1=100.0, D=OLUA.K*0.8)

plot(times, out1.F[ix], lab = "F - with_carb_ben_2a", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith carb ben (maxD PF < maxD sf)", legend=:topright);
plot!(times, out2.F[ix], lab = "F - with_carb_ben_2b", linecolor="darkgreen", ls=:dot);
plot!(times, out1.S[ix], lab = "S - with_carb_ben_2a", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out2.S[ix], lab = "S - with_carb_ben_2b", linecolor="darkseagreen3", ls=:dot);
plot!(times, out1.A[ix], lab = "A - with_carb_ben_2a", linewidth = 2, linecolor="sienna");
plot!(times, out2.A[ix], lab = "A - with_carb_ben_2b", linecolor="sienna", ls=:dot)
#-
plot(times, base.V[ix]./ base.S[ix],   lab = "with_carb_ben_2a", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times, out.V[ix] ./ out.S[ix],   lab = "with_carb_ben_2b", linecolor="darkseagreen3", ls=:dot)

# #####  `with_carb_ben_3`:  maxD PF W maxD sf

out1 = luc_model(bc_seq_c1=000.0, D=OLUA.K*1.2)
out2 = luc_model(bc_seq_c1=100.0, D=OLUA.K*1.2)

plot(times, out1.F[ix], lab = "F - with_carb_ben_3a", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith carb ben (maxD PF > maxD sf)", legend=:topright);
plot!(times, out2.F[ix], lab = "F - with_carb_ben_3b", linecolor="darkgreen", ls=:dot);
plot!(times, out1.S[ix], lab = "S - with_carb_ben_3a", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out2.S[ix], lab = "S - with_carb_ben_3b", linecolor="darkseagreen3", ls=:dot);
plot!(times, out1.A[ix], lab = "A - with_carb_ben_3a", linewidth = 2, linecolor="sienna");
plot!(times, out2.A[ix], lab = "A - with_carb_ben_3b", linecolor="sienna", ls=:dot)
#-
plot(times, base.V[ix]./ base.S[ix],   lab = "with_carb_ben_3a", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times, out.V[ix] ./ out.S[ix],   lab = "with_carb_ben_3b", linecolor="darkseagreen3", ls=:dot)

# **Take home**
# Interesting, because in our default parameterization the maximum density of secondary forest is (not much) higher than primary forests, considering carbon value for timber sequestration would accellerate deforestation in favour of secondary forests
# In other words, as feared by some, carbon payments for forest sequestration could indeed favour monospecific plantations
# However this depends on the relative maximum density between PF and SF (and the realtive  profitability): as max D on PF approach or overpass SF the opposite is true, and less deforestation of PF happens as remain convenient to keep the carbon in the PF

# ------------------------------------------------------------------------------
# #### Scen 3: `with_carb_ben_grp`: Carbon benefits (storage) accounted for (with growing carbon price)
out = luc_model(bc_seq_c1=100.0,bc_seq_c2=(OLUA.σ-0.005)) # setting bc_seq_c2=OLUA.σ doesn't solve

plot(times, base.F[ix], lab = "F - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith increasing carb benefits");
plot!(times, out.F[ix], lab = "F - with_carb_ben_grp", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - with_carb_ben_grp", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - with_carb_ben_grp", linecolor="sienna", ls=:dot)
#-
plot(times, base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times, out.V[ix] ./ out.S[ix],   lab = "with_carb_ben_grp", linecolor="darkseagreen3", ls=:dot)

# **Take home**
# Not what I did expected. I did expect that, becasue when you have to repay it is more expensive, it doesn't become appropriate to storage carbon.
# Instead it looks like it is just an increase version of the `with_carb_ben_2` scenario

# ------------------------------------------------------------------------------
# #### Scen 5: `restricted_pf_harv`: Restricted primary forest harvesting
out = luc_model(chpf_c1=OLUA.chpf_c1*5)

plot(times, base.F[ix], lab = "F - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith PF harv restrictions");
plot!(times, out.F[ix], lab = "F - restricted_pf_harv", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - restricted_pf_harv", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - restricted_pf_harv", linecolor="sienna", ls=:dot)
#-
plot(times, base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times, out.V[ix] ./ out.S[ix],   lab = "restricted_pf_harv", linecolor="darkseagreen3", ls=:dot)
#-
plot(times, base.h[ix] .* base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF hV");
plot!(times, out.h[ix] .* out.V[ix] ./ out.S[ix],   lab = "restricted_pf_harv", linecolor="darkseagreen3", ls=:dot)



# ------------------------------------------------------------------------------
# #### Scen 4: `incr_timber_demand``: Increased benefits of timber resources
out = luc_model(bwood_c1=OLUA.bwood_c1*5)

plot(times, base.F[ix], lab = "F - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith increasing timber demand");
plot!(times, out.F[ix], lab = "F - incr_timber_demand", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - incr_timber_demand", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - incr_timber_demand", linecolor="sienna", ls=:dot)
#-
plot(times, out.h[ix] .* out.V[ix] ./ out.S[ix], lab = "hV", linewidth = 2, linecolor="darkseagreen3", xlabel="years", title="Harvested volumes (increased wood demand)")
plot(times, base.h[ix] .* base.V[ix] ./ base.S[ix], lab = "hV", linewidth = 2, linecolor="darkseagreen3", xlabel="years", title="Harvested volumes (increased wood demand)")
#-
plot(times[2:end], base.V[ix[2:end]]./ base.S[ix[2:end]],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times[2:end], out.V[ix[2:end]] ./ out.S[ix[2:end]],   lab = "incr_timber_demand", linecolor="darkseagreen3", ls=:dot)
#-
plot(times, base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times, out.V[ix] ./ out.S[ix],   lab = "incr_timber_demand", linecolor="darkseagreen3", ls=:dot)
#-
plot(times, base.h[ix] , lab = "h - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Harvested and regeneration area (increased wood demand)");
plot!(times, base.r[ix] , lab = "r - base", linecolor="darkgreen", ls=:dot);
plot!(times, out.h[ix] , lab = "h - incr_timber_demand", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.r[ix] , lab = "r - incr_timber_demand", linecolor="darkseagreen3", ls=:dot)
# Shadow prices..
plot(times[2:end],  out.pF[ix[2:end]], lab = "pF: primary forest area price", linecolor="darkgreen", xlabel="years", title="Land shadow prices (increased wood demand)");
plot!(times[2:end], out.pS[ix[2:end]], lab = "pS: secondary forest area price", linecolor="darkseagreen3");
plot!(times[2:end], out.pA[ix[2:end]], lab = "pA: agricultural area price",linecolor="sienna");
plot!(times[2:end], out.pV[ix[2:end]], lab = "pV: secondary forest timber volumes price",linecolor="lightgreen")
#-


# **Take home**
# Increased timber demand would favour the switch from PF to SF
# Both base and `incr_timber_demand` leads to the same SF density equilibrium (this depends from the discount rate) but in base the density is gradually reduced, while in `incr_timber_demand` there is a first large harvesting, followed by small but  gradual 


# ------------------------------------------------------------------------------
# #### Scen 6: `lower_disc_rate`: Decreased discount rate
out = luc_model(σ=OLUA.σ-0.01)

plot(times, base.F[ix], lab = "F - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation under lower discount rate");
plot!(times, out.F[ix], lab = "F - lower_disc_rate", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - lower_disc_rate", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - lower_disc_rate", linecolor="sienna", ls=:dot)
#-
plot(times, base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times, out.V[ix] ./ out.S[ix],   lab = "lower_disc_rate", linecolor="darkseagreen3", ls=:dot)

# **Take home**
# Lower discount rate surprisily has a positive effect in reducing deforestation, primarily by reducing allocation to agricultural areas. Even with less avilable area from the less deforestation, SF area increases.
# Also, we note that the eq density of SF, as the discount rate decrease, move up toward the MSY. We find again a well know principle in capital standard capital theory and nat res economics: as harvesting SF doesn't depend from SF stocks, the intertemporal equilibrium is obtained when the rate of biological growth (dG/dV) equals the discount rate (PERMAN, Natural resource and environmental economics, 4th ed., p. 583) , and when this is zero this corresponds to the MSY, as we saw in the validation section.

# ------------------------------------------------------------------------------
# #### Scen 7: `cc_effect_pf`: Climate change effects (reduced PF density)
out_cc_effect_pf = luc_model(D = OLUA.D * 0.8)
out = out_cc_effect_pf
# CC effect under the hipothesis of PFbeing more involved

plot(times, base.F[ix], lab = "F - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith cc impact of PF");
plot!(times, out.F[ix], lab = "F - cc_effect_pf", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - cc_effect_pf", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - cc_effect_pf", linecolor="sienna", ls=:dot)
#-
plot(times, base.h[ix] .* base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF hV");
plot!(times, out.h[ix] .* out.V[ix] ./ out.S[ix],   lab = "cc_effect_pf", linecolor="darkseagreen3", ls=:dot)

# **Take home**
# Even thought the modelled effect concerns only the secondary forest, we see that deforestation of primary forest is also impacted, most likely this is a spillover effect due to the need to compensate the lower harvesting volumes coming from secondary forests. 

# ------------------------------------------------------------------------------
# #### Scen 8: `cc_effect_sf`: Climate change effects (reduced forest growth rate)
out_cc_effect_sf = luc_model(γ = OLUA.γ-0.02)
out = out_cc_effect_sf

# In our simple model we model CC effects on the forest (mortality and growth) as a reduction of the growth rate of the overall secondary forest aggregate density. We do not consider here the effects on the primary forest nor on the agricultural sector.

plot(times, base.F[ix], lab = "F - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith cc impact on SF");
plot!(times, out.F[ix], lab = "F - cc_effect_sf", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - cc_effect_sf", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - cc_effect_sf", linecolor="sienna", ls=:dot)
#-
plot(times, base.h[ix] .* base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF hV");
plot!(times, out.h[ix] .* out.V[ix] ./ out.S[ix],   lab = "cc_effect_sf", linecolor="darkseagreen3", ls=:dot)

# **Take home**
# Even thought the modelled effect concerns only the secondary forest, we see that deforestation of primary forest is also impacted, most likely this is a spillover effect due to the need to compensate the lower harvesting volumes coming from secondary forests. 

# ------------------------------------------------------------------------------
# #### Scen 8: `cc_effect_ag`: Climate change effects (reduced forest growth rate)
out_cc_effect_ag = luc_model(bagr_c1 = OLUA.bagr_c1*0.8)
out = out_cc_effect_sf

plot(times, base.F[ix], lab = "F - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith cc impact on ag");
plot!(times, out.F[ix], lab = "F - cc_effect_ag", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - cc_effect_ag", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - cc_effect_ag", linecolor="sienna", ls=:dot)
#-
plot(times, base.h[ix] .* base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF hV");
plot!(times, out.h[ix] .* out.V[ix] ./ out.S[ix],   lab = "cc_effect_ag", linecolor="darkseagreen3", ls=:dot)

# **Take home**
#


plot(times, base.h[ix] .* base.V[ix] ./ base.S[ix], lab = "base", linewidth = 2, linecolor="grey", xlabel="years", title="Harvesing volumes in SF under different scenarios")
plot!(times, out_cc_effect_pf.h[ix] .* out_cc_effect_pf.V[ix] ./ out_cc_effect_pf.S[ix], lab = "cc_effect_pf", linewidth = 2, linecolor="darkgreen")
plot!(times, out_cc_effect_sf.h[ix] .* out_cc_effect_sf.V[ix] ./ out_cc_effect_sf.S[ix], lab = "cc_effect_sf", linewidth = 2, linecolor="lightgreen")
plot!(times, out_cc_effect_ag.h[ix] .* out_cc_effect_ag.V[ix] ./ out_cc_effect_ag.S[ix], lab = "cc_effect_ag", linewidth = 2, linecolor="darkred")
savefig("hV_by_cc_scenario.png")


# ------------------------------------------------------------------------------
# Carbon sequestration vs carbon substitution


out_seq = luc_model(bc_seq_c1=100.0)
out_sub = luc_model(bc_sub_c1=100.0)
out_totc = luc_model(bc_seq_c1=100.0, bc_sub_c1=100.0)


timesi = times100
ixs    = ixs100

marginal_arrays =  [base.co2_seq[ixs], base.co2_sub[ixs], base.co2_seq[ixs] .+ base.co2_sub[ixs],
                    out_seq.co2_seq[ixs], out_seq.co2_sub[ixs], out_seq.co2_seq[ixs] .+ out_seq.co2_sub[ixs],
                    out_sub.co2_seq[ixs], out_sub.co2_sub[ixs], out_sub.co2_seq[ixs] .+ out_sub.co2_sub[ixs],
                    out_totc.co2_seq[ixs], out_totc.co2_sub[ixs], out_totc.co2_seq[ixs] .+ out_totc.co2_sub[ixs]]
cumulative_arrays = [cumsum(a .* step) for a in marginal_arrays]

marg_lims = (minimum(minimum.(marginal_arrays)), maximum(maximum.(marginal_arrays)))
cum_lims = (minimum(minimum.(cumulative_arrays)), maximum(maximum.(cumulative_arrays)))

# Yeary values...
p_base_m = plot(timesi, base.co2_seq[ixs], lab = "Sequestered carbon", linewidth = 1, linecolor="grey", title="No carbon evaluation", ls=:dot, ylabel="Mt CO₂ / y", ylims=marg_lims,xformatter=Returns(""), titlefontsize=12);
plot!(p_base_m, timesi, base.co2_sub[ixs], lab = "Substituted carbon", linewidth = 1, linecolor="grey", ls=:dash);
plot!(p_base_m, timesi, base.co2_seq[ixs] .+ base.co2_sub[ixs], lab = "Total", linewidth = 2, linecolor="grey", ls=:solid);

p_seq_m = plot(timesi, out_seq.co2_seq[ixs], lab = nothing, linewidth = 1, linecolor="darkgreen", ls=:dot, ylims=marg_lims,title="Seq evaluation",formatter=Returns(""), titlefontsize=12);
plot!(p_seq_m, timesi, out_seq.co2_sub[ixs], lab = nothing, linewidth = 1, linecolor="darkgreen", ls=:dash);
plot!(p_seq_m, timesi, out_seq.co2_seq[ixs] .+ out_seq.co2_sub[ixs], lab = nothing, linewidth = 2, linecolor="darkgreen", ls=:solid);

p_sub_m = plot(timesi, out_sub.co2_seq[ixs], lab = nothing, linewidth = 1, linecolor="blue", ls=:dot, ylims=marg_lims,title="Sub evaluation",formatter=Returns(""), titlefontsize=12);
plot!(p_sub_m, timesi, out_sub.co2_sub[ixs], lab = nothing, linewidth = 1, linecolor="blue", ls=:dash);
plot!(p_sub_m, timesi, out_sub.co2_seq[ixs] .+ out_sub.co2_sub[ixs], lab = nothing, linewidth = 2, linecolor="blue", ls=:solid);

p_totc_m = plot(timesi, out_totc.co2_seq[ixs], lab = nothing, linewidth = 1, linecolor="sienna", ls=:dot, ylims=marg_lims,title="Seq+sub evaluation",formatter=Returns(""), titlefontsize=12);
plot!(p_totc_m, timesi, out_totc.co2_sub[ixs], lab = nothing, linewidth = 1, linecolor="sienna", ls=:dash);
plot!(p_totc_m, timesi, out_totc.co2_seq[ixs] .+ out_totc.co2_sub[ixs], lab = nothing, linewidth = 2, linecolor="sienna", ls=:solid);

# Cumulative values...
p_base_c = plot(timesi, cumsum(base.co2_seq[ixs] .* step), lab = nothing, linewidth = 1, linecolor="grey", ls=:dot, ylabel="Cum Mt CO₂", ylims=cum_lims);
plot!(p_base_c, timesi, cumsum(base.co2_sub[ixs] .* step), lab = nothing, linewidth = 1, linecolor="grey", ls=:dash);
plot!(p_base_c, timesi, cumsum((base.co2_seq[ixs] .+ base.co2_sub[ixs]) .* step), lab = nothing, linewidth = 2, linecolor="grey", ls=:solid);

p_seq_c = plot(timesi, cumsum(out_seq.co2_seq[ixs] .* step), lab = nothing, linewidth = 1, linecolor="darkgreen", xlabel="years", ls=:dot, ylims=cum_lims,yformatter=Returns(""));
plot!(p_seq_c, timesi, cumsum(out_seq.co2_sub[ixs] .* step), lab = nothing, linewidth = 1, linecolor="darkgreen", ls=:dash, ylims=cum_lims);
plot!(p_seq_c, timesi, cumsum((out_seq.co2_seq[ixs] .+ out_seq.co2_sub[ixs]) .* step), lab = nothing, linewidth = 2, linecolor="darkgreen", ls=:solid, ylims=cum_lims);

p_sub_c = plot(timesi, cumsum(out_sub.co2_seq[ixs] .* step), lab = nothing, linewidth = 1, linecolor="blue", ls=:dot, ylims=cum_lims,yformatter=Returns(""));
plot!(p_sub_c, timesi, cumsum(out_sub.co2_sub[ixs] .* step), lab = nothing, linewidth = 1, linecolor="blue", ls=:dash);
plot!(p_sub_c, timesi, cumsum((out_sub.co2_seq[ixs] .+ out_sub.co2_sub[ixs]) .* step), lab = nothing, linewidth = 2, linecolor="blue", ls=:solid);

p_totc_c = plot(timesi, cumsum(out_totc.co2_seq[ixs] .* step), lab = nothing, linewidth = 1, linecolor="sienna", ls=:dot, ylims=cum_lims,yformatter=Returns(""));
plot!(p_totc_c, timesi, cumsum(out_totc.co2_sub[ixs] .* step), lab = nothing, linewidth = 1, linecolor="sienna", ls=:dash);
plot!(p_totc_c, timesi, cumsum((out_totc.co2_seq[ixs] .+ out_totc.co2_sub[ixs]) .* step), lab = nothing, linewidth = 2, linecolor="sienna", ls=:solid);

plot(p_base_m,p_seq_m,p_sub_m,p_totc_m,p_base_c,p_seq_c,p_sub_c,p_totc_c, layout = (2,4),size=(1200,600),left_margin = [20px -10px], bottom_margin = [30px -30px])

savefig("carbon_out.svg")

# The same, but now only energy substitution:


out_seq = luc_model(bc_seq_c1=100.0,co2sub=0.264)
out_sub = luc_model(bc_sub_c1=100.0,co2sub=0.264)
out_totc = luc_model(bc_seq_c1=100.0, bc_sub_c1=100.0),co2sub=0.264


timesi = times100
ixs    = ixs100

marginal_arrays =  [base.co2_seq[ixs], base.co2_sub[ixs], base.co2_seq[ixs] .+ base.co2_sub[ixs],
                    out_seq.co2_seq[ixs], out_seq.co2_sub[ixs], out_seq.co2_seq[ixs] .+ out_seq.co2_sub[ixs],
                    out_sub.co2_seq[ixs], out_sub.co2_sub[ixs], out_sub.co2_seq[ixs] .+ out_sub.co2_sub[ixs],
                    out_totc.co2_seq[ixs], out_totc.co2_sub[ixs], out_totc.co2_seq[ixs] .+ out_totc.co2_sub[ixs]]
cumulative_arrays = [cumsum(a .* step) for a in marginal_arrays]

marg_lims = (minimum(minimum.(marginal_arrays)), maximum(maximum.(marginal_arrays)))
cum_lims = (minimum(minimum.(cumulative_arrays)), maximum(maximum.(cumulative_arrays)))

# Yeary values...
p_base_m = plot(timesi, base.co2_seq[ixs], lab = "Sequestered carbon", linewidth = 1, linecolor="grey", title="No carbon evaluation", ls=:dot, ylabel="Mt CO₂ / y", ylims=marg_lims,xformatter=Returns(""), titlefontsize=12);
plot!(p_base_m, timesi, base.co2_sub[ixs], lab = "Substituted carbon", linewidth = 1, linecolor="grey", ls=:dash);
plot!(p_base_m, timesi, base.co2_seq[ixs] .+ base.co2_sub[ixs], lab = "Total", linewidth = 2, linecolor="grey", ls=:solid);

p_seq_m = plot(timesi, out_seq.co2_seq[ixs], lab = nothing, linewidth = 1, linecolor="darkgreen", ls=:dot, ylims=marg_lims,title="Seq evaluation",formatter=Returns(""), titlefontsize=12);
plot!(p_seq_m, timesi, out_seq.co2_sub[ixs], lab = nothing, linewidth = 1, linecolor="darkgreen", ls=:dash);
plot!(p_seq_m, timesi, out_seq.co2_seq[ixs] .+ out_seq.co2_sub[ixs], lab = nothing, linewidth = 2, linecolor="darkgreen", ls=:solid);

p_sub_m = plot(timesi, out_sub.co2_seq[ixs], lab = nothing, linewidth = 1, linecolor="blue", ls=:dot, ylims=marg_lims,title="Sub evaluation",formatter=Returns(""), titlefontsize=12);
plot!(p_sub_m, timesi, out_sub.co2_sub[ixs], lab = nothing, linewidth = 1, linecolor="blue", ls=:dash);
plot!(p_sub_m, timesi, out_sub.co2_seq[ixs] .+ out_sub.co2_sub[ixs], lab = nothing, linewidth = 2, linecolor="blue", ls=:solid);

p_totc_m = plot(timesi, out_totc.co2_seq[ixs], lab = nothing, linewidth = 1, linecolor="sienna", ls=:dot, ylims=marg_lims,title="Seq+sub evaluation",formatter=Returns(""), titlefontsize=12);
plot!(p_totc_m, timesi, out_totc.co2_sub[ixs], lab = nothing, linewidth = 1, linecolor="sienna", ls=:dash);
plot!(p_totc_m, timesi, out_totc.co2_seq[ixs] .+ out_totc.co2_sub[ixs], lab = nothing, linewidth = 2, linecolor="sienna", ls=:solid);

# Cumulative values...
p_base_c = plot(timesi, cumsum(base.co2_seq[ixs] .* step), lab = nothing, linewidth = 1, linecolor="grey", ls=:dot, ylabel="Cum Mt CO₂", ylims=cum_lims);
plot!(p_base_c, timesi, cumsum(base.co2_sub[ixs] .* step), lab = nothing, linewidth = 1, linecolor="grey", ls=:dash);
plot!(p_base_c, timesi, cumsum((base.co2_seq[ixs] .+ base.co2_sub[ixs]) .* step), lab = nothing, linewidth = 2, linecolor="grey", ls=:solid);

p_seq_c = plot(timesi, cumsum(out_seq.co2_seq[ixs] .* step), lab = nothing, linewidth = 1, linecolor="darkgreen", xlabel="years", ls=:dot, ylims=cum_lims,yformatter=Returns(""));
plot!(p_seq_c, timesi, cumsum(out_seq.co2_sub[ixs] .* step), lab = nothing, linewidth = 1, linecolor="darkgreen", ls=:dash, ylims=cum_lims);
plot!(p_seq_c, timesi, cumsum((out_seq.co2_seq[ixs] .+ out_seq.co2_sub[ixs]) .* step), lab = nothing, linewidth = 2, linecolor="darkgreen", ls=:solid, ylims=cum_lims);

p_sub_c = plot(timesi, cumsum(out_sub.co2_seq[ixs] .* step), lab = nothing, linewidth = 1, linecolor="blue", ls=:dot, ylims=cum_lims,yformatter=Returns(""));
plot!(p_sub_c, timesi, cumsum(out_sub.co2_sub[ixs] .* step), lab = nothing, linewidth = 1, linecolor="blue", ls=:dash);
plot!(p_sub_c, timesi, cumsum((out_sub.co2_seq[ixs] .+ out_sub.co2_sub[ixs]) .* step), lab = nothing, linewidth = 2, linecolor="blue", ls=:solid);

p_totc_c = plot(timesi, cumsum(out_totc.co2_seq[ixs] .* step), lab = nothing, linewidth = 1, linecolor="sienna", ls=:dot, ylims=cum_lims,yformatter=Returns(""));
plot!(p_totc_c, timesi, cumsum(out_totc.co2_sub[ixs] .* step), lab = nothing, linewidth = 1, linecolor="sienna", ls=:dash);
plot!(p_totc_c, timesi, cumsum((out_totc.co2_seq[ixs] .+ out_totc.co2_sub[ixs]) .* step), lab = nothing, linewidth = 2, linecolor="sienna", ls=:solid);

plot(p_base_m,p_seq_m,p_sub_m,p_totc_m,p_base_c,p_seq_c,p_sub_c,p_totc_c, layout = (2,4),size=(1200,600),left_margin = [20px -10px], bottom_margin = [30px -30px])

savefig("carbon_out_sub_onlyenergy.svg")


# ------------------------------------------------------------------------------

# Verification of the steady state

# verification of h equilibrium
rhs = @. OLUA.γ*(base.S * OLUA.K - base.V) / OLUA.K
lhs = base.h

plot(times, rhs[ix])
plot!(times,base.h[ix])

isapprox(base.h[findfirst(t-> t==200,times)], rhs[findfirst(t-> t==200,times)], rtol=0.001)

base.h[findfirst(t-> t==200,times)]

OLUA.crsf_c1

base.pS
base.pA

pricediff = base.pA .- base.pS

plot(times, base.pS[ix])
plot!(times, base.pA[ix])
plot(times, pricediff[ix])