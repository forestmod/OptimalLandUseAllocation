# # Forest dynamics numerical model

# Literate.markdown("harvest_dynamics.jl", "."; flavor=Literate.CommonMarkFlavor(), execute=true) #src

# ### Setting up the environment / packages...
cd(@__DIR__)
using Pkg
Pkg.activate(".")
#Pkg.resolve()
#Pkg.Instantiate()

using Plots
using Markdown
using Revise
using Test

# ### Load the model structure (function "luc_model")
include("model.jl")
using .OLUA


# ------------------------------------------------------------------------------
# ### Compute the "base" optimisation...
base     = luc_model() # compute the optimization for the "base" scenario
ix       = 1:(length(base.support)-300) # index for plotting (we discard the distant future)
times    = base.support[ix] 
ntpoints = length(times)
step     = times[end] / (ntpoints-1) ; # We choosen to discretize every 5 years

# ------------------------------------------------------------------------------
# ### Model validation...
# As the overall model is analytically too complex, we opted by a validation "by parts", that is, we set the parameters of the models such that some components are disables (for example setting wood benefits to zero disable harvesting) so that the remaining parts have predictable characteristics.
# We performed 4 validation exercises: in the first one we test the effect of the specific time discretization emploied; in the second test we verify that with 0 discount rate, the SF is harvested such to reach in the long run the MSY; in the third test we verify that when we have benefits only from harvesting the PF we end up to a "eating the cake" model, and in particular the marginal penefit of wood (net price) follows the Hotelling's rule. Finally in the last validation test we check that without harvesting the SF follows the Verhulst (logistic) model.
# Each validation is performed with a in-code test that the results follow the assumptions and by plotting a chart of the relevant outputs.  

# ------------------------------------------------------------------------------
# #### Test 1: discretization choice (using less or more time points) doesn't influence much the results
out_dense    = luc_model(ns=1001,opt_options = Dict("max_cpu_time" => 60.0))
ix_dense     = 1:(length(out_dense.support)-750) # index for plotting
times_dense  = out_dense.support[ix_dense] # every 2 years
out_sparce   = luc_model(ns=201)
ix_sparce    = 1:(length(out_sparce.support)-150) # index for plotting
times_sparce = out_sparce.support[ix_sparce] # every 10 years
@testset "Number of discretization points" begin
    @test isapprox(base.r[findfirst(t-> t==80,times)],out_dense.r[findfirst(t-> t==80,times_dense)],rtol=0.1)
    @test isapprox(base.r[findfirst(t-> t==80,times)],out_sparce.r[findfirst(t-> t==80,times_sparce)],rtol=0.15)
    @test isapprox(base.F[findfirst(t-> t==80,times)],out_dense.F[findfirst(t-> t==80,times_dense)],rtol=0.05)
    @test isapprox(base.F[findfirst(t-> t==80,times)],out_sparce.F[findfirst(t-> t==80,times_sparce)],rtol=0.05)
end

# Graphically...
plot(times, base.r[1:101], xlabel="years", label="Base time point density (5 y)", title="SF reg area (flow) under different time densities");
plot!(times_dense, out_dense.r[1:251], label="Dense time point density (2 y)");
plot!(times_sparce, out_sparce.r[1:51], label="Sparce time point density (10 y)")
#-
plot(times, base.F[1:101], xlabel="years", label="Base time point density (5 y)", title="Forest prim area (stock) under different time densities");
plot!(times_dense, out_dense.F[1:251], label="Dense time point density (2 y)");
plot!(times_sparce, out_sparce.F[1:51], label="Sparce time point density (10 y)")

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
out2 = luc_model(benv_c1=0.0,bagr_c1=0,chpf_c3=0,chpf_c1=10,chpf_c2=1.2,crsf_c1=1000,chsf_c1=10000,γ=0,d₀=1.13) # with timber benefits and costs

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

out = luc_model(bwood_c1=0.01,h₀=0)
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
# ### Base Case Analysis

# Areas..
plot(times,  base.F[ix], lab = "F: primary forest area", linecolor="darkgreen", xlabel="years", title="Land Areas (base scen)");
plot!(times, base.S[ix], lab = "S: secondary forest area", linecolor="darkseagreen3");
plot!(times, base.A[ix], lab = "A: agricultural area",linecolor="sienna")

# Volumes...
plot(times, base.F[ix] .* OLUA.D, lab = "V (pf): primary forest volumes",linecolor="darkgreen", title="Growing Volumes (base scen)", xlabel="title");
plot!(times, base.V[ix], lab = "V (sf): secondary forest volumes", linecolor="darkseagreen3");
plot!(times, base.F[ix] .* OLUA.D .+ base.V[ix],  lab = "TOT V: total forest volumes", linecolor="brown")

# Welfare analysis...
plot(times,  base.ben_env[ix], lab = "Environmental benefits", linecolor="darkgreen", title="Welfare balance (base scen)"); 
plot!(times, base.ben_agr[ix], lab = "Agr benefits",linecolor="sienna") ;
plot!(times, base.ben_wood[ix], lab = "Wood use benefits", linecolor="darkseagreen3");
plot!(times, base.ben_carbon[ix], lab = "Carbon benefits", linecolor="grey");
plot!(times, .- base.cost_pfharv[ix], lab = "PF harvesting costs", linecolor="darkgreen");
plot!(times, .- base.cost_sfharv[ix], lab = "SF harvesting costs", linecolor="darkseagreen3");
plot!(times, .- base.cost_sfreg[ix], lab = "SF regeneration costs",linecolor="sienna")
#-
plot(times[1:101], base.welfare[1:101], lab = "Total welfare (base scen)", ylims=(0,1E6) )

# ------------------------------------------------------------------------------
# ### Scenario analysis

# ------------------------------------------------------------------------------
# #### Scen 1: `no_env_ben`: PF environmental benefits not considered

out = luc_model(benv_c1=0.0)

plot(times, base.F[ix], lab = "F - base", linecolor="darkgreen", linewidth = 2, xlabel="years", title="Land allocation under no env benefits");
plot!(times, out.F[ix], lab = "F - no_env_ben", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - no_env_ben", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - no_env_ben", linecolor="sienna", ls=:dot)

# **Take home**
# Non considering environmental benefits would, as expect, largelly faster the deforestation of primary forests. Interesting, the increased freed areas would largelly be allocated to agriculture, as the deforestation would imply larger supply of timber that would not benefits the secondary forests (our timber benefits function is concave, and the timber from primary and secondary forest is a completelly homogeneous product)

# ------------------------------------------------------------------------------
# #### Scen 2: Scen `with_carb_ben_1``: Carbon benefits (storage) accounted for (constant carbon price)
out = luc_model(bc_c1=100.0)

plot(times, base.F[ix], lab = "F - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith fixed carb benefits");
plot!(times, out.F[ix], lab = "F - with_carb_ben_1", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - with_carb_ben_1", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - with_carb_ben_1", linecolor="sienna", ls=:dot)
#-
plot(times, base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times, out.V[ix] ./ out.S[ix],   lab = "with_carb_ben_1", linecolor="darkseagreen3", ls=:dot)

# **Take home**
# Interesting, because the maximum density of secondary forest is (not much) higher than primary forests, considering carbon value for timber sequestration would accellerate deforestation in favour of secondary forests
# In other words, as feared by some, carbon payments for forest sequestration could indeed favour monospecific plantations

# ------------------------------------------------------------------------------
# #### Scen 3: `with_carb_ben_2`: Carbon benefits (storage) accounted for (with growing carbon price)
out = luc_model(bc_c1=100.0,bc_c2=(OLUA.σ-0.005)) # setting bc_c2=OLUA.σ doesn't solve

plot(times, base.F[ix], lab = "F - base", linewidth = 2, linecolor="darkgreen", xlabel="years", title="Land allocation whith increasing carb benefits");
plot!(times, out.F[ix], lab = "F - with_carb_ben_2", linecolor="darkgreen", ls=:dot);
plot!(times, base.S[ix], lab = "S - base", linewidth = 2, linecolor="darkseagreen3");
plot!(times, out.S[ix], lab = "S - with_carb_ben_2", linecolor="darkseagreen3", ls=:dot);
plot!(times, base.A[ix], lab = "A - base", linewidth = 2, linecolor="sienna");
plot!(times, out.A[ix], lab = "A - with_carb_ben_2", linecolor="sienna", ls=:dot)
#-
plot(times, base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times, out.V[ix] ./ out.S[ix],   lab = "with_carb_ben_2", linecolor="darkseagreen3", ls=:dot)

# **Take home**
# Not what I did expected. I did expect that, becasue when you have to repay it is more expensive, it doesn't become appropriate to storage carbon.
# Instead it looks like it is just an increase version of the `with_carb_ben_2` scenario

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
plot(times, base.V[ix]./ base.S[ix],   lab = "base", linewidth = 2, linecolor="darkseagreen3", title="SF density");
plot!(times, out.V[ix] ./ out.S[ix],   lab = "incr_timber_demand", linecolor="darkseagreen3", ls=:dot)

# **Take home**
# Increased timber demand would favour the switch from PF to SF

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
# #### Scen 3: `lower_disc_rate`: Decreasaed discount rate
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

