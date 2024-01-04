# Forest dynamics numerical model

cd(@__DIR__)
using Pkg
Pkg.activate(".")
#Pkg.resolve()
#Pkg.Instantiate()

using Plots
using Markdown
using Revise
using Test

# Load the model structure (function "luc_model")
includet("model.jl")
using .OLUA

# ------------------------------------------------------------------------------
# Symbols
md"""

Abbreviations:
- `PF`: primary forests
- `SF`: secondary forests
- `AG`: agricultural area

## Assumptions
- h1: Primary forests are in stationary state, i.e. their volume variation is only given by clear cutting and its density `D` is fixed
- h2: Harvested area of $PF$ can be allocated to either $SF$ or $AG$
- h3: Harvesting of $SF$ doesn't change its area, i.e. if clear-cut, the area remains $SF$
- h4: At each moment in time, benefits for the society depend from area of $PF$ (environmental benefits), sum of harvesting volumes of $PF$ and $SF$ (i.e. indifferentiated wood), area of $AG$
- h5: Costs of harvesting $PF$ are inversely proportional to the area of $PF$
- h6: The society can "control" the harvesting of primary forest $d$, the harvesting of secondary forest $h$ and the allocation of harvested primary forest
"""

# Loading the "base" optimisation...
base     = luc_model()
ix       = 1:(length(base.support)-300) # index for plotting
times    = base.support[ix] # Every 5 years
ntpoints = length(times)
step     = times[end] / (ntpoints-1)

# ------------------------------------------------------------------------------
# Test 1: using less or more points doesn't influence much the results
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

plot(times, base.F[1:101], xlabel="years", label="Base time point density (5 y)", title="Forest prim area (stock) under different time densities");
plot!(times_dense, out_dense.F[1:251], label="Dense time point density (2 y)");
plot!(times_sparce, out_sparce.F[1:51], label="Sparce time point density (10 y)")

# ------------------------------------------------------------------------------
# Test 2: setting no interest rate leads to MSY in secondary forest (density = half the maximum density)

out = luc_model(σ=0)

Dsf = out.V ./ out.S # Secondary forest density
@testset "No discount leads to MSY in SF" begin
  @test isapprox(Dsf[findfirst(t-> t==100,times)], OLUA.K/2, rtol=0.01)
  @test isapprox(Dsf[findfirst(t-> t==200,times)], OLUA.K/2, rtol=0.01)
end

# Graphically...
plot(times, out.V[ix] ./ out.S[ix], lab = "Secondary forest density", xlabel="years", linecolor="darkseagreen3", title= "Secondary forest density under no discount");
plot!(times, fill(OLUA.K,ntpoints), lab="Max density")

# ------------------------------------------------------------------------------
# Test 3: setting leads to the "eating the cake" model

# The first version is with only benefits, the second one with some harvesting costs
# The second version has some numerical instability at the beginning
out1 = luc_model(benv_c1=0.0,bagr_c1=0,chpf_c3=0,chpf_c1=0,crsf_c1=1000,chsf_c1=10000,γ=0) # only benefits
out2 = luc_model(benv_c1=0.0,bagr_c1=0,chpf_c3=0,chpf_c1=10,chpf_c2=1.2,crsf_c1=1000,chsf_c1=10000,γ=0,d₀=1.13) # with costs

dbw_dV(V,c1=OLUA.bwood_c1,c2=OLUA.bwood_c2) = c1*c2*V^(c2-1) # marginal benetits
dcw_dV(V,c1=10,c2=1.2) = c1*c2*V^(c2-1) 

hv_1 = @. out1.d * OLUA.D + out1.h * out1.V / out1.S
hv_2 = @. out2.d * OLUA.D + out2.h * out2.V / out2.S

np1 = dbw_dV.(hv_1) # net price
np2 = dbw_dV.(hv_2) .- dcw_dV.(hv_2) # net price

r1a = log(np1[28]/np1[27])/step # growth rate of the price
r1b = log(np1[25]/np1[20])/(step*5) 
r2a = log(np2[28]/np2[27])/step # growth rate of the price
r2b = log(np2[25]/np2[20])/(step*5) 

@testset "Checking net price growth following Hotelling rule" begin
  @test isapprox(r1a,OLUA.σ,rtol=0.01)
  @test isapprox(r1b,OLUA.σ,rtol=0.01)
  @test isapprox(r2a,OLUA.σ,rtol=0.05)
  @test isapprox(r2b,OLUA.σ,rtol=0.05)
end

# Graphically...
plot(times, np1[ix], lab = "Considering only benefits", xlabel="years", title= "Net price of the timber resource")
plot!(times, np2[ix], lab = "Considering harvesting costs")

# ------------------------------------------------------------------------------
# Test 4: setting no Sf area change and no harvesting leads SF volumes to the Verhulst model

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

plot(times,out.V[ix], label="V: model output logistic", xlabel="years", title= "Sec. forest volumes under no harv and no reg");

# manually computing the V increase by loop
var_vol_test(S,V;γ=OLUA.γ,K=OLUA.K) = (V)*γ*(1-(V / (S * K) ))
ts2 = 0:step:times[end]
Vtest = copy(OLUA.V₀)
v_by_step = zeros(length(ts2))
for (it,t) in enumerate(ts2)
  v_by_step[it] = Vtest
  Vtest += var_vol_test(OLUA.S₀,Vtest)*step
end
Vtest
plot!(ts2, collect(v_by_step[i] for i in 1:length(ts2)), label="V: discrete steps logistic" );

# using the true logistic function...
logistic(x;k,r,x0) = k/(1+((k-x0)/x0)*exp(-r*x))
plot!(x->logistic(x,k=OLUA.K*OLUA.S₀,r=OLUA.γ,x0=OLUA.V₀),0,times[end], label = "V: true logistic funcion")


# ------------------------------------------------------------------------------
# Base Case Analysis









plot(times, out.S[ix], lab = "S: secondary forest area", linecolor="darkseagreen3")
plot(times, out.V[ix], lab = "V: secondary forest volumes", linecolor="darkseagreen3")
plot(times, Dsf[ix],   lab = "D: secondary forest density", linecolor="darkseagreen3")
plot(times, out.h[ix], lab = "h", linecolor="darkseagreen3")
plot(times, out.r[ix], lab = "r", linecolor="darkseagreen3")




# Test not passed.. this is strange, check why the hell it goes just a bit over K/2 instead of going to full carrying capacity
# It isn't influenced by the time density


plot(times,  out.F[ix], lab = "F: primary forest area", linecolor="darkgreen", title="Land Areas")
plot!(times, out.S[ix], lab = "S: secondary forest area", linecolor="darkseagreen3")
plot!(times, out.A[ix], lab = "A: agricultural area",linecolor="sienna")

plot(times,  out.d[ix], lab = "d: primary forest harvested area", linecolor="darkgreen", title="Land Variations")
plot!(times, out.r[ix], lab = "r: secondary forest regenerated area", linecolor="darkseagreen3")
plot!(times, out.d[ix] - out.r[ix], lab = "new agricultural area", linecolor="sienna")


plot(times, out.F[ix] .* OLUA.D, lab = "V (pf): primary forest volumes",linecolor="darkgreen", title="Growing Volumes")
plot!(times, out.V[ix], lab = "V: secondary forest volumes", linecolor="darkseagreen3")
# 20230913: Checked by setting dA_sf = 0 and dV_sf only to the natural dynamics (removing harvesting) that the path of V_opt_sf follows the logistic function with the given parameters 
plot(times, out.F[ix] .* OLUA.D .+ out.V[ix],  lab = "TOT V: total forest volumes", linecolor="darkseagreen3")




plot(times,  out.d[ix], lab = "d: primary forest harvested area", linecolor="darkgreen", title="Distribution of Harvested areas")
plot!(times, out.h[ix], lab = "h: secondary forest harvested area",linecolor="darkseagreen3")



plot(times,  out.d[ix] .* OLUA.D, lab = "hV_pf: primary forest harvested volumes", linecolor="darkgreen", title="Distribution of Harvested volumes")
plot!(times, out.h[ix] .* out.V[ix] ./ out.S[ix], lab = "hV_sf: secondary forest harvested volumes", linecolor="darkseagreen3")

plot!(times, out.d[ix] .* OLUA.D .+ out.h[ix] .* out.V[ix] ./ out.S[ix], lab = "hV: total forest harvested volumes", linecolor="darkseagreen3")

plot(times, out.V[ix] ./ out.S[ix], lab = "Secondary forest density", linecolor="darkseagreen3", title= "Secondary forest density")


plot(times,  out.ben_env[ix], lab = "Environmental benefits", linecolor="darkgreen", title="Benefits") 
plot!(times, out.ben_agr[ix], lab = "Agr benefits",linecolor="sienna") 
plot!(times, out.ben_wood[ix], lab = "Wood use benefits", linecolor="darkseagreen3")
plot!(times, out.ben_carbon[ix], lab = "Carbon benefits", linecolor="grey")
plot!(times, .- out.cost_pfharv[ix], lab = "PF harvesting costs", linecolor="darkgreen")
plot!(times, .- out.cost_sfharv[ix], lab = "SF harvesting costs", linecolor="darkseagreen3")
plot!(times, .- out.cost_sfreg[ix], lab = "SF regeneration costs",linecolor="sienna")

plot(times, out.welfare[ix], lab = "Total welfare") 

out_msy = luc_model(σ  = 0)
plot(times,  out.V[ix] ./ out.S[ix], lab = "BAU", linecolor="darkseagreen3", title="D: secondary forest density")
plot!(times, out_msy.V[ix] ./ out_msy.S[ix], lab = "Increased env benefits", linecolor="darkseagreen3", ls=:dot)



# Increase of Environmental Benefits

out_benv = luc_model(benv_c1 = OLUA.benv_c1*1.5)
plot(times,  out.F[ix], lab = "BAU", linecolor="darkgreen", title="F: primary forest area")
plot!(times,  out_benv.F[ix], lab = "Increased env benefits", linecolor="darkgreen",ls=:dot )

plot(times, out.S[ix], lab = "BAU", linecolor="darkseagreen3", title="S: secondary forest area")
plot!(times, out_benv.S[ix], lab = "Increased env benefits", linecolor="darkseagreen3", ls=:dot)

plot(times, out.A[ix], lab = "BAU", linecolor="sienna", title="A: agricultural area")
plot!(times, out_benv.A[ix], lab = "Increased env benefits", linecolor="sienna", ls=:dot)


plot(times, out.h[ix] .* out.V[ix] ./ out.S[ix], lab = "BAU", linecolor="darkseagreen3", title="hV_sf: secondary forest harvested volumes")
plot!(times, out_benv.h[ix] .* out_benv.V[ix] ./ out_benv.S[ix], lab = "Increased env benefits", linecolor="darkseagreen3", ls=:dot)

plot(times,  out.V[ix] ./ out.S[ix], lab = "BAU", linecolor="darkseagreen3", title="D: secondary forest density")
plot!(times, out_benv.V[ix] ./ out_benv.S[ix], lab = "Increased env benefits", linecolor="darkseagreen3", ls=:dot)



plot(times,  out.F[ix], lab = "F: primary forest area")
plot!(times,  out_fb.F[ix], lab = "F: primary forest area (increased forest benefits)")

plot(times, out.S[ix], lab = "S: secondary forest area")
plot!(times, out_fb.S[ix], lab = "S: secondary forest area (increased forest benefits)")

plot(times, out.A[ix], lab = "A: agricultural area")
plot!(times, out_fb.A[ix], lab = "A: agricultural area (increased forest benefits)")





ben_env_temp(F;benv_c1=benv_c1,benv_c2=benv_c2) = benv_c1*F^benv_c2 

σt = 0.03

ix = 1:(length(out.support)-180) # index for plotting
times = out.support[ix]
sum( ben_env_temp(out.F[i])*exp(-times[i]*σt ) for i in ix )
sum( ben_env_temp(out.F[1])*exp(-times[i]*σt ) for i in ix )

[out.F[i] >= out.F[1] for i in ix]

