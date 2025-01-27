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
push!(LOAD_PATH,joinpath(@__DIR__,"lib"))
using OLUA

bc_seq_c1 = 100.0
bc_sub_c1 = 100.0
σ = OLUA.σ


# ### Compute some indexes..
base     = luc_model()
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

# ------------------------------------------------------------------------------------------------------
# LAND USE CHANGE ONLY FROM PF TO SF or A
fvars = Dict("a"=>0.0,"r_A"=>0.0)

# Other things to try:
# - concerning fvars:
# fvars = Dict("a"=>0.0,"r_F"=>0.0,"r_A"=>0.0,"d"=>0.0,"h"=>0.0) # no no land use change and no harvest
# fvars = Dict("a"=>0.0,"r_F"=>0.0,"r_A"=>0.0,"d"=>0.0) # no land use change
# fvars=Dict() # all land use changes possible
# - concerning bc_seq_c1 and bc_sub_c1:
# bc_seq_c1 = 0.0; bc_sub_c1 = 0.0 # No carbon benefits


# ### Base Case Analysis
base     = luc_model(bc_seq_c1=bc_seq_c1 , bc_sub_c1 =bc_sub_c1, fvars=fvars) # compute the optimization for the scenario without damage


# Damage scenario
damage_rate = 0.5
tdamage     = 40 # for simplicity, always rounded to the nearest 5 year

fully_anticipated  = luc_model(bc_seq_c1=bc_seq_c1, bc_sub_c1 =bc_sub_c1, damage_rate=damage_rate, tdamage=tdamage, fvars=fvars) # compute the optimization for the scenario without damage

# Computing new state var at time of domage...
ix_domage = findfirst(t -> t >= tdamage, times)+1

F₀          = base.F[ix_domage]  # Initial primary-forest area
S₀          = base.S[ix_domage]  # Initial secondary forest area
A₀          = base.A[ix_domage]  # Initial agricultural area
V₀          = base.V[ix_domage]* damage_rate  # Initial secondary forest volumes

adaptation = luc_model(bc_seq_c1=bc_seq_c1, bc_sub_c1 =bc_sub_c1, F₀=F₀, S₀=S₀, A₀=A₀, V₀=V₀, fvars=fvars) 

range_part1 = 1:ix_domage-1
range_part2 = 1:length(base.support)-ix_domage+1

merged = (
  F = vcat(base.F[range_part1],adaptation.F[range_part2]),
  S = vcat(base.S[range_part1],adaptation.S[range_part2]),
  A = vcat(base.A[range_part1],adaptation.A[range_part2]),
  V = vcat(base.V[range_part1],adaptation.V[range_part2]),
  d = vcat(base.d[range_part1],adaptation.d[range_part2]),
  h = vcat(base.h[range_part1],adaptation.h[range_part2]),
  r_F = vcat(base.r_F[range_part1],adaptation.r_F[range_part2]),
  r_A = vcat(base.r_A[range_part1],adaptation.r_A[range_part2]),
  a   = vcat(base.a[range_part1],adaptation.a[range_part2]),
  ben_env = vcat(base.ben_env[range_part1],adaptation.ben_env[range_part2]),
  ben_agr = vcat(base.ben_agr[range_part1],adaptation.ben_agr[range_part2]),
  ben_wood = vcat(base.ben_wood[range_part1],adaptation.ben_wood[range_part2]),
  ben_carbon_seq = vcat(base.ben_carbon_seq[range_part1],adaptation.ben_carbon_seq[range_part2]),
  ben_carbon_sub = vcat(base.ben_carbon_sub[range_part1],adaptation.ben_carbon_sub[range_part2]),
  cost_pfharv = vcat(base.cost_pfharv[range_part1],adaptation.cost_pfharv[range_part2]),
  cost_sfharv = vcat(base.cost_sfharv[range_part1],adaptation.cost_sfharv[range_part2]),
  cost_sfreg = vcat(base.cost_sfreg[range_part1],adaptation.cost_sfreg[range_part2]),
  welfare = vcat(base.welfare[range_part1],adaptation.welfare[range_part2]),
  co2_seq = vcat(base.co2_seq[range_part1],adaptation.co2_seq[range_part2]),
  co2_sub = vcat(base.co2_sub[range_part1],adaptation.co2_sub[range_part2]),
  pF= vcat(base.pF[range_part1],adaptation.pF[range_part2]),
  pS= vcat(base.pS[range_part1],adaptation.pS[range_part2]),
  pA= vcat(base.pA[range_part1],adaptation.pA[range_part2]),
  pV= vcat(base.pV[range_part1],adaptation.pV[range_part2]),
)


# As in merged the damage happens out of the otimization, we need to subtract from welfare, ben_carbon_seq and co2_seq the values that are left because of the storm
v_damaged = merged.V[ix_domage-1]*damage_rate
co2_released = (v_damaged * OLUA.co2seq) / step # this because then wer are going to multiply each point value by the step
bc_seq_to_remove = bc_seq_c1*exp(OLUA.bc_seq_c2 * times[ix_domage-1]) * co2_released 
merged.co2_seq[ix_domage] -= co2_released
merged.ben_carbon_seq[ix_domage] -= bc_seq_to_remove
merged.welfare[ix_domage] -= bc_seq_to_remove


plot(times,base.welfare[ix],title="Welfare",label="base");
plot!(times,fully_anticipated.welfare[ix], label="fully anticipated");
plot!(times,merged.welfare[ix], label="ex-post adapt")

# Computation of total welfare...
full_ix      = 2:length(base.support)
w_base = 0.0
w_fully_anticipated = 0.0
w_merged = 0.0
for ix in full_ix
  t = base.support[ix]
  w_base += step*base.welfare[ix]*(1+σ)^-t
  w_fully_anticipated += step*fully_anticipated.welfare[ix]*(1+σ)^-t
  w_merged += step*merged.welfare[ix]*(1+σ)^-t
end

println("Welfare analysis:")
println("- w_base:\t\t\t $w_base")
println("- w_fully_anticipated:\t $w_fully_anticipated ($(100*w_fully_anticipated/w_base)%)")
println("- w_merged:\t\t $w_merged ($(100*w_merged/w_base)%)")

# Some plots...
plot(times,base.V[ix], title="Volumes", label="base");
plot!(times,fully_anticipated.V[ix], label="fully anticipated");
plot!(times,merged.V[ix], label="ex-post adapt")

plot(times,base.ben_carbon_seq[ix], title= "Benefit from carbon sequestration",label="base")
plot!(times,fully_anticipated.ben_carbon_seq[ix], label="fully anticipated")
plot!(times,merged.ben_carbon_seq[ix], label="ex-post adapt")

plot(times,base.ben_wood[ix], title= "Benefit from wood harvesting",label="base")
plot!(times,fully_anticipated.ben_wood[ix], label="fully anticipated")
plot!(times,merged.ben_wood[ix], label="ex-post adapt")

plot(times,(base.r_A[ix] .+ base.r_F[ix]), title= "SF area IN",label="base");
plot!(times,(fully_anticipated.r_A[ix] .+ fully_anticipated.r_F[ix]), label="fully anticipated");
plot!(times,(merged.r_A[ix] .+ merged.r_F[ix]), label="ex-post adapt")

plot(times,(base.a[ix]), title= "SF area OUT",label="base");
plot!(times,(fully_anticipated.a[ix]), label="fully anticipated");
plot!(times,(merged.a[ix]), label="ex-post adapt")

plot(times,(base.r_A[ix] .+ base.r_F[ix] .- base.a[ix]), title= "SF net area change",label="base");
plot!(times,(fully_anticipated.r_A[ix] .+ fully_anticipated.r_F[ix] - fully_anticipated.a[ix]), label="fully anticipated");
plot!(times,(merged.r_A[ix] .+ merged.r_F[ix] - merged.a[ix]), label="ex-post adapt")

