# Forest dynamics numerical model

cd(@__DIR__)
using Pkg
Pkg.activate(".")
#Pkg.resolve()
#Pkg.Instantiate()

using Plots
using Markdown
using Revise

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

out = luc_model()
ix = 1:(length(out.support)-300) # index for plotting
times = out.support[ix]

# Load the "base" optimisation...
#out = luc_model(ns=1001,opt_options = Dict("max_cpu_time" => 60.0))

#ix = 1:(length(out.support)-750) # index for plotting
#times = out.support[ix]

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

