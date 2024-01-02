
# ------------------------------------------------------------------------------
# BASE PARAMETERS/SETTINGS

# Init values of state variables...
# Sources:
# - Forest area (M ha): FAO FRA 2020 report, Brazil (https://www.fao.org/3/ca9976en/ca9976en.pdf), p 24
# - Forest volumes (M m³): idem, p 44
# - Agricultural area (M ha): World bank 2020, https://data.worldbank.org/indicator/AG.LND.AGRI.K2?locations=BR
TOT_LAND = (485_396_011 + 11_223_609 + 2_387_471 * 100) / 1_000_000
F₀ = 485_396_011.0 / 1_000_000   # Init Primary Forest area Source: 
S₀ = 11_223_609 / 1_000_000       # Init Secondary Forest area
A₀ = 2_387_471 * 100  / 1_000_000 # Init Agricultural area [M ha]
V₀ = 3_056_570_000 / 1_000_000    # Init timber volumes of secondary forests [M m^3]
d₀ = F₀/100000        # Init prim forest harvesting area (used only to compute parameter and start value in the sense of optimisation start) [M ha]
h₀ = S₀/100           # Init sec forest harvesting area (used only to compute parameter and start value in the sense of optimisation start)  [M ha]
r₀ = 0.5 * d₀         # Init sec forest regeneration area (used only to compute parameter and start value in the sense of optimisation start)  [M ha]

# Parameters
# Sources:
# - bwood_c1(USD/m^3): production and export value from FAOSTAT 2020 https://www.fao.org/faostat/en/#data/FO
# - bagr_c1 (USD/ (ha  y⁻¹) ): FAOSTAT, Gross Production Value (current US$), https://www.fao.org/faostat/en/#data/QV
# - benv_c1 (USD/ (ha y⁻¹) ): Costanza 1997, from https://www.fair-and-precious.org/en/news/349/tropical-forests-the-facts-and-figures
_tot_vol = 121_503_650_000 / 1_000_000 # Total forest tibmer volumes [M m^3] (not used in the model, only to compute parameters)
_h_vol   = 290_586_490 / 1_000_000 # Total "produced" timber volumes [M m^3] (not used in the model, only to compute parameters)
_agr_unit_value = 144_160_490 *1000  / (A₀* 1_000_000) # $/ha from the value of agr production in K$
D          = (_tot_vol - V₀) / F₀  # Density of the primary forest (constant) [m^3/ha]
σ          = 0.04 # Discount rate 
K          = 600  # Maximum density of the secondary forest (i.e. carrying capacity)
γ          = 0.1 # Growth rate of the logistic function in terms of density
benv_c2    = 1/2  # Power of the environmental benefit
"""
Multiplier of the environmental benefits

benv_c1 is calibrated such that:
ben_env = (0.3*2700) * F₀
That is, the environmental benefits is equal to the €/hectar value from literature multiplied the initial area (in Mha),
Hence it is such that ben_env is in M\$
"""
benv_c1    = 0.3*2700 * F₀ /(F₀^(benv_c2))
bagr_c2    = 1/2  # Power of the agricultural benefits
"""
Multiplier of the agricultural benefits

bagr_c1 is calibrated such that:
ben_agr = (603) * A₀
That is, the agricultural benefits is equal to the €/hectar value from value of the agr production multiplied the initial area (in Mha),
Hence it is such that ben_agr is in M\$
"""
bagr_c1    = _agr_unit_value * A₀  / (A₀^bagr_c2) # Multiplier of the agricultural benefits
bwood_c2   = 1/2  # Power of the wood-use benefits
"""
Multiplier of the wood benefits

bwood_c1 is calibrated such that:
ben_wood = (300) * hV
That is, the wood benefits is equal to the \$/m^3 value from export value multiplied the initial harvested values (in Mm^3),
Hence it is such that ben_wood is in M\$
"""
bwood_c1   = 3*100.756304548882 * (d₀*D + h₀*V₀/S₀) / ((d₀*D + h₀*V₀/S₀)^bwood_c2 ) # Multiplier of the wood-use benefits

bc_c2      = 0             # Carbon price growth rate
"""
Multiplier of the carbon benefits

bc_c1 is calibrated such that:
bc_c1 = (100) * ΔV
That is, the carbon benefits is equal to the \$/m^3 value of (assumed) carbon value multiplied the initial ΔV, the total
delta forest volumes, i.e. natural growt of secondary forest less harvested secondary forest less harvested primary
forest volumes (in Mm^3)
Hence it is such that ben_carbon is in M\$
Note that here we have:
bc_c1*exp(bc_c2 * t) * ΔV  [M\$] =  100 * ΔV [M\$]
So ΔV cancel it out and being interested in start time, so does exp(bc_c2 * t)
"""
bc_c1      = 100   # Carbon price init price
chpf_c2    = 2  # Power of the harvesting costs of primary forest (harvested area)
chpf_c3    = -16   # Power of the harvesting costs of primary forest (primary forest area)
"""
Multiplier of the primary forest harvesting costs

chpf_c1 is calibrated such that:
cost_pfharv = (100) * hv_pf
That is, the timber harvesting cost of primary forests is equal to the \$/m^3 cost (assumed) multiplied by the harvested
volumes of primary forest (in Mm^3)
Hence it is such that cost_pfharv is in M\$
"""
chpf_c1    = 100.0 * (d₀*D) / ((d₀*D)^chpf_c2 * F₀^chpf_c3)   # Multiplier of the harvesting costs of primary forest
chsf_c2    = 0    # Power of the harvesting costs of secondary forest
"""
Multiplier of the secondary forest harvesting costs

chsf_c1 is calibrated such that:
cost_sfharv = (30) * hV_sf
That is, the timber harvesting cost of secondary forests is equal to the \$/m^3 cost (assumed) multiplied by the harvested
volumes of secondary forest (in Mm^3)
Hence it is such that cost_sfharv is in M\$
"""
chsf_c1    = 30.0 * (h₀*V₀/S₀) / ((h₀*V₀/S₀)^chsf_c2)   # Multiplier of the harvesting costs of secondary forest
crsf_c2    = 0    # Power of the regeneration costs of secondary forest
"""
Multiplier of the secondary forest regeneration costs

crsf_c1 is calibrated such that:
cost_sfreg = 100 * r
That is, the forest regeneration cost of secondary forests is equal to the \$/ha cost (assumed) multiplied by the
regeneration area of secondary forest (in Mha)
Hence it is such that cost_sfreg is in M\$
"""
crsf_c1    = 100 * r₀ / (r₀^crsf_c2) # Multiplier of the regeneration costs of secondary forest

# Options
optimizer   = Ipopt.Optimizer  # Desired optimizer (solver)
opt_options = Dict("max_cpu_time" => 20.0)
T           = 2000             # Time horizon (years)
ns          = 201;             # Number of points in the time grid - seems not to influence much the results (good!)


