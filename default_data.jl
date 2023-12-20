
# ------------------------------------------------------------------------------
# BASE PARAMETERS/SETTINGS

# Init values of state variables...
# Sources:
# - Forest area (ha): FAO FRA 2020 report, Brazil (https://www.fao.org/3/ca9976en/ca9976en.pdf), p 24
# - Forest volumes (m³): idem, p 44
# - Agricultural area (ha): World bank 2020, https://data.worldbank.org/indicator/AG.LND.AGRI.K2?locations=BR
TOT_LAND = (485_396_011 + 11_223_609 + 2_387_471 * 100) / 1_000_000
F₀ = 485_396_011.0 / 1_000_000   # Init Primary Forest area Source: 
S₀ = 11_223_609 / 1_000_000       # Init Secondary Forest area
A₀ = 2_387_471 * 100  / 1_000_000 # Init Agricultural area [ha]
V₀ = 3_056_570_000 / 1_000_000    # Init timber volumes of secondary forests [m^3]
d₀ = F₀/100000        # Init prim forest harvesting area (used only to compute parameter and start value in the sense of optimisation start)
h₀ = S₀/100           # Init sec forest harvesting area (used only to compute parameter and start value in the sense of optimisation start) 
r₀ = 0.5 * d₀         # Init sec forest regeneration area (used only to compute parameter and start value in the sense of optimisation start) 

V₀/S₀
# Parameters
# Sources:
# - bwood_c1(USD/m^3): production and export value from FAOSTAT 2020 https://www.fao.org/faostat/en/#data/FO
# - bagr_c1 (USD/ (ha  y⁻¹) ): FAOSTAT, Gross Production Value (current US$), https://www.fao.org/faostat/en/#data/QV
# - benv_c1 (USD/ (ha y⁻¹) ): Costanza 1997, from https://www.fair-and-precious.org/en/news/349/tropical-forests-the-facts-and-figures
_tot_vol = 121_503_650_000 / 1_000_000 # Total forest tibmer volumes [m^3] (not used in the model, only to compute parameters)
_h_vol   = 290_586_490 / 1_000_000# Total "produced" timber volumes [m^3] (not used in the model, only to compute parameters)
D   = (_tot_vol - V₀) / F₀  # Density of the primary forest (constant) [m^3/ha]
σ          = 0.04 # Discount rate 
K          = 600  # Maximum density of the secondary forest (i.e. carrying capacity)
γ          = 0.1 # Growth rate of the logistic function in terms of density
benv_c2    = 1/2  # Power of the environmental benefit
benv_c1    = 0.3*2700 * F₀ /(F₀^(benv_c2))  # Multiplier of the environmental benefits
bagr_c2    = 1/2  # Power of the agricultural benefits
bagr_c1    = (144_160_490 *1000 * A₀ / (2_387_471 * 100)) / (A₀^bagr_c2) # Multiplier of the agricultural benefits
bwood_c2   = 1/2  # Power of the wood-use benefits
bwood_c1   = 3*100.756304548882 * ((d₀*D + h₀*V₀/S₀)) / ((d₀*D + h₀*V₀/S₀)^bwood_c2 ) # Multiplier of the wood-use benefits
chpf_c2    = 1  # Power of the harvesting costs of primary forest (harvested area)
chpf_c3    = -2   # Power of the harvesting costs of primary forest (primary forest area)
chpf_c1    = 100.0 * d₀ / (d₀*F₀^chpf_c3)   # Multiplier of the harvesting costs of primary forest
chsf_c2    = 0    # Power of the harvesting costs of secondary forest
chsf_c1    = 30.0 * h₀ / (h₀^chsf_c2)   # Multiplier of the harvesting costs of secondary forest
crsf_c2    = 0    # Power of the regeneration costs of secondary forest
crsf_c1    = 100 * r₀ / (r₀^crsf_c2) # Multiplier of the regeneration costs of secondary forest


# Options
opt = Ipopt.Optimizer  # Desired optimizer (solver)
T   = 2000             # Time horizon (years)
ns  = 201;             # Number of points in the time grid - seems not to influence much the results (good!)


