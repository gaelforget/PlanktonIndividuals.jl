####################################### Parameters ########################################
PCmax = [1.4, 1.4] ./86400   # Maximum primary production rate (per second)
PC_b = [0.6, 0.6]            # Shape parameter for size
Chl2N = 3.0                  # Maximum Chla:N ratio in phytoplankton
R_NC  = 16/106               # N:C ratio in cell biomass, should be lower than 16:106
α = 0.030                    # Irradiance absorption coeff (m^2/gChl)
katten_w = 0.046             # PAR attenuation (/m)
katten_c = 0.04              # PAR attenuation (/mgChl/m^3/m)
Cquota = [1.8e-11, 1.8e-10]  # Average C quota in cell (mmolC).
Grz_P  = 2000                # phyt.size/Grz_P is the probability to be grazed
dvid_size = 0.9              # relative cell size a cell can start divide
a_dvi = [0.35, 0.35]         # shape parameter for division probability
b_dvi = [2.1, 2.1]           # shape parameter for division probability
death_age = 100.0             # average life time of a cell (hour)
a_death = 0.05               # shape parameter for natural death
b_death = 0.5                # shape parameter for natural death

Tempref = 293.15   # reference temperature in K
TempAe = -4000.0   # Arrenhius equation
TempCoeff = 0.8    # Arrenhius equation

VNmax = [0.8, 0.8] ./ 86400 # Maximum N uptake rate (mmol N/mmol C/second)
Nqmax = 0.5                 # Maximum N quota in cell (mmol N/mmol C)
Nqmin = 0.13                # Minimum N quota in cell (mmol N/mmol C)
KsatN = [0.01, 0.015]       # Half-saturation coeff
VN_b = [0.6, 0.6]           # Shape parameter for size

a_β = 2.8       # scale parameter for metabolic partitioning of biosynthesis
b_β = -2.5      # shape parameter for metabolic pratitioning of biosynthesis
k_mtb = 0.15    # metabolic rate (per hour)

respir_ex= 3.0e-4    # Extra cost of C for biosynthesis
respir_b = 0.13      # Shape parameter for size

grazFracC = 0.7 # Fraction goes into dissolved organic pool
grazFracN = 0.7 # Fraction goes into dissolved organic pool
mortFracC = 0.5 # Fraction goes into dissolved organic pool
mortFracN = 0.5 # Fraction goes into dissolved organic pool

k_sink = 0.01/86400 # Sink rates for agents (m/s)
kDOC   = 1/40/86400 # remineralization rate for DOC, turn over time: 40 days (per second)
kDON   = 1/30/86400 # remineralization rate for DON, turn over time: a month (per second)
kPOC   = 1/30/86400 # remineralization rate for POC, turn over time: a month (per second)
kPON   = 1/30/86400 # remineralization rate for PON, turn over time: a month (per second)

κh = 0
κv = 0
####################################### end parameters ##################################