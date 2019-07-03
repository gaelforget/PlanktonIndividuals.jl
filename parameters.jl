####################################### Parameters ########################################
PCmax = [1.3, 1.6] ./86400   # Maximum primary production rate (per second)
PC_b = [0.15, -0.2]         # Shape parameter
Chl2N = 3.0                  # Maximum Chla:N ratio in phytoplankton
R_NC  = 16/106               # N:C ratio in cell biomass, should be lower than 16:106
Cmin  = 5.0e-13              # Minimum C quota (mmolC)
dvdcount = 0
α = 0.030                    # Irradiance absorption coeff (m^2/gChl)
katten_w = 0.046             # PAR attenuation (/m)
katten_c = 0.04              # PAR attenuation (/mgChl/m^3/m)
Cquota = [1.8e-11, 1.8e-10]  # Average C quota in cell (mmolC).
Grz_P  = 2000                # phyt.size/Grz_P is the probability to be grazed
dvid_size = 0.5              # relative cell size a cell can start divide
a_dvi = [0.20, 0.20]         # shape parameter for division probability
b_dvi = [2.3, 2.3]           # shape parameter for division probability
Dor_P = 0.1                  # probability for a cell to do nothing in a time step

Tempref = 293.15   # reference temperature in K
TempAe = -4000.0   # Arrenhius equation
TempCoeff = 0.8    # Arrenhius equation

VNmax = [0.8, 0.9] ./ 86400 # Maximum N uptake rate (per second)
Nqmax = 0.5               # Maximum N quota in cell (mmol N/mmol C)
Nqmin = 0.13              # Minimum N quota in cell (mmol N/mmol C)
KsatN = 0.01                # Half-saturation coeff
VN_b = [0.15, -0.2]        # Shape parameter

respir_a = 5.0e-9  # Respiration ratio 
respir_ex= 3.0e-11 # Extra cost of C for biosynthesis
respir_b = 0.93    # Shape parameter

grazFracC = 0.7 # Fraction goes into dissolved organic pool
grazFracN = 0.7 # Fraction goes into dissolved organic pool
mortFracC = 0.5 # Fraction goes into dissolved organic pool
mortFracN = 0.5 # Fraction goes into dissolved organic pool
FracExuC  = 0.7 # Fraction of extra fixed carbon to be exuded

k_sink = 0.01/86400 # Sink rates for agents (m/s)
kDOC   = 1/40/86400 # remineralization rate for DOC, turn over time: 40 days (per second)
kDON   = 1/30/86400 # remineralization rate for DON, turn over time: a month (per second)
kPOC   = 1/30/86400 # remineralization rate for POC, turn over time: a month (per second)
kPON   = 1/30/86400 # remineralization rate for PON, turn over time: a month (per second)

κh = 0
κv = 0
####################################### end parameters ##################################
