# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.1
#   kernelspec:
#     display_name: Julia 1.1.0
#     language: julia
#     name: julia-1.1
# ---

using DataFrames, NetCDF, Printf, CSV, Serialization
using Random
using Distributions
cd("/nobackup1b/users/zhenwu/ABPM_3D/")
include("parameters.jl")
include("model_setup.jl")
include("model_struct.jl")
include("phyt_process.jl")
include("utils.jl")
include("agent_div.jl")
include("dst3fl.jl")
include("nutrient_processes.jl")
include("2nd_adv_diffu.jl")
# remove old files
isfile("results/cons_C.txt") && rm("results/cons_C.txt");
isfile("results/cons_N.txt") && rm("results/cons_N.txt");
isfile("results/cons_DIN.txt") && rm("results/cons_DIN.txt");
isfile("results/B1.bin") && rm("results/B1.bin");
isfile("results/B2.bin") && rm("results/B2.bin");
isfile("results/output.bin") && rm("results/output.bin");
isfile("results/output1.bin") && rm("results/output1.bin");
isfile("results/output2.bin") && rm("results/output2.bin");
isfile("results/grid.bin") && rm("results/grid.bin");
isfile("results/IR.bin") && rm("results/IR.bin");
isfile("results/VD1.bin") && rm("results/VD1.bin");
isfile("results/VD2.bin") && rm("results/VD2.bin");
isfile("results/HD1.bin") && rm("results/HD1.bin");
isfile("results/HD2.bin") && rm("results/HD2.bin");
# Read input files
nTime = 240 # number of time steps
ΔT = 3600 # time step: 3600 for 1 hour
temp,IR = read_input("T_IR.csv",trunc(Int,nTime*ΔT/3600));

# grid selected : [251:750,1251:1800]: 27.3242N to 36.2605N, 167.802W to 157.406W
fieldroot = "/nobackup1b/users/jahn/hinpac/grazsame3/run/run.0354/";
g = grid_offline(fieldroot);

# deal with time steps of offline velocityfields
itvalLo = 144;
itvalHi = 687888;
itList = collect(itvalLo:144:itvalHi);
tN = 4056; # starting time
vfroot = "/nobackup1b/users/jahn/hinpac/grazsame3/run/run.0354/offline-0604/"; # directory of velocity fields

N = 100000   # Number of initial individuals of each species
Nsp = 2     # Number of species
Nn = Int(1e15)   # Number of cells one super-agent represents
B=setup_agents(N,Cquota,Nn,1.1,0.18,g) # Normal distribution with mean and variance
# model initialization
# create output file
output = create_output(B);
nut = [2.0, 0.05, 20.0, 0.0, 0.0, 0.0] #DIC, DIN, DOC, DON, POC, PON, mmol/m3
nutrients= setup_nutrients(g,nut)
remin = rem(kDOC,kDON,kPOC,kPON);

for t in 1:nTime
    phyts_a = copy(B[t]) # read data from last time step
    phyts_b,dvid_ct,graz_ct,death_ct,consume=phyt_update(t, ΔT, g, phyts_a, nutrients, IR, temp)
    velᵇ = read_offline_vels(vfroot,itList,tN,trunc(Int,t*ΔT/3600));
    velᵈ = double_grid(velᵇ,g)
    agent_move(phyts_b,velᵈ,g,ΔT) 
    push!(B,phyts_b)
    write_output(t,phyts_b,dvid_ct,graz_ct,death_ct,output)
    agent_num = size(phyts_b,1)
    F = compute_nut_biochem(nutrients, remin)
    gtr = compute_source_term(nutrients, velᵇ, g, F)
    nutₜ = nut_update(nutrients, consume, g, gtr, ΔT)
    write_nut_nc(g, nutₜ, t)
    write_nut_cons(g, gtr, nutₜ, velᵇ, agent_num, t)
    global nutrients = nutₜ;
end

B1 = []; B2 = [];
for i in 1:size(B,1)
    sort_species(B[i], B1, B2)
end

HD1 = []; HD2 = [];
for i in 1:size(B,1)
    HD_1 = count_horizontal_num(B1[i],g);
    push!(HD1,HD_1)
    HD_2 = count_horizontal_num(B2[i],g);
    push!(HD2,HD_2)
end

for i in 1:size(B,1)
    convert_coordinates(B1[i],g) # convert grids to lon, lat and depth
    convert_coordinates(B2[i],g) # convert grids to lon, lat and depth
end

VD1 = []; VD2 = [];
for i in 1:size(B,1)
    VD_1 = count_vertical_num(B1[i]);
    push!(VD1,VD_1)
    VD_2 = count_vertical_num(B2[i]);
    push!(VD2,VD_2)
end

output1, output2 = compute_mean_species(B1, B2, nTime);

# save model output
open("results/B1.bin", "w") do io
    serialize(io, B1)
end
open("results/B2.bin", "w") do io
    serialize(io, B2)
end
open("results/grid.bin", "w") do io
    serialize(io, g)
end
open("results/output.bin", "w") do io
    serialize(io, output)
end
open("results/output1.bin", "w") do io
    serialize(io, output1)
end
open("results/output2.bin", "w") do io
    serialize(io, output2)
end
open("results/VD1.bin", "w") do io
    serialize(io, VD1)
end
open("results/VD2.bin", "w") do io
    serialize(io, VD2)
end
open("results/HD1.bin", "w") do io
    serialize(io, HD1)
end
open("results/HD2.bin", "w") do io
    serialize(io, HD2)
end
open("results/IR.bin", "w") do io
    serialize(io, IR)
end

#=t = 1
phyts_a = copy(B[t]) # read data from last time step
phyts_b,dvid_ct,graz_ct,death_ct,consume=phyt_update(t, ΔT, g, phyts_a, nutrients, IR, temp)
velᵇ = read_offline_vels(vfroot,itList,tN,trunc(Int,t*ΔT/3600));
velᵈ = double_grid(velᵇ,g)
agent_move(phyts_a,velᵈ,g,ΔT) 
push!(B,phyts_b)
write_output(t,phyts_b,dvid_ct,graz_ct,death_ct,output)
agent_num = size(phyts_b,1)
F = compute_nut_biochem(nutrients, remin)
gtr = compute_source_term(nutrients, velᵇ, g, F)
nutₜ = nut_update(nutrients, consume, g, gtr, ΔT)
global nutrients = nutₜ;

maximum(consume.DON)

findall(isequal(minimum(.POC)),F.POC)

Btest=DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[], sp=Int64[])
Btest1=DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[], sp=Int64[])
Btest2=DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[], sp=Int64[]);

for i in findall(x -> 187 > x >184,phyts_a.x)
    push!(Btest,phyts_a[i,:])
end

for i in findall(x -> 262 > x > 260,Btest.y)
    push!(Btest1,Btest[i,:])
end

agent_move(Btest1,velᵈ,g,ΔT) 

phyts_a = copy(B[7])
Num_phyt = size(phyts_a,1)
chl_num = count_chl(phyts_a, g)
cumsum_chl = cumsum(chl_num, dims = 3);

phyt =phyts_a[615,:]

z = trunc(Int, phyt.z); x = trunc(Int, phyt.x); y = trunc(Int, phyt.y);
DIN = max(0.0, nutrients.DIN[x, y, z])
#compute probabilities of grazing and division
P_graz = rand(Bernoulli(exp(Num_phyt/N*Nsp)*phyt.size/Grz_P))
# Hypothesis: the population of grazers is large enough to graze on phytoplanktons
P_dvi=max(0.0,phyt.size-dvid_size)*1.0e5*rand(Bernoulli(phyt.size/Dvid_P))
PAR = PAR_cal(IR[trunc(Int,t*ΔT/3600)], g.zF[z], cumsum_chl[x, y, z])
PP = PC(PAR,temp[trunc(Int,t*ΔT/3600)],phyt)*phyt.Cq2*ΔT
VN = min(DIN*g.V[x,y,z]/10.0,Nuptake(DIN,phyt)*phyt.Cq2*ΔT)
Dmd_NC = (1+respir_extra(phyt.Cq1))*VN/R_NC
Res2 = k_respir(phyt.size)*phyt.Cq2*ΔT
ρ_chl = chl_sync(phyt,PC(PAR,temp[trunc(Int,t*ΔT/3600)],phyt),IR[trunc(Int,t*ΔT/3600)])

if PP > Dmd_NC
                    ExuC = (PP - Dmd_NC)*FracExuC
                    CostC= Dmd_NC
                    SynC = VN/R_NC
                else
                    CostC= PP
                    SynC = PP/(1+respir_extra(phyt.Cq1))
                    ExuC = 0.0
                    VN   = SynC*R_NC
                end #exudation

VN



DIN

dCq1 = PP - CostC - ExuC
dCq2 = SynC - Res2
dNq  = VN - Res2*R_NC
dsize= dCq2/phyt.Cq2

phyt.Cq2+phyt.Cq1 ≥ Cmin*Nn

consume.DIC[x, y, z] = consume.DIC[x, y, z] + Res2 + CostC - SynC
                consume.DIN[x, y, z] = consume.DIN[x, y, z] - VN + Res2*R_NC
                consume.DOC[x, y, z] = consume.DOC[x, y, z] + ExuC

consume.DIC[x, y, z]=#


