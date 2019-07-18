###########################################
# define functions of pysiology processes #
###########################################
respir_extra(size) = respir_ex*size^respir_b # extra cost for biosynthesis, related to Cq1
function k_respir(age)
    k_res = respir_a*age^respir_ab # respiration rate, unit: per second
    return k_res
end

function daynight(t::Int64, IR)
    if IR[trunc(Int,t*ΔT/3600)] < 5.0
        return false
    else
        return true
    end
end

function PAR_cal(I, z, cumsum_chl)
    atten = (katten_w + katten_c * cumsum_chl)*(-z)
    PAR = α*I*(1.0 - exp(-atten))/atten
    return PAR
end

function PC(PAR, Temp, phyt)  
    Tempstd = exp(TempAe*(1.0/(Temp+273.15)-1.0/Tempref))
    photoTempFunc = TempCoeff*max(1.0e-10,Tempstd)
    PCm = PCmax[phyt.sp]*phyt.size^PC_b[phyt.sp]
    PC = PCm*photoTempFunc*(1-exp(-PAR*phyt.chl/(phyt.Cq2*PCm*photoTempFunc)))
    PS = PC*Cquota[phyt.sp]*Nn
    return PS # unit: mmol C/second/individual
end

function Nuptake(Nit, phyt)
    Qn = phyt.Nq/(phyt.Cq1+phyt.Cq2)
    #In-Cell N uptake limitation
    regQ = max(0.0,min(1.0,(Nqmax-Qn)/(Nqmax-Nqmin)))
    VNm = VNmax[phyt.sp]*phyt.size^VN_b[phyt.sp]
    Nuptake = VNm*Nit/(Nit+KsatN)*regQ
    VNcell = Nuptake*Cquota[phyt.sp]*Nn
    return VNcell # unit: mmol N/second/individual
end

function chl_sync(phyt,PP,I)
    if I > 0
        ρ_chl = PP/(α*I*phyt.chl/phyt.Cq2)
    else
        ρ_chl = 0.0
    end
    return ρ_chl
end

function divide(phyt::DataFrameRow)
    phytops = DataFrame(x=Float64[0.0,0.0], y=Float64[0.0,0.0], z=Float64[0.0,0.0], gen=Int64[1,1], size=Float64[0.0,0.0], Cq1=Float64[0.0,0.0], Cq2=Float64[0.0,0.0], Nq=Float64[0.0,0.0], chl=Float64[0.0,0.0], sp=Int64[0,0], age=Float64[0.0,0.0])  # initialize new cell
    # NOT all C and N can turn into new cells
    phytops[1,:].x = phyt.x
    phytops[1,:].y = phyt.y
    phytops[1,:].z = phyt.z
    phytops[1,:].gen = phyt.gen + 1
    phytops[1,:].Cq1 = phyt.Cq1 * 0.5
    phytops[1,:].Cq2 = phyt.Cq2 * 0.45
    phytops[1,:].Nq  = phyt.Nq  * 0.5
    phytops[1,:].size= phyt.size* 0.5
    phytops[1,:].chl = phyt.chl * 0.5
    phytops[1,:].sp = phyt.sp
    phytops[1,:].age = 1.0

    phytops[2,:].x = phyt.x
    phytops[2,:].y = phyt.y
    phytops[2,:].z = phyt.z
    phytops[2,:].gen = phyt.gen + 1
    phytops[2,:].Cq1 = phyt.Cq1 * 0.5
    phytops[2,:].Cq2 = phyt.Cq2 * 0.45
    phytops[2,:].Nq  = phyt.Nq  * 0.5
    phytops[2,:].size= phyt.size* 0.5
    phytops[2,:].chl = phyt.chl * 0.5
    phytops[2,:].sp = phyt.sp
    phytops[2,:].age = 1.0

    return phytops
end
###################################
# model update for phytoplanktons #
###################################
function phyt_update(t::Int64, ΔT::Int64, g, phyts_a, nutrients, IR, temp)
    # load nutrients
    dvid_ct = 0; graz_ct = 0; death_ct = 0
    Num_phyt = size(phyts_a,1)
    chl_num = count_chl(phyts_a, g)
    cumsum_chl = cumsum(chl_num, dims = 3)
    #set up a dataframe to record all updated agents
    phyts_b = DataFrame(x=Float64[], y=Float64[], z=Float64[], gen=Int64[], size=Float64[], Cq1=Float64[], Cq2=Float64[], Nq=Float64[], chl=Float64[], sp=Int64[], age=Float64[])
    #
    consume = nutrient_fields(zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz), zeros(g.Nx, g.Ny, g.Nz))
    # iterate phytoplankton agents
    #
    for i in 1:size(phyts_a,1)
        phyt = phyts_a[i,:]
        sp = phyt.sp
        z = trunc(Int, phyt.z); x = trunc(Int, phyt.x); y = trunc(Int, phyt.y);
        DIN = max(0.0, nutrients.DIN[x, y, z])
        # Probability for a cell to do nothing in a time step.
        P_dormant = false #rand(Bernoulli(Dor_P))
        # compute probabilities of grazing and division
        P_graz = rand(Bernoulli(exp(Num_phyt/N*Nsp)*phyt.size/Grz_P))
        # Hypothesis: the population of grazers is large enough to graze on phytoplanktons
        reg_size = max(0.0, phyt.size-dvid_size)
        P_dvi=rand(Bernoulli((a_dvi[sp]*reg_size)^b_dvi[sp]/(1+(a_dvi[sp]*reg_size)^b_dvi[sp]))) 
        PAR = PAR_cal(IR[trunc(Int,t*ΔT/3600)], g.zF[z], cumsum_chl[x, y, z])
        PP = PC(PAR,temp[trunc(Int,t*ΔT/3600)],phyt)*ΔT
        VN = min(DIN*g.V[x,y,z]/10.0,Nuptake(DIN,phyt)*ΔT)
        Dmd_NC = (1+respir_extra(phyt.Cq1))*VN/R_NC
        Res = k_respir(phyt.age)*phyt.Cq2*ΔT
        ρ_chl = chl_sync(phyt,PC(PAR,temp[trunc(Int,t*ΔT/3600)],phyt),IR[trunc(Int,t*ΔT/3600)])
        if P_graz == false #not grazed
            if P_dormant == false # not dormant
                if daynight(t,IR) # day
                    if PP > Dmd_NC
                        ExuC = (PP - Dmd_NC)*FracExuC
                        CostC= Dmd_NC
                        SynC = VN/R_NC
                    else
                        CostC= PP
                        SynC = PP/(1+respir_extra(phyt.size))
                        ExuC = 0.0
                    end #exudation
                else # night
                    CostC= min(0.2*phyt.Cq1,Dmd_NC)
                    SynC = min(VN/R_NC,0.2*phyt.Cq1/(1+respir_extra(phyt.size)))
                    ExuC = 0.0
                end # day night?
                dCq1 = PP - CostC - ExuC
                dCq2 = SynC - Res
                dNq  = VN
                dsize= dCq2/phyt.Cq2
                phyt.Cq1 = max(Cmin[phyt.sp]*Nn/10.0, phyt.Cq1 + dCq1)
                phyt.Cq2 = max(Cmin[phyt.sp]*Nn/10.0, phyt.Cq2 + dCq2)
                phyt.Nq  = max(Cmin[phyt.sp]*Nn*R_NC/10.0, phyt.Nq + dNq)
                phyt.size= max(0.0,phyt.size+dsize)
                phyt.chl = phyt.chl + ρ_chl*VN*Chl2N
                phyt.age = phyt.age + 1.0*(ΔT/3600)
                if (phyt.Cq2+phyt.Cq1 ≥ Cmin[phyt.sp]*Nn) & (phyt.size > 0) # not natural death
                    if P_dvi < 1 # not divide
                        push!(phyts_b,phyt)
                    else # divide
                        dvid_ct += 2
                        global dvdcount += 2
                        phyts2 = divide(phyt)
                        append!(phyts_b,phyts2)
                        consume.DIC[x, y, z] = consume.DIC[x, y, z] + phyt.Cq2*0.1 # consume C when cell is divided
                    end # divide
                else # natural death
                    consume.DOC[x, y, z] = consume.DOC[x, y, z] + (phyt.Cq1+phyt.Cq2)*mortFracC
                    consume.DON[x, y, z] = consume.DON[x, y, z] + phyt.Nq*mortFracN
                    consume.POC[x, y, z] = consume.POC[x, y, z] + (phyt.Cq1+phyt.Cq2)*(1.0 - mortFracC)
                    consume.PON[x, y, z] = consume.PON[x, y, z] + phyt.Nq*(1.0 - mortFracN)
                    death_ct += 1
                end # naturan death
                consume.DIC[x, y, z] = consume.DIC[x, y, z] + Res + CostC - SynC
                consume.DIN[x, y, z] = consume.DIN[x, y, z] - VN
                consume.DOC[x, y, z] = consume.DOC[x, y, z] + ExuC
            else # dormant, do nothing in this time step
                dCq2 = -Res
                dsize= dCq2/phyt.Cq2
                phyt.Cq2 = max(Cmin[phyt.sp]*Nn/50.0,phyt.Cq2 + dCq2)
                phyt.size= max(0.0,phyt.size+dsize)
                phyt.age = phyt.age + 1.0*(ΔT/3600)
                push!(phyts_b,phyt)
            end # dormant
        else #grazed
            graz_ct += 1
            consume.DOC[x, y, z] = consume.DOC[x, y, z] + (phyt.Cq1+phyt.Cq2)*grazFracC*0.5
            consume.DON[x, y, z] = consume.DON[x, y, z] + phyt.Nq*grazFracN*0.5
            consume.POC[x, y, z] = consume.POC[x, y, z] + (phyt.Cq1+phyt.Cq2)*(1.0 - grazFracC)*0.5
            consume.PON[x, y, z] = consume.PON[x, y, z] + phyt.Nq*(1.0 - grazFracN)*0.5
        end # graze
    end # while loop to traverse the array of agents
    return phyts_b,dvid_ct,graz_ct,death_ct,consume
end # for loop of time
