##### temperature function
# @inline function tempFunc(temp, p)
#     TFunc = max(1.0e-30, exp(p.TAe * (1.0/(temp+273.15) - (1.0/p.TRef)))) * p.TCoeff
#     return TFunc
# end
@inline function tempFunc(temp, p)
    k = exp(-p.Ea/(8.3145*(temp+273.15)))*(1.0-exp(temp - p.T⁺))
    k = max(0.0, k)
    OGT_rate = exp(-p.Ea/(8.3145*(p.T⁺+273.15-2)))
    return k/OGT_rate
end

##### calculate photosynthesis rate (mmolC/individual/second)
@inline function calc_PS(par, temp, Chl, Bm, p)
    αI  = par * p.α * p.Φ
    PCm = p.PCmax * tempFunc(temp, p)
    PS  = PCm * (1.0 - exp(-αI / max(1.0e-30, PCm) * Chl / max(1.0e-30, Bm))) * Bm
    return PS
end

##### calculate nutrient uptake rate (mmolN/individual/second)
@inline function calc_NP_uptake(NH4, NO3, PO4, temp, Cq, Nq, Pq, Bm, p, ac)
    Qn = (Nq + Bm * p.R_NC)/max(1.0e-30, Bm + Cq)
    Qp = (Pq + Bm * p.R_PC)/max(1.0e-30, Bm + Cq)
    fQn = max(0.0, min(1.0, 1.0 - Qn / p.Nqmax))
    fQp = max(0.0, min(1.0, 1.0 - Qp / p.Pqmax))
    regQN = fQn^4 / (1.0e-4 + fQn^4)
    regQP = fQp^4 / (1.0e-4 + fQp^4)
    VNH4 = p.VNH4max * regQN * NH4/max(1.0e-30, NH4+p.KsatNH4) * tempFunc(temp, p) * Bm * ac
    VNO3 = p.VNO3max * regQN * NO3/max(1.0e-30, NO3+p.KsatNO3) * tempFunc(temp, p) * Bm * ac
    VPO4 = p.VPO4max * regQP * PO4/max(1.0e-30, PO4+p.KsatPO4) * tempFunc(temp, p) * Bm * ac
    return VNH4, VNO3, VPO4
end

@inline function calc_DOC_uptake(DOC, temp, Cq, Bm, p)
    Qc = Cq/max(1.0e-30, Bm + Cq)
    fQc = max(0.0, min(1.0, 1.0 - Qc / p.Cqmax))
    regQ = fQc^4 / (1.0e-4 + fQc^4)
    VN = p.VDOCmax * regQ * DOC/max(1.0e-30, DOC+p.KsatDOC) * tempFunc(temp, p) * Bm
    return VN
end

@kernel function calc_inorganic_uptake_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.PS[i] = calc_PS(nuts.par[i], nuts.T[i], plank.Chl[i], plank.Bm[i], p) * plank.ac[i]

    @inbounds plank.VNH4[i], plank.VNO3[i], plank.VPO4[i] =
                        calc_NP_uptake(nuts.NH4[i], nuts.NO3[i], nuts.PO4[i], nuts.T[i],
                                        plank.Cq[i], plank.Nq[i], plank.Pq[i],
                                        plank.Bm[i], p, plank.ac[i])
end
function calc_inorganic_uptake!(plank, nuts, p, arch::Architecture)
    kernel! = calc_inorganic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p)
    return nothing
end

##### update C, N, P quotas for the first time of each time step
@kernel function update_quotas_1_kernel!(plank, ΔT)
    i = @index(Global)
    @inbounds plank.Cq[i] += plank.PS[i]   * ΔT
    @inbounds plank.Pq[i] += plank.VPO4[i] * ΔT
    @inbounds plank.Nq[i] +=(plank.VNH4[i] + plank.VNO3[i]) * ΔT
end
function update_quotas_1!(plank, ΔT, arch)
    kernel! = update_quotas_1_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT)
    return nothing
end

##### calculate DOC uptake rate (mmolC/individual/second)
##### DOC uptake needs support of photosynthesis for at least 1% of total C acquisition.
@kernel function calc_organic_uptake_kernel!(plank, nuts, p)
    i = @index(Global)
    @inbounds plank.VDOC[i] = calc_DOC_uptake(nuts.DOC[i], nuts.T[i], plank.Cq[i], plank.Bm[i], p) * plank.ac[i]

    @inbounds plank.VDOC[i] = plank.VDOC[i] * isless(0.01, plank.PS[i]/(plank.VDOC[i]+plank.PS[i]))
end
function calc_organic_uptake!(plank, nuts, p, arch::Architecture)
    kernel! = calc_organic_uptake_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, p)
    return nothing
end

##### calculate ρChl
@kernel function calc_ρChl_kernel!(plank, par, p)
    i = @index(Global)
    @inbounds plank.ρChl[i] = plank.PS[i]/max(1.0e-30, plank.Bm[i]) * p.Chl2N / 
                             max(1.0e-30, par[i] * p.α * p.Φ * plank.Chl[i]/plank.Bm[i]) *
                             isless(1.0e-1, par[i]) * plank.ac[i]
end
function calc_ρChl!(plank, par, p, arch)
    kernel! = calc_ρChl_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, par, p)
    return nothing
end

##### calculate respiration (mmolC/individual/second)
@kernel function calc_respir_kernel!(plank, T, p)
    i = @index(Global)
    @inbounds plank.resp[i] = p.respir_a * plank.Bm[i] * tempFunc(T[i], p) * plank.ac[i]
end
function calc_respir!(plank, T, p, arch)
    kernel! = calc_respir_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, T, p)
    return nothing
end


##### update C, N, P quotas for the second time of each time step
@kernel function update_quotas_2_kernel!(plank, ΔT, p)
    i = @index(Global)
    @inbounds plank.Cq[i] = plank.Cq[i] + (plank.VDOC[i] - plank.resp[i]) * ΔT
    @inbounds plank.Cq[i] = plank.Cq[i] + max(0.0, (0.0 - plank.Cq[i]))
    @inbounds plank.Bm[i] = plank.Bm[i] - max(0.0, (0.0 - plank.Cq[i]))
    @inbounds plank.Nq[i] = plank.Nq[i] + max(0.0, (0.0 - plank.Cq[i])) * p.R_NC
    @inbounds plank.Pq[i] = plank.Pq[i] + max(0.0, (0.0 - plank.Cq[i])) * p.R_PC
end
function update_quotas_2!(plank, ΔT, p, arch)
    kernel! = update_quotas_2_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, ΔT, p)
    return nothing
end

##### calculate biosynthesis and exudation (mmolC/individual/second)
@kernel function calc_BS_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.BS[i] = min(plank.Cq[i], plank.Nq[i]/p.R_NC, plank.Pq[i]/p.R_PC) * 
                           p.k_mtb * plank.ac[i]
    @inbounds plank.exu[i]= max(0.0, plank.Cq[i] - min(plank.Cq[i], plank.Nq[i]/p.R_NC, plank.Pq[i]/p.R_PC)) *
                           p.k_mtb * plank.ac[i]
end
function calc_BS!(plank, p, arch)
    kernel! = calc_BS_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### update C, N, P quotas, biomass, Chla
@kernel function update_biomass_kernel!(plank, p, ΔT)
    i = @index(Global)
    @inbounds plank.Bm[i]  += ΔT * plank.BS[i]
    @inbounds plank.Cq[i]  -= ΔT *(plank.BS[i] + plank.exu[i])
    @inbounds plank.Nq[i]  -= ΔT * plank.BS[i] * p.R_NC
    @inbounds plank.Pq[i]  -= ΔT * plank.BS[i] * p.R_PC
    @inbounds plank.Chl[i] += ΔT * plank.BS[i] * p.R_NC * plank.ρChl[i]
    @inbounds plank.age[i] += ΔT / 3600.0 * plank.ac[i]
end
function update_biomass!(plank, p, ΔT, arch)
    kernel! = update_biomass_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p, ΔT)
    return nothing
end

##### update cell size
##### only calculate functional biomass
@kernel function update_cellsize_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.Sz[i] = plank.Bm[i] / (p.Cquota * p.Nsuper)
end
function update_cellsize!(plank, p, arch)
    kernel! = update_cellsize_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end
