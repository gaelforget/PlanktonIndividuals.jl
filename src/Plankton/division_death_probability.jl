abstract type division_type end
struct Sizer <: division_type end
struct Adder <: division_type end
struct Timer <: division_type end
struct Sizer_Timer <: division_type end
struct Adder_Timer <: division_type end

##### calculate probability of cell division
@inline calc_division(::Sizer, Sz, iS, t, stp, reg, stp2, reg2, P) = P * (tanh(stp*(Sz - reg)) + 1.0)

@inline calc_division(::Adder, Sz, iS, t, stp, reg, stp2, reg2, P) = P * (tanh(stp*(Sz - iS - reg)) + 1.0)

@inline calc_division(::Timer, Sz, iS, t, stp, reg, stp2, reg2, P) = P * (tanh(stp*(t%86400/3600 - reg)) + 1.0)

@inline calc_division(::Sizer_Timer, Sz, iS, t, stp, reg, stp2, reg2, P) = 
                        P * (tanh(stp*(Sz - reg)) + 1.0) * (tanh(stp2*(t%86400/3600 - reg2)) + 1.0)

@inline calc_division(::Adder_Timer, Sz, iS, t, stp, reg, stp2, reg2, P) = 
                        P * (tanh(stp*(Sz - iS - reg)) + 1.0) * (tanh(stp2*(t%86400/3600 - reg2)) + 1.0)

@inline function divide_type(dvid_type)
    if dvid_type == 1
        return Sizer()
    elseif dvid_type == 2
        return Adder()
    elseif dvid_type == 3
        return Timer()
    elseif dvid_type == 4
        return Sizer_Timer()
    elseif dvid_type == 5
        return Adder_Timer()
    else
        throw(ArgumentError("Wrong cell division type, must be in 1 to 5"))
    end
end

@kernel function calc_dvid_kernel!(plank, dvid_type, p, t)
    i = @index(Global)
    @inbounds plank.dvid[i] = calc_division(dvid_type, plank.Sz[i], plank.iS[i], t, 
                                            p.dvid_stp, p.dvid_reg, p.dvid_stp2, p.dvid_reg2, p.dvid_P) * 
                                            isless(2.0 * p.Cquota * p.Nsuper, plank.Bm[i]) * plank.ac[i]
end
function calc_dvid!(plank, dvid_type, p, t, arch)
    kernel! = calc_dvid_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, dvid_type, p, t)
    return nothing
end

@kernel function calc_MM_dvid_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.dvid[i] = p.dvid_P * isless(2.0, plank.DNA[i]/(p.C_DNA*p.Nsuper))
end
function calc_MM_dvid!(plank, p, arch)
    kernel! = calc_MM_dvid_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate the probability of grazing
##### quadratic grazing
@kernel function calc_graz_quadratic_kernel!(plank, nuts, P)
    i = @index(Global)
    @inbounds plank.graz[i] = nuts.pop[i] * P * plank.ac[i]
end
function calc_graz_quadratic!(plank, nuts, P, arch)
    kernel! = calc_graz_quadratic_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, nuts, P)
    return nothing
end

##### calculate the probability of mortality
@kernel function calc_mort_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.mort[i] = p.mort_P * (tanh(20.0*(p.mort_reg - plank.Sz[i])) + 1.0) * plank.ac[i]
end
function calc_mort!(plank, p, arch)
    kernel! = calc_mort_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate the probability of mortality caused by thermal exposure
@kernel function calc_thermal_mort_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.mort[i] = p.mort_P * (tanh(p.Th_stp*(plank.Th[i] - p.Th_reg)) + 1.0) * plank.ac[i]
end
function calc_thermal_mort!(plank, p, arch)
    kernel! = calc_thermal_mort_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### calculate the probability of mortality for MMM
@kernel function calc_MM_mort_kernel!(plank, p)
    i = @index(Global)
    @inbounds plank.mort[i] = p.mort_P * (tanh(1.0*(plank.age[i] - p.mort_reg)) + 1.0) * plank.ac[i]
end
function calc_MM_mort!(plank, p, arch)
    kernel! = calc_MM_mort_kernel!(device(arch), 256, (size(plank.ac,1)))
    kernel!(plank, p)
    return nothing
end

##### generate the random results from probabilities of grazing, mortality and cell division
function get_probability!(plank, rnd, ΔT, arch)
    ##### generate random numbers (0,1) 
    rand!(rng_type(arch), rnd.x)
    rand!(rng_type(arch), rnd.y)
    rand!(rng_type(arch), rnd.z)

    ##### compare the random number with the given probability (per time step or per hour whichever is shorter)
    ##### return 1 if random number is smaller
    @inbounds plank.graz .= isless.(rnd.x, plank.graz .* ΔT .* min(10,3600÷ΔT)) .* plank.ac
    @inbounds plank.mort .= isless.(rnd.y, plank.mort .* ΔT .* min(10,3600÷ΔT)) .* plank.ac
    @inbounds plank.dvid .= isless.(rnd.z, plank.dvid .* ΔT .* min(10,3600÷ΔT)) .* plank.ac
    return nothing
end
