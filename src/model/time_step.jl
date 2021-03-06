"""
    PI_TimeStep!(model, ΔT, resultpath)
Update physiology processes and nutrient field of `PI_Model` one time step forward.

Keyword Arguments
=================
- `model`: `PI_Model` to be updated one time step forward
- `ΔT`: The length of a time step
- `resultpath` (optional): The file path to store model output. 
"""
function PI_TimeStep!(model::PI_Model, ΔT, resultspath::String)
    model.t = model.t+ΔT

    @inbounds model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5
    @inbounds model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5
    @inbounds model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5

    zero_fields!(model.timestepper.plk)
    @inbounds model.timestepper.chl .= 0.0
    @inbounds model.timestepper.pop .= 0.0  # may add an option of self grazing, besides shared grazing
    ##### plankton advection and diffusion
    for sp in keys(model.individuals.phytos)
        ##### RK4
        plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                               model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)
        ##### AB2 only for 1 species setup
        # plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos,
        #                     model.grid, 0.1, model.timestepper.vel₁, ΔT, model.arch)
        ##### Euler-forward
        # plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos,
        #                     model.grid, model.timestepper.vel₁, ΔT, model.arch)
        ##### Diffusion
        plankton_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd,
                            model.bgc_params["κhP"], ΔT, model.grid, model.arch)

        #### calculate accumulated chla quantity (not concentration)
        find_inds!(model.individuals.phytos[sp].data, model.grid, model.arch)
        acc_counts!(model.timestepper.chl, model.timestepper.pop,
                    model.individuals.phytos[sp].data.chl, model.individuals.phytos[sp].data.ac, 
                    model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi, 
                    model.individuals.phytos[sp].data.zi, model.arch)
    end

    ##### calculate PAR
    calc_par!(model.timestepper.par, model.arch, model.timestepper.chl, model.timestepper.PARF,
              model.grid, model.bgc_params["kc"], model.bgc_params["kw"])

    ##### plankton physiology
    for sp in keys(model.individuals.phytos)
        find_NPT!(model.timestepper.nuts, model.individuals.phytos[sp].data.xi, 
                  model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                  model.individuals.phytos[sp].data.ac, model.nutrients.NH4.data, 
                  model.nutrients.NO3.data, model.nutrients.PO4.data, model.nutrients.DOC.data,
                  model.timestepper.par, model.timestepper.temp, model.timestepper.pop,
                  model.individuals.phytos[sp].p)

        plankton_physiology!(model.individuals.phytos[sp].data, model.timestepper.nuts,
                             model.individuals.phytos[sp].proc, model.individuals.phytos[sp].p,
                             model.timestepper.plk, model.diags.spcs[sp], ΔT, model.t, model.arch)
    end
    write_species_dynamics(model.t, model.individuals.phytos, resultspath)

    ##### diagnostics for nutrients
    @inbounds model.diags.tr.PAR .+= model.timestepper.par
    for key in keys(model.diags.tr)
        if key in keys(model.nutrients)
            @inbounds model.diags.tr[key] .+= model.nutrients[key].data
        end
    end

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.nut_temp, model.arch,
                model.grid, model.bgc_params, model.timestepper.vel₁, model.timestepper.plk, ΔT)

    write_nut_cons(model.grid, model.timestepper.Gcs, model.nutrients, model.t, resultspath)

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data

    return nothing
end

function PI_TimeStep!(model::PI_Model, ΔT)
    model.t = model.t+ΔT

    @inbounds model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5
    @inbounds model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5
    @inbounds model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5

    zero_fields!(model.timestepper.plk)
    @inbounds model.timestepper.chl .= 0.0
    @inbounds model.timestepper.pop .= 0.0  # may add an option of self grazing, besides shared grazing

    ##### plankton advection
    for sp in keys(model.individuals.phytos)
        ##### RK4
        plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos, model.grid,
                               model.timestepper.vel₀, model.timestepper.vel½, model.timestepper.vel₁, ΔT, model.arch)
        ##### AB2 only for 1 species setup
        # plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos,
        #                     model.grid, 0.1, model.timestepper.vel₁, ΔT, model.arch)
        ##### Euler-forward
        # plankton_advection!(model.individuals.phytos[sp].data, model.timestepper.velos,
        #                     model.grid, model.timestepper.vel₁, ΔT, model.arch)
        ##### Diffusion
        plankton_diffusion!(model.individuals.phytos[sp].data, model.timestepper.rnd,
                            model.bgc_params["κhP"], ΔT, model.grid, model.arch)

        #### calculate accumulated chla quantity (not concentration)
        find_inds!(model.individuals.phytos[sp].data, model.grid, model.arch)
        acc_counts!(model.timestepper.chl, model.timestepper.pop,
                    model.individuals.phytos[sp].data.chl, model.individuals.phytos[sp].data.ac, 
                    model.individuals.phytos[sp].data.xi, model.individuals.phytos[sp].data.yi, 
                    model.individuals.phytos[sp].data.zi, model.arch)
    end

    ##### calculate PAR
    calc_par!(model.timestepper.par, model.arch, model.timestepper.chl, model.timestepper.PARF,
              model.grid, model.bgc_params["kc"], model.bgc_params["kw"])

    ##### plankton physiology
    for sp in keys(model.individuals.phytos)
        find_NPT!(model.timestepper.nuts, model.individuals.phytos[sp].data.xi, 
                  model.individuals.phytos[sp].data.yi, model.individuals.phytos[sp].data.zi, 
                  model.individuals.phytos[sp].data.ac, model.nutrients.NH4.data, 
                  model.nutrients.NO3.data, model.nutrients.PO4.data, model.nutrients.DOC.data,
                  model.timestepper.par, model.timestepper.temp, model.timestepper.pop,
                  model.individuals.phytos[sp].p)

        plankton_physiology!(model.individuals.phytos[sp].data, model.timestepper.nuts,
                             model.individuals.phytos[sp].proc, model.individuals.phytos[sp].p,
                             model.timestepper.plk, model.diags.spcs[sp], ΔT, model.t, model.arch)
    end

    nut_update!(model.nutrients, model.timestepper.Gcs, model.timestepper.nut_temp, model.arch,
                model.grid, model.bgc_params, model.timestepper.vel₁, model.timestepper.plk, ΔT)

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data

    return nothing
end