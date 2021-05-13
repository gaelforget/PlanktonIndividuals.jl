##### update nutrient fields
function nut_update!(nutrients, Gcs, nut_temp, arch::Architecture, g::AbstractGrid, params, vel, consume, ΔT, iter)
    ##### compute advection tendency
    nut_advection!(nutrients, nut_temp, Gcs, vel, g, ΔT, arch)

    ##### compute nutrient diffusion,for each time step
    nut_diffusion!(Gcs, arch, g, nutrients, params["κh"], ΔT)

    ##### compute biogeochemical forcings of nutrients,for each time step
    zero_fields!(nut_temp)
    nut_forcing!(Gcs, nut_temp, nutrients, params, ΔT)

    ##### apply boundary conditions
    apply_bcs!(Gcs, nutrients, g, iter, arch)

    ##### apply diffusion and forcing tendency
    for name in nut_names
        nutrients[name].data .= nutrients[name].data .+ Gcs[name].data .* ΔT .+ consume[name].data ./ g.V
    end

    fill_halo_nut!(nutrients, g)
end

function nut_update!(nutrients, Gcs, nut_temp, arch::Architecture, g::AbstractGrid, params, consume, ΔT, iter)
    # compute biogeochemical forcings of nutrients,for each time step
    nut_forcing!(Gcs, nut_temp, nutrients, params, ΔT)

    ##### apply boundary conditions
    apply_bcs!(Gcs, nutrients, g, iter, arch)

    # apply forcing tendency
    for name in nut_names
        nutrients[name].data .= nutrients[name].data .+ Gcs[name].data .* ΔT .+ consume[name].data ./ g.V
    end

    fill_halo_nut!(nutrients, g)
end
