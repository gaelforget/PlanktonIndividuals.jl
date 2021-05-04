mutable struct PI_Model
    arch::Architecture          # architecture on which models will run
    t::Int64
    individuals::individuals    # initial individuals generated by `setup_agents`
    nutrients::NamedTuple       # initial nutrient fields
    grid::RegularRectilinearGrid# grid information
    bgc_params::Dict            # biogeochemical parameter set
    timestepper::timestepper    # operating Tuples and arrays for timestep
    diags::Diagnostics          # diagnostics
end

"""
    PI_Model(arch::Architecture, grid::RegularRectilinearGrid;
            individual_size = (Nsp = 1, N = 1024, cap = 8),
            bgc_params = bgc_params_default(), 
            phyt_params = phyt_params_default(),
            nut_source = default_nut_init(),
            diag_ntrs = (:PAR, :DOC, :NH4, :NO3),
            diag_nprocs = (:num, :graz, :mort, :dvid),
            t = 0.0,
            )
Generate the `PI_Model` struct on `grid`. 

Keyword Arguments
=================
- `arch` (required): `CPU()` or `GPU()`. The computer architecture used to time-step `model`.
- `grid` (required): The resolution and discrete geometry on which `model` is solved.
- `individual_size` (optional): `NamedTuple` used to set number of species `Nsp`, number of individuals `N`,
                                and max individual capacity `cap`.
- `bgc_params` (optional): Parameter set for biogeochemical processes modeled in the model.
- `phyt_params` (optional): Parameter set for physiological processes of individuals modeled in the model.
- `nut_source` (optional): The source of initial conditions of nutrient fields, should be either a `NamedTuple` 
                           or a `Dict` containing the file paths pointing to the files of nutrient initial conditions.
- `diag_ntrs` (optional): a `Tuple` containing the names of nutrient fields to be diagnosed.
- `diag_nprocs` (optional): a `Tuple` containing the names of physiological processes to be diagnosed.
- `t` (optional): Model time, start from 0 by default, in second.
- `mask` (optional): Mask out the individuals and tracers generated out of the domain, a 3D array with size `(Nx, Ny, Nz)`.
"""
function PI_Model(arch::Architecture, grid::RegularRectilinearGrid;
                  individual_size = (Nsp = 1, N = 1024, cap = 8),
                  bgc_params = bgc_params_default(), 
                  phyt_params = phyt_params_default(),
                  nut_source = default_nut_init(),
                  diag_ntrs = (:PAR, :DOC, :NH4, :NO3),
                  diag_nprocs = (:num, :graz, :mort, :dvid),
                  t = 0.0,
                  mask = nothing,
                  )

    if arch == GPU() && !has_cuda()
        throw(ArgumentError("Cannot create a GPU model. No CUDA-enabled GPU was detected!"))
    end

    inds = individuals(phyt_params, arch, individual_size.Nsp, individual_size.N, individual_size.cap)

    for plank in inds.phytos
        gen_individuals!(plank, individual_size.N, grid, arch; mask = mask)
    end

    nutrients = generate_nutrients(arch, grid, nut_source; mask = mask)

    ts = timestepper(arch, grid, individual_size.N, individual_size.cap)

    diags = diags_setup(diag_ntrs, diag_nprocs, grid, individual_size.Nsp, arch)

    model = PI_Model(arch, t, inds, nutrients, grid, bgc_params, ts, diags)

    return model
end

import Base: show

function show(io::IO, model::PI_Model)
    Nsp = length(model.individuals.phytos)
    N = Int(dot(model.individuals.phytos.sp1.data.ac,model.individuals.phytos.sp1.data.ac))
    cap = length(model.individuals.phytos.sp1.data.ac)
    print(io, "grid: Nx = $(model.grid.Nx), Ny = $(model.grid.Ny), Nz = $(model.grid.Nz)\n",
              "individuals: $(Nsp) phytoplankton species each with $(N) individuals\n",
              "capacity of individuals: $(cap) per species\n",
              "diagnostics of tracers: $(keys(model.diags.tr))\n",
              "diagnostics of individuals: $(keys(model.diags.spcs.sp1))")
end