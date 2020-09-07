module PlanktonIndividuals

using NCDatasets, Serialization
using Random, Distributions, Statistics
using Printf
using CUDA, KernelAbstractions

src=""

include("$src"*"architecture.jl")
include("$src"*"model/grids.jl")
include("$src"*"params/param_default.jl")
include("$src"*"params/param_update.jl")
include("$src"*"nut/fields.jl")
include("$src"*"nut/diffusivity.jl")
include("$src"*"nut/forcing.jl")
include("$src"*"nut/third_order_DSTFL.jl")
include("$src"*"nut/multi_dim_adv.jl")
include("$src"*"nut/halo_regions.jl")
include("$src"*"nut/gen_nut_fields.jl")
include("$src"*"nut/nutrient_processes.jl")
include("$src"*"plankton/gen_plankton.jl")
include("$src"*"plankton/advection/velocity_interpolations.jl")
include("$src"*"plankton/advection/advection_operations.jl")
include("$src"*"plankton/advection/plankton_advection.jl")
include("$src"*"plankton/advection/plankton_diffusion.jl")
include("$src"*"plankton/phyt_process.jl")
include("$src"*"plankton/zooplankton.jl")
include("$src"*"utils.jl")
include("$src"*"output/output_writers.jl")
include("$src"*"output/diagnostics.jl")
include("$src"*"model/models.jl")
include("$src"*"model/time_step.jl")


export
    # model structures
    PI_Model, Grids, nutrient_fields, velocities,
    RunOptions, RunParams, read_Ogrids, gen_Grid,
    Architecture, GPUs, CPUs,

    # read input functions
    read_IR_input, read_temp_input,
    update_params!, param_default,
    PrepRunDir, generate_vel_itp, diags_setup,

    # initialize nutrient field and individual sets
    gen_agents, gen_nutrients, load_nut_initials,

    # Run the model
    RunParam, RunOption, PI_TimeStep!,
    PI_advectRK4!, PI_advect!,

    # write output functions
    write_nut_nc_alltime, write_nut_nc_each_step,
    count_vertical_num, count_horizontal_num,
    write_species_dynamics, write_nut_cons,
    sort_species, write_output
end # module
