mutable struct Model_Input
    temp::AbstractArray{Float64,4}      # temperature
    PARF::AbstractArray{Float64,3}      # PARF
end

mutable struct Model_Struct
    arch::Architecture          # architecture on which models will run
    t::Int64
    individuals::individuals    # initial individuals generated by `setup_agents`
    nutrients::NamedTuple       # initial nutrient fields
    grid::Grids                 # grid information
    input::Model_Input          # model input, temp and PAR
    params::Dict                # biogeochemical parameter set
    diags::Diagnostics          # diagnostics
    timestepper::timestepper    # operating Tuples and arrays for timestep
end

"""
    PI_model(grid, RunParam)
Generate the model structure for time step
Default distribution of individuals is Normal distribution with 1.0 as mean and 0.25 as SD
Default PAR and temp are from ../samples
"""
function PI_Model(arch::Architecture, grid, RunParam;
                  t = 1-RunParam.ΔT,
                  individuals = individuals(RunParam.params, arch),
                  nutrients,
                  PARF = read_IR_input(RunParam.ΔT, grid),
                  temp = read_temp_input(RunParam.ΔT, grid),
                  params = RunParam.params,
                  diags = diags_setup(arch, RunParam.nTime, RunParam.ΔT, grid,
                                     RunParam.params["diag_freq"], RunParam.params["Nsp"], 5),
                  )

    if arch == GPUs() && !has_cuda()
        throw(ArgumentError("Cannot create a GPU model. No CUDA-enabled GPU was detected!"))
    end

    if arch == GPUs()
        input = Model_Input(CuArray(temp), CuArray(PARF))
    else
        input = Model_Input(temp,PARF)
    end

    for plank in individuals.phytos
        gen_individuals!(plank, params["Nind"], grid, arch)
    end

    ts = timestepper(arch, grid, params["Nind"])

    model = Model_Struct(arch, t, individuals, nutrients, grid, input, params, diags, ts)

    if RunParam.Zoo == true
        model.params["grz_P"] = 0
    else
        nothing
    end
    return model
end
