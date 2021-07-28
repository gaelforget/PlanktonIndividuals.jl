using PlanktonIndividuals, Serialization

grid = RegularRectilinearGrid(size = (16, 16, 1), spacing = (2, 2, 2), halo = (2, 2, 2))

model = PlanktonModel(CPU(), grid) 

TP = sum((interior(model.nutrients.PO4.data, grid) .+ 
          interior(model.nutrients.DOP.data, grid) .+ 
          interior(model.nutrients.POP.data, grid)) .* grid.Δx .* grid.Δy .* grid.Δz)
TP = TP + sum(model.individuals.phytos.sp1.data.Pq .+ 
              model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC)

uvel = zeros(16,16,1,11)
vvel = zeros(16,16,1,11)
wvel = zeros(16,16,2,11)

for i in 1:11
    uvel[:,:,:,i] .= randn(16,16,1) .* 1e-4
    vvel[:,:,:,i] .= randn(16,16,1) .* 1e-4
    wvel[:,:,:,i] .= randn(16,16,2) .* 1e-4
end

sim = PlanktonSimulation(model, ΔT = 60, iterations = 10, vels=(u=uvel, v=vvel, w=wvel)) 

update!(sim)

TPt = sum((interior(model.nutrients.PO4.data, grid) .+ 
          interior(model.nutrients.DOP.data, grid) .+ 
          interior(model.nutrients.POP.data, grid)) .* grid.Δx .* grid.Δy .* grid.Δz)
TPt = TPt + sum(model.individuals.phytos.sp1.data.Pq .+ 
                model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC)
