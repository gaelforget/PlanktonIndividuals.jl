using PlanktonIndividuals, Serialization

grid = RectilinearGrid(size = (1, 1, 1), x = (0,32), y = (0,32), z = (0,-32))

model = PlanktonModel(CPU(), grid;
                      N_species = 5,
                      N_individual = 1024,
                      max_individuals = 1024*10)

function tot_mass(nut, g)
    mass = zeros(g.Nx, g.Ny, g.Nz)
    for i in 1:g.Nx
        for j in 1:g.Ny
            for k in 1:g.Nz
                mass[i,j,k] = nut[i+g.Hx, j+g.Hy, k+g.Hz] * PlanktonIndividuals.Grids.volume(i+g.Hx, j+g.Hy, k+g.Hz, g)
            end
        end
    end
    return sum(mass)
end

TP = tot_mass(model.nutrients.PO4.data, grid) +
     tot_mass(model.nutrients.DOP.data, grid) +
     tot_mass(model.nutrients.POP.data, grid)
TP = TP + sum(model.individuals.phytos.sp1.data.Pq .+ 
              model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC) + 
          sum(model.individuals.phytos.sp2.data.Pq .+ 
              model.individuals.phytos.sp2.data.Bm .* model.individuals.phytos.sp2.p.R_PC) +
          sum(model.individuals.phytos.sp3.data.Pq .+ 
              model.individuals.phytos.sp3.data.Bm .* model.individuals.phytos.sp3.p.R_PC) +
          sum(model.individuals.phytos.sp4.data.Pq .+ 
              model.individuals.phytos.sp4.data.Bm .* model.individuals.phytos.sp4.p.R_PC) +
          sum(model.individuals.phytos.sp5.data.Pq .+ 
              model.individuals.phytos.sp5.data.Bm .* model.individuals.phytos.sp5.p.R_PC)


sim = PlanktonSimulation(model, ΔT = 60, iterations = 10)

update!(sim)

TPt = tot_mass(model.nutrients.PO4.data, grid) +
      tot_mass(model.nutrients.DOP.data, grid) +
      tot_mass(model.nutrients.POP.data, grid)
TPt=TPt + sum(model.individuals.phytos.sp1.data.Pq .+ 
              model.individuals.phytos.sp1.data.Bm .* model.individuals.phytos.sp1.p.R_PC) +
          sum(model.individuals.phytos.sp2.data.Pq .+ 
              model.individuals.phytos.sp2.data.Bm .* model.individuals.phytos.sp2.p.R_PC) +
          sum(model.individuals.phytos.sp3.data.Pq .+ 
              model.individuals.phytos.sp3.data.Bm .* model.individuals.phytos.sp3.p.R_PC) +
          sum(model.individuals.phytos.sp4.data.Pq .+ 
              model.individuals.phytos.sp4.data.Bm .* model.individuals.phytos.sp4.p.R_PC) +
          sum(model.individuals.phytos.sp5.data.Pq .+ 
              model.individuals.phytos.sp5.data.Bm .* model.individuals.phytos.sp5.p.R_PC)
