module MyDropletVectorJ!
include("mydropletvectorjcomponents.jl")
include("myparticlepoints.jl")

using .MyDropletVectorJComponents
using LinearAlgebra
using DifferentialEquations
export dropvecj!
using .ParticlePoints

function dropvecj!(ρ, k, L , γ, D, timeend, ICmatrix, boolean_osm, tstepextract)

ICvector = reshape(ICmatrix, :, 1);
# display(ICvector)
# Define some of the model parameters 
dimensions = (size(ICmatrix,2)-2)/2;
dimensions = Int(dimensions);
dropletcount = Int(size(ICmatrix,1));
tspan = (0.0, timeend);

# The 3 value is the number of dimensions, x, y, z 
dradii = zeros(dropletcount,1)
dosmolarities = zeros(dropletcount,1)
dvelocities = zeros(3*dropletcount,1)
maxradius::Vector{Float64} = zeros(1)


# typeof(droplets)
# display(droplets)




# Create the current residual that makes the osmolarity activate after the droplets settle


p = [k, L, γ, dropletcount, D, boolean_osm, dradii, dosmolarities, dvelocities, ρ, maxradius] 
prob = ODEProblem(mydropfun!, ICvector, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, saveat = 0.1)
# graphics = mySplineConstructor(sol, dropletcount, dimensions, tstepextract)

W = myOvitoPrint!(sol, dropletcount, dimensions)
outfile = "E:/Research Scripts and Functions/OVITO_Files/visual_model.ovito"
f = open(outfile, "w")
for i in eachindex(W)   
    println(f, W[i])
end
close(f)



end # Function end


end # Module end