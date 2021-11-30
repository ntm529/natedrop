#= Let's create a matrix that has two droplets next to each other on the x-axis. These droplets should oscillate for a bit, and then diffusion will occur when velocity reaches zero. 
Additionally, let the droplets be initially at rest when this simulation begins: =#
include("E:/Research Scripts and Functions/Julia Scripts/mydropletvectorjcomponents.jl")
using .MyDropletVectorJComponents
using LinearAlgebra
using DifferentialEquations
L = 0.8 # The natural length proportional to the radii that droplets will settle at
k = 1000. # The spring force between each droplet
rho = 0.2 # The density of the droplet, used to calculate the mass of each droplet within the system as a function of radius.
gamma = 30. # The simplified version of both fluid response and spring damping in the system
D = 2. *10^-2  # The rate at which diffusion will happen in the system
boolean_osm = [0.0] # The "on" switch for when the diffusion in the sytem will activate
timeend = 250. # The time at which the ODE will end

# Now create the matrix for the initial radius, osmolarity, positions, and velocities for these two droplets on the x-axis
ICmatrix = [3 0.1 0 0 0 0 0 0; 3 1 6 0 0 0 0 0]

# Now call the function dropvecj!, which will solve an ordinary differential equation for the given timespan, and will print out an xyz file 
dropvecj!(rho, k, L , gamma, D, timeend, ICmatrix, boolean_osm)

