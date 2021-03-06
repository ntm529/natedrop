module XX

include("myLoopMatrix.jl")
include("myStableLoopMatrix.jl")


include("mysheetmatrix.jl")

# include("E:/Research Scripts and Functions/Julia Scripts/mydropletvectorj!.jl")
include("E:/Research Scripts and Functions/Julia Scripts/mydropletvectorjcomponents.jl")
include("Droplet_Sphere_sheet.jl")

# using .MyDropletVectorJ!
using .MyDropletVectorJComponents
using LinearAlgebra
using DifferentialEquations
using Plots
using Interpolations
# using DataFrames
using StatsPlots

#Characters I might use for some initialized matrixes
# ℵ ℶ ℷ ℸ
# For a more "theatrical" effect which exaggerates the momentum and gives the visualizations more charisma, set ρ and k higher and set γ lower. 

# volume is given as 10
radii = 3
# (3/4 * 1/pi * 10)^(1/3) for the real radii 
dropletcount =  5
osmolaritybase = 1.
osmolaritytop = 0.1


# Note: For future reference, the method should be used on Float64 values for type specification, right now Any is used.
ρ = 0.2 # The density of the droplet, used to calculate the mass of each droplet within the system as a function of radius. To wit Eqn = mass(rho,volume(radius))
# 0.2/10 for the real ρ
k = 1000. # The spring force between each droplet 
L = 0.8 # The natural length proportional to the radii that droplets will settle at
γ = 30. # The simplified version of both fluid response and spring damping in the system
timeend = 250. # The time at which the ODE will end
D = 2. *10^-2  # The rate at which diffusion will happen in the system
boolean_osm = [0.0] # The "on" switch for when the diffusion in the sytem will activate
tstepextract = Int(timeend * 10) 





# ICmatrix = myLoopMatrix(radii , dropletcount, osmolaritybase, osmolaritytop); # This is another u0 for testing
# function myCircleMat(droplet_radius, bottom_sheet_radius_in_droplets, top_sheet_radius_in_droplets, bottom_osm, top_osm, L)

# ICmatrix = myCircleMat(radii, 6, 5, osmolaritybase, osmolaritytop, L)


# ICmatrix = myStableLoopMatrix(radii, L, dropletcount, osmolaritybase, osmolaritytop); # This is another u0 for testing

# ICmatrix = [3 0.1 0 0 0 0 0 0; 3 1 0 0 6 0 0 0] # The incorrect z coords 
# ICmatrix = [3 0.1 0 0 0 0 0 0; 3 1 0 6 0 0 0 0] # The y coords



# ICmatrix = [3 0.1 0 0 0 0 0 0; 3 1 6 0 0 0 0 0] # The x coords




# ICmatrix = [3 0.1 0 0 0 0 0 0; 3 1 4 3 2 0 0 0; 3 1 -4 1 1 0 0 0] # The hybrid coords
# a = @elapsed dropvecj!(ρ, k, L , γ, D, timeend, ICmatrix, boolean_osm, tstepextract)
# ICmatrix = myStableSheetMatrix(radii, dropletcount, L)
# display(ICmatrix)
# a = @time dropvecj!(ρ, k, L , γ, D, timeend, ICmatrix, boolean_osm)

# dropvecj!(ρ, k, L , γ, D, timeend, ICmatrix, boolean_osm)


O = []
drops = []
for i in 22:52

    ICmatrix = myStableSheetMatrix(radii, i, L)
    push!(drops, (i^2 + (i-1)^2))
    display(drops)
    a = @elapsed dropvecj!(ρ, k, L , γ, D, timeend, ICmatrix, boolean_osm)
    push!(O,a)
    display(O)
    zed = plot(drops, O, title = "Test results vs nlog(ϕ,n) time" , bg = :black, legend =:topleft, label = "Results")
    plot!(drops, drops .* log.(100,drops), label = "nLog_{100}n complexity")

    display(zed)
    

end # end loop


# if dropletcount > 4    
#     xinner = graphics[1]
#     yinner = graphics[2]
#     xsinner = graphics[3]
#     ysinner = graphics[4]
#     xouter = graphics[5]
#     youter = graphics[6]
#     xsouter = graphics[7]
#     ysouter = graphics[8]
#     xcenter = graphics[9]
#     ycenter = graphics[10]
#     xscenter = graphics[11]
#     yscenter = graphics[12]




#     # display(scatter(xinner, yinner, label="outer knots", legend=:topleft))
#     display(plot(xsinner, ysinner, label="outer spline", legend=:topleft))

#     # display(plot!(xsinner, ysinner, label="outer spline"))
#     # display(scatter!(xouter, youter, label="inner knots"))
#     display(plot!(xsouter, ysouter, label="inner spline"))
#     display(plot!(xscenter, yscenter, label="center spline"))

# end




println("done")



end
