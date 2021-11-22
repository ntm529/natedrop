function myCircleMat(droplet_radius, bottom_sheet_radius_in_droplets, top_sheet_radius_in_droplets, bottom_osm, top_osm, L)
# droplet_radius =  3 #droplet radius 
# droplet_osmolarity = 0.1
# bottom_sheet_radius_in_droplets = 1 # The number of droplets that make up the radius of this circular sheet in terms of droplets


droplet_locations_x = []
droplet_locations_y = []
droplet_locations_z = [] 
droplet_osmolarities = []





range = bottom_sheet_radius_in_droplets*droplet_radius
println("range is ", range)

x = collect(-range: droplet_radius*2: range)
display(x)
for i in 1:length(x)
    ψ = sqrt(abs(x[i]^2 - (bottom_sheet_radius_in_droplets*droplet_radius)^2))
    ψ = (round(ψ/(2*droplet_radius)))
    ψ = ψ * (2*droplet_radius)
    for j in 0:droplet_radius*2*L:ψ*L
        append!(droplet_locations_x, x[i]*L)
        append!(droplet_locations_y, j)
        append!(droplet_locations_z, 0)
        append!(droplet_osmolarities, bottom_osm)
    end
end


range = bottom_sheet_radius_in_droplets*droplet_radius- 2*droplet_radius 
println("range is ", range)
x = collect(-range: droplet_radius*2: range)
display(x)

for i in 1:length(x)
    ψ = sqrt(abs(x[i]^2 - (bottom_sheet_radius_in_droplets*droplet_radius)^2))
    ψ = (round(ψ/(2*droplet_radius)))
    ψ = ψ * (2*droplet_radius)
    for j in droplet_radius*2*L:droplet_radius*2*L:ψ*L
        append!(droplet_locations_x, x[i]*L)
        append!(droplet_locations_y, -j)
        append!(droplet_locations_z, 0)
        append!(droplet_osmolarities, bottom_osm)
    end
end


# <><><><><><><><><><> NOW GET THE TOP LAYER OF THE CIRCULAR MATRIX <><><><><><><><><><>
range = top_sheet_radius_in_droplets*droplet_radius
println("range is ", range)

x = collect(-range: droplet_radius*2: range)
display(x)
for i in 1:length(x)
    ψ = sqrt(abs(x[i]^2 - (top_sheet_radius_in_droplets*droplet_radius)^2))
    ψ = (round(ψ/(2*droplet_radius)))
    ψ = ψ * (2*droplet_radius)
    for j in 0:droplet_radius*2*L:ψ*L
        append!(droplet_locations_x, x[i]*L)
        append!(droplet_locations_y, j)
        append!(droplet_locations_z, norm([droplet_radius*L droplet_radius*L droplet_radius*L]))
        append!(droplet_osmolarities, top_osm)
    end
end

range = top_sheet_radius_in_droplets*droplet_radius- 2*droplet_radius 
println("range is ", range)
x = collect(-range: droplet_radius*2: range)
display(x)
for i in 1:length(x)
    ψ = sqrt(abs(x[i]^2 - (top_sheet_radius_in_droplets*droplet_radius)^2))
    ψ = (round(ψ/(2*droplet_radius)))
    ψ = ψ * (2*droplet_radius)
    for j in droplet_radius*2*L:droplet_radius*2*L:ψ*L
        append!(droplet_locations_x, x[i]*L)
        append!(droplet_locations_y, -j)
        append!(droplet_locations_z, norm([droplet_radius*L droplet_radius*L droplet_radius*L]))
        append!(droplet_osmolarities, top_osm)
    end
end

# ======================================= Finish it off =============================================== # 

velocities = zeros(length(droplet_locations_x),3)
osmolarities = bottom_osm * ones(length(droplet_locations_x))
osmolarities = top_osm * ones(length(droplet_locations_x))

radii = droplet_radius * ones(length(droplet_locations_x))

sheetmat = [radii droplet_osmolarities droplet_locations_x droplet_locations_y droplet_locations_z velocities]
sheetmat = Float64.(sheetmat)

return sheetmat

end # <<<<<<<< End of function
