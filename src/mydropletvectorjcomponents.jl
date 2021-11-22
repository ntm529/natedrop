module MyDropletVectorJComponents
include("myparticlepoints.jl")
# include("myquadtreequery.jl")
using Interpolations
using LinearAlgebra


using RegionTrees
using StaticArrays: SVector
import RegionTrees: AbstractRefinery, needs_refinement, refine_data
using .ParticlePoints
using LinearAlgebra

# Continue working through this 
using RegionTrees
using StaticArrays: SVector
import RegionTrees: AbstractRefinery, needs_refinement, refine_data, adaptivesampling!


# The two major files that need to imported that were written out by me 
using .ParticlePoints
# using .QuadTreeQuery

# import RegionTrees: AbstractRefinery, needs_refinement, refine_data


# The exports of this particular module 
export mydropfun!, myOvitoPrint!


# using StaticArrays: SVector




function mydropfun!(du, u, p, t)
    # Put in the parameters that are needed, very important line of code: #=========================================================================
    k, L, γ, dropletcount, D, boolean_osm, dradii, dosmolarities, dvelocities, ρ = p;
    # ================================================================================================ # 
    # println("NEW FRAME # -------------------------------------------------------------------------------------------------------------- #")
    # display(radii)
    # printstyled("The maximum radius type  is:   " ; color= :green)\
    # printstyled("The current position is:" ; color = :magenta)
    # position = u[Int(2*dropletcount+1) : Int(5*dropletcount) ]
    # display(position)
    maxradius = maximum(view(u, 1:dropletcount))
    
    # display(typeof(maxradius))
    # display(radii)
    # display(osmolarities)
    # display(positions)
    # printstyled("The velocities being put in are: :" ; color = :magenta)
       
    # display(dposition)

    # display(velocities)
    normvelocity = norm(view(u, 5*dropletcount+1:8*dropletcount))

    # printstyled("The type of  SVECTOR being put in are: :" ; color = :magenta)
    # typeof(SVector)
    # display(SVector)

    if (t > 0.1 && normvelocity/dropletcount <= 1e-4)  
        boolean_osm[1] = 1.0
    end # End conditional

    droplets::Vector{Droplet} = Droplet[] # Change this to a vector that needs to be imported in 
    #Create the requisite list of droplets that will then be queried on: 
    for i in 1:dropletcount
        thisradius::Float64 = u[i]
        # printstyled("The this radius type  is:   " ; color= :green)
        # display(typeof(thisradius)) 
        thisosmolarity::Float64 = u[Int(i+dropletcount)]
        thisx::Float64 = u[Int(i+2*dropletcount)]
        thisy::Float64 = u[Int(i+3*dropletcount)]
        thisz::Float64 = u[Int(i+4*dropletcount)]
        thisxvelocity::Float64 =  u[Int(i+5*dropletcount)]
        thisyvelocity::Float64 =  u[Int(i+6*dropletcount)]
        thiszvelocity::Float64 =  u[Int(i+7*dropletcount)]
        droplet::Droplet = Droplet(thisradius, thisosmolarity, thisx, thisy, thisz, thisxvelocity, thisyvelocity, thiszvelocity, 0, 0, 0, maxradius + thisradius)
        push!(droplets, droplet)
    end # End the for loop

    


  
    xvec = view(u, Int(2*dropletcount+1):Int(3*dropletcount))
    xorigin::Float64 = minimum(xvec) - 0.1
    yvec = view(u, Int(3*dropletcount+1):Int(4*dropletcount))
    yorigin::Float64 = minimum(yvec) - 0.1
    zvec = view(u, Int(4*dropletcount+1):Int(5*dropletcount))
    zorigin::Float64 = minimum(zvec) - 0.1

    xwidth::Float64 = (maximum(xvec) - xorigin) * 1.1
    ywidth::Float64 = (maximum(yvec) - yorigin) * 1.1
    zwidth::Float64 = (maximum(zvec) - zorigin) * 1.1

    



    # printstyled("xorigin is:   " ; color= :green)
    # display(xorigin)
    # printstyled("yorigin is:   " ; color= :green)
    # display(yorigin)
    # printstyled("zorigin is:   " ; color= :green)
    # display(zorigin)

    # printstyled("xwidth is:   " ; color= :green)
    # display(xwidth)
    # printstyled("ywidth is:   " ; color= :green)
    # display(ywidth)
    # printstyled("zwidth is:   " ; color= :green)
    # display(zwidth)

    # ===================================================== BEGIN AND END THE COMMENT BLOCK HERE ============================================================== # 
    r = MyRefinery(1.0)
    root::Cell = Cell(SVector(xorigin, yorigin, zorigin), SVector(xwidth, ywidth, zwidth), droplets)
    myquadtree = adaptivesampling!(root, r)
    # println("the dradii at the start is:")
    # display(dradii)
    


    xyzaccel = zeros(3,1)
    

    for i in 1:length(droplets)
        found::Vector{Droplet} = Droplet[]
        # display(i)
        founddrops::Vector{Droplet} = query(myquadtree, droplets[i], found)
        # printstyled("the droplet being examined is: $(droplets[i])  " ; color= :green)
        # display(founddrops)
        myJ = 0 
        for j in 1:length(founddrops)
            # First deal with the update to acceleration
            # Oh, I also need to add mass to this later... 
            xdist = founddrops[j].x - droplets[i].x
            
            xyzaccel[1] = founddrops[j].x - droplets[i].x
            xyzaccel[2] = founddrops[j].y - droplets[i].y
            xyzaccel[3] = founddrops[j].z - droplets[i].z
            # printstyled("The xdist is going to be given as such:" ; color = :blue)
            # display(xdist)
            ydist = founddrops[j].y - droplets[i].y
            zdist = founddrops[j].z - droplets[i].z
            # printstyled("The zdist is going to be given as such:" ; color = :green)
            # display(zdist)

            directionx = xyzaccel[1]/sqrt(xyzaccel[1]^2 + xyzaccel[2]^2 + xyzaccel[3]^2)
            

            # printstyled("The direction X is going to be given as such:" ; color = :green)
            # display(directionx)
            directiony = xyzaccel[2]/sqrt(xyzaccel[1]^2 + xyzaccel[2]^2 + xyzaccel[3]^2)
            directionz = xyzaccel[3]/sqrt(xyzaccel[1]^2 + xyzaccel[2]^2 + xyzaccel[3]^2)
            # printstyled("The direction Z is going to be given as such:" ; color = :green)
            # display(directionz)
            mybilayerdist = L * (droplets[i].radius + founddrops[j].radius)

            #myaccelx# 
            # Theoretical ugly code implementation of a full line: (founddrops[j].x - droplets[i].x)/sqrt((founddrops[j].x - droplets[i].x)^2 + (founddrops[j].y - droplets[i].y)^2 + (founddrops[j].z - droplets[i].z)^2) * -k * (abs((L * (droplets[i].radius + founddrops[j].radius) ) * (founddrops[j].x - droplets[i].x)/sqrt((founddrops[j].x - droplets[i].x)^2 + (founddrops[j].y - droplets[i].y)^2 + (founddrops[j].z - droplets[i].z)^2))  -abs(founddrops[j].x - droplets[i].x) )
            
            xyzaccel[1] = directionx * -k * (abs(mybilayerdist * directionx) - abs(xyzaccel[1]))
            xyzaccel[2] = directiony * -k * (abs(mybilayerdist * directiony) - abs(xyzaccel[2])) 
            xyzaccel[3] = directionz * -k * (abs(mybilayerdist * directionz) - abs(xyzaccel[3])) 
            droplets[i].xacceleration += xyzaccel[1]
            droplets[i].yacceleration += xyzaccel[2]
            droplets[i].zacceleration += xyzaccel[3]


            # xyzaccel[1] =                            -   abs(founddrops[j].x - droplets[i].x)
            # The distance part 
            #x# (founddrops[j].x - droplets[i].x)
            #y# (founddrops[j].y - droplets[i].y)
            #z# (founddrops[j].z - droplets[i].z)
            # The directional part    
            #x# (founddrops[j].x - droplets[i].x)/sqrt((founddrops[j].x - droplets[i].x)^2 + (founddrops[j].y - droplets[i].y)^2 + (founddrops[j].z - droplets[i].z)^2)
            #y# (founddrops[j].y - droplets[i].y)/sqrt((founddrops[j].x - droplets[i].x)^2 + (founddrops[j].y - droplets[i].y)^2 + (founddrops[j].z - droplets[i].z)^2)
            #z# (founddrops[j].z - droplets[i].z)/sqrt((founddrops[j].x - droplets[i].x)^2 + (founddrops[j].y - droplets[i].y)^2 + (founddrops[j].z - droplets[i].z)^2)
            # The bilayer distance part 
            ## (L * (droplets[i].radius + founddrops[j].radius))


            # =================== The osmolarity part of the code ================================= # 
            intersectcase1 =  droplets[i].radius + founddrops[j].radius - sqrt(xdist^2 + ydist^2 + zdist^2) >= 0
            intersectcase2 = sqrt(xdist^2 + ydist^2 + zdist^2) + droplets[i].radius  - founddrops[j].radius >= 0
            intersectcase3 = sqrt(xdist^2 + ydist^2 + zdist^2) + founddrops[j].radius  - droplets[i].radius >= 0

            if boolean_osm[1] >= 1.0 
                if (intersectcase1 && intersectcase2 && intersectcase3)
                    component1 = (sqrt(xdist^2 + ydist^2 + zdist^2)^2 - founddrops[j].radius^2 + droplets[i].radius^2)^2
                    component2 = 4 * sqrt(xdist^2 + ydist^2 + zdist^2)^2 * droplets[i].radius^2
                    component3 = component2 - component1
                    component4 = 1/(2*sqrt(xdist^2 + ydist^2 + zdist^2)) * sqrt(component3)
                    mycircular_area = pi * component4^2
                    osmdiff = droplets[i].osmolarity - founddrops[j].osmolarity
                    myJ += D * mycircular_area * osmdiff
                end # End the conditional
            end # End the conditional
        end # End the for loop

        myvolume = 4/3 * pi * droplets[i].radius^3
        dvelocities[i] = (droplets[i].xacceleration - γ * droplets[i].xvelocity) / (ρ * myvolume) 
        dvelocities[i+dropletcount] = (droplets[i].yacceleration - γ * droplets[i].yvelocity) / (ρ * myvolume) 
        dvelocities[i+2*dropletcount] = (droplets[i].zacceleration - γ * droplets[i].zvelocity) / (ρ * myvolume) 

        mymoles = droplets[i].osmolarity * myvolume
        myvolume = myvolume + myJ 
        dradii[i] = (3/(4*pi) * myvolume)^(1/3) - droplets[i].radius
        dosmolarities[i] = mymoles/myvolume - droplets[i].osmolarity
        
    end # End the for loop

    # ================================ BEGIN AND END THE COMMENT BLOCK HERE =================================================== # 
    # printstyled("The dposition is going to be given as such:" ; color = :cyan)
    # display(dposition)



    
    # printstyled("The dvelocities are going to be given as such:" ; color = :green)
    # display(dvelocities)
    # printstyled("So the change in position will be something like:" ; color = :green)
    # display(dvelocities .+ dposition)

    
    du[1:dropletcount] = dradii
    du[Int(dropletcount+1):Int(2*dropletcount)] = dosmolarities
    du[Int(2*dropletcount+1) : Int(5*dropletcount) ] = view(u, 5*dropletcount+1:8*dropletcount)
    du[Int(5*dropletcount + 1) : Int(8*dropletcount)] = dvelocities
    
    # ============================================ # The current non-essential functions that are fixing the beam in place
    # ℵ = dropletcount+1 # the variable used to fix the first droplet in space
    # ℶ = dropletcount + 1 + dropletcount # the variable used to fix the second droplet in space
    # ℵ = 2*dropletcount+1
    # ℶ = 2*dropletcount + 1 + dropletcount/2
    # du[Int(ℵ)] = 0.
    # du[Int(ℵ+dropletcount)] = 0.
    # du[Int(ℵ+2*dropletcount)] = 0.
    # du[Int(ℶ)] = 0.
    # du[Int(ℶ+dropletcount)] = 0.
    # du[Int(ℶ+2*dropletcount)] = 0.
    end # End the main function
# ===========================================================================================================================================================================



function myOvitoPrint!(sol, dropletcount, dimensions)

    # This function will recast the solution of our ordinary differential equation into a printed file that can be visualized in ovito.
    
    V = convert(Array,sol) # Convert the ODE solution into an array that can be worked on 
    V = dropdims(V, dims = 2) # Dimension 2 is a singleton and should be eliminated
    # display(V)
    V = V'
    tsteps = size(V,1)
    omicron = .~I(dimensions);
    omicron = repeat(omicron, dropletcount, 1)
    omicron = reshape(omicron, dimensions, :)
    omicron = omicron'
    omicron = repeat(omicron, tsteps, 1)
    omicron = reshape(omicron, dimensions*dropletcount, :)
    
    # Shave off the velocities for now 
    W = V[:,2*dropletcount+1:2*dropletcount + dimensions*dropletcount]; ##ADJUST##
    # Perform a sequence of transformations to re-arrange W to what we want
    Wtrue = [W[:,1:2]; W[:,3:4] ; W[:,5:6]]
    Wtransptrue = Wtrue'
    W = W'
    W = repeat(W, 1, dimensions)
    W = W[.~omicron]
    W = reshape(W, dropletcount, :)
    W = W'
    W = reshape(W, tsteps, :)
    W = W'
    W = reshape(W, dimensions, :)
    W = W'
    
    radvector = reshape(V[:,1:dropletcount]',1,:)';
    osmvector = reshape(V[:,dropletcount+1:2*dropletcount]',1,:)'
    
    W = [W radvector osmvector]
    # Now create a unified string vector 
    α = []
    ω = size(W,1)
    
    W = string.(W)
    for i = 1:ω
        push!(α,join(W[i,:]," "))
    end
    timetitle = ["ITEM: TIMESTEP"];
    timetitle = repeat(timetitle, Int(size(W,1)/dropletcount), 1);
    
    timenumber = convert(Array, sol.t)
    timenumber = string.(timenumber)
    
    dropnumtitle = ["ITEM: NUMBER OF ATOMS"]
    dropnumtitle = repeat(dropnumtitle,Int(size(W,1)/dropletcount),1 )
    
    dropnumnumber = [string.(dropletcount)]
    dropnumnumber = repeat(dropnumnumber, Int(size(W,1)/dropletcount), 1)
    
    
    typetitle = ["ITEM: ATOMS x y z radius osmolarity"]
    typetitle = repeat(typetitle, Int(size(W,1)/dropletcount), 1)
    
    assembly = permutedims([timetitle timenumber dropnumtitle dropnumnumber typetitle])
    assembly = reshape(assembly, :, 1);
    # timenumber dropnumtitle dropnumnumber typetitle]
    
    boolean_me = repeat([["boolean_me"]; repeat(["boolean_me2"],Int(dropletcount-1) ,1)],tsteps,5);
    W = α
    W = [W boolean_me]
    W = permutedims(W) # Seems good to this point
    W = reverse(W, dims = 1)
    W = permutedims(reshape(W,1,:))
    jj = W .== "boolean_me2"
    W= W[.~jj]
    boolean_me = [repeat(["boolean_me2"],4,1); ["boolean_me"]];
    assembly = [assembly repeat(boolean_me,tsteps,dropletcount)]
    assembly = permutedims(assembly)
    assembly = reshape(assembly,1,:)
    assembly = permutedims(assembly)
    
    kk = assembly .== "boolean_me2"
    assembly = assembly[.~kk]
    W = [W assembly]
    W = permutedims(W)
    W = reshape(W,1,:)
    ii = W .== "boolean_me"
    W = W[.~ii]
    
    return W
    
    end # Function end 







# =============================================== All the functions that really should go inside their own module ========================# 
#=========================================================================================================================================#
# =============================================== The refinement structure that creates our quad tree ====================================# 
struct MyRefinery <: AbstractRefinery
    tolerance::Float64
end # End struct


function needs_refinement(r::MyRefinery, cell)
    pointvector::Vector{Droplet} = cell.data
    # display(typeof(pointvector))
#     display(pointvector)
    list::Vector{Droplet} = Droplet[]
    for i in 1:length(pointvector)
        # Adress an edge case where the point sits directly in top of the line of division
        edgecorrectionx = 0. # Intialize this as zero every loop and then make it "budge up or down" the tree structure if it is on an edge of the cell
        edgecorrectiony = 0. # Intialize this as zero every loop and t|hen make it "budge up or down" the tree structure if it is on an edge of the cell
        edgecorrectionz = 0. # Intialize this as zero every loop and t|hen make it "budge up or down" the tree structure if it is on an edge of the cell
        
        edgecasex1::Bool = pointvector[i].x == cell.boundary.origin[1] 
        edgecasex2::Bool = pointvector[i].x == cell.boundary.origin[1] + cell.boundary.widths[1] 
        edgecasey1::Bool = pointvector[i].y == cell.boundary.origin[2] 
        edgecasey2::Bool = pointvector[i].y == cell.boundary.origin[2] + cell.boundary.widths[2] 
        edgecasez1::Bool = pointvector[i].z == cell.boundary.origin[3] 
        edgecasez2::Bool = pointvector[i].z == cell.boundary.origin[3] + cell.boundary.widths[3] 
#         println("Cell.boundary.origin 3 is $(cell.boundary.origin[3])")
#         println("Cell.boundary.widths 3 is $(cell.boundary.widths[3])")   
        
        edgecasex1 ? edgecorrectionx::Float64 += (cell.boundary.widths[1]+10)^(-1/cell.boundary.widths[1]) : nothing 
        edgecasex2 ? edgecorrectionx::Float64 += (cell.boundary.widths[1]+10)^(-1/cell.boundary.widths[1]) : nothing
        edgecasey1 ? edgecorrectiony::Float64 += (cell.boundary.widths[2]+10)^(-1/cell.boundary.widths[2]) : nothing
        edgecasey2 ? edgecorrectiony::Float64 += (cell.boundary.widths[2]+10)^(-1/cell.boundary.widths[2]) : nothing
        edgecasez1 ? edgecorrectionz::Float64 += (cell.boundary.widths[3]+10)^(-1/cell.boundary.widths[3]) : nothing
        edgecasez2 ? edgecorrectionz::Float64 += (cell.boundary.widths[3]+10)^(-1/cell.boundary.widths[3]) : nothing
#         println("Not getting an error yet")
      
        case1 = cell.boundary.origin[1] < pointvector[i].x + edgecorrectionx < cell.boundary.origin[1] + cell.boundary.widths[1]  
        case2 = cell.boundary.origin[2] < pointvector[i].y + edgecorrectiony < cell.boundary.origin[2] + cell.boundary.widths[2] 
        case3 = cell.boundary.origin[3] < pointvector[i].z + edgecorrectionz < cell.boundary.origin[3] + cell.boundary.widths[3] 
        case1 && case2 && case3 ? push!(list,pointvector[i]) : nothing
        
        cell.data = list
#         display(list)
    end # End the for loop
    length(list) > r.tolerance
end # End function




function refine_data(r::MyRefinery, cell::Cell, indices)
    cell.data 
end # End function

# =============================================== The querying operators that return intersecting droplets ==================================== # 


function intersects(cell, mydroplet)
    #     println("cell is this: $cell") ; println("myrect is this: $myrect")
        case1x = mydroplet.x - mydroplet.searchradius + 2 * mydroplet.searchradius <= cell.boundary.origin[1]
        case2x = mydroplet.x - mydroplet.searchradius >= cell.boundary.origin[1] + cell.boundary.widths[1] 
        casex = case1x || case2x
    #     println("casex is this: $casex")
        # printstyled("The value and type of the droplet.x is:   " ; color= :green)
        # display(mydroplet.x)
        # display(typeof(mydroplet.x))
        # printstyled("The value and type of the cell.boundary.origin is:   " ; color= :green)
        # display(cell.boundary.origin[1])
        # display(typeof(cell.boundary.origin[1]))

        case1y = mydroplet.y - mydroplet.searchradius + 2 * mydroplet.searchradius <= cell.boundary.origin[2]
        case2y = mydroplet.y - mydroplet.searchradius >= cell.boundary.origin[2] + cell.boundary.widths[2] 
        casey = case1y || case2y
        
        
        case1z = mydroplet.z - mydroplet.searchradius + 2 * mydroplet.searchradius <= cell.boundary.origin[3]
        case2z = mydroplet.z - mydroplet.searchradius >= cell.boundary.origin[3] + cell.boundary.widths[3] 
        casez = case1z || case2z
    #     println("casey is this: $casey")
        
        !(casex || casey || casez)
    
    end # End function



    function contains(point, mydroplet)
        xdist = mydroplet.x - point.x
        ydist = mydroplet.y - point.y
        zdist = mydroplet.z - point.z
        sqrt(xdist^2 + ydist^2 + zdist^2) <= mydroplet.radius+point.radius && (sqrt(xdist^2 + ydist^2 + zdist^2) != 0)
    end  # End function 



    function query(cell, mydroplet, found) # IMPORTANT: Instead of returning a list to get my acceleration values from, I bet I can make query just return the acceleration values itself by tweaking calculations. If I can integrate those two together I save BIG memory allocations
        #     println("found is given to us as $found")
            if intersects(cell, mydroplet) == false
                return found
            else 
                pointvector = cell.data
                if isleaf(cell)
                    for i in 1 : length(pointvector) 
                        contains(pointvector[i], mydroplet) ? push!(found, pointvector[i]) : nothing
                    end # End the for loop
                end # End the conditional


                # The recursive part of the equation ------------- # 
                
                if !isleaf(cell) 
                    query(cell[1,1,1], mydroplet, found)
                    query(cell[1,1,2], mydroplet, found)
                    query(cell[1,2,1], mydroplet, found)
                    query(cell[1,2,2], mydroplet, found)
                    
                    query(cell[2,2,2], mydroplet, found)
                    query(cell[2,2,1], mydroplet, found)
                    query(cell[2,1,2], mydroplet, found)
                    query(cell[2,1,1], mydroplet, found)
                end # End conditional
                
            end # End conditional    
            return found
        end # End function 







        function mygetacceleration(fardrop::Droplet, neardrop::Droplet, accelerationvec::Matrix{Float64})       
            xdist = fardrop.x - neardrop.X
            ydist = fardrop.y - neardrop.y 
            zdist = fardrop.z - neardrop.Z
            directionx = xdist/sqrt(xdist^2 + ydist^2 + zdist^2)
            directiony = ydist/sqrt(xdist^2 + ydist^2 + zdist^2)
            directionz = zdist/sqrt(xdist^2 + ydist^2 + zdist^2)
            mybilayerdist = L * (neardrop.radius + fardrop.radius)
            xaccel = directionx * -k * (abs(mybilayerdist * directionx) - abs(xdist))
            yaccel = directiony * -k * (abs(mybilayerdist * directiony) - abs(ydist)) 
            zaccel = directionz * -k * (abs(mybilayerdist * directionz) - abs(zdist))   
            accelerationvec = [xaccel; yaccel; zaccel]















            # directionx = xdist/sqrt(xdist^2 + ydist^2 + zdist^2)
            

            # # printstyled("The direction X is going to be given as such:" ; color = :green)
            # # display(directionx)
            # directiony = ydist/sqrt(xdist^2 + ydist^2 + zdist^2)
            # directionz = zdist/sqrt(xdist^2 + ydist^2 + zdist^2)
            # # printstyled("The direction Z is going to be given as such:" ; color = :green)
            # # display(directionz)
            # mybilayerdist = L * (droplets[i].radius + founddrops[j].radius)

            # #myaccelx# 
            # # Theoretical ugly code implementation of a full line: (founddrops[j].x - droplets[i].x)/sqrt((founddrops[j].x - droplets[i].x)^2 + (founddrops[j].y - droplets[i].y)^2 + (founddrops[j].z - droplets[i].z)^2) * -k * (abs((L * (droplets[i].radius + founddrops[j].radius) ) * (founddrops[j].x - droplets[i].x)/sqrt((founddrops[j].x - droplets[i].x)^2 + (founddrops[j].y - droplets[i].y)^2 + (founddrops[j].z - droplets[i].z)^2))  -abs(founddrops[j].x - droplets[i].x) )
            # myaccely = directiony * -k * (abs(mybilayerdist * directiony) - abs(ydist)) 
            # myaccelz = directionz * -k * (abs(mybilayerdist * directionz) - abs(zdist)) 
            # myaccelx = directionx * -k * (abs(mybilayerdist * directionx) - abs(xdist))
            # myaccely = directiony * -k * (abs(mybilayerdist * directiony) - abs(ydist)) 
            # myaccelz = directionz * -k * (abs(mybilayerdist * directionz) - abs(zdist)) 
            # droplets[i].xacceleration += myaccelx
            # droplets[i].yacceleration += myaccely
            # droplets[i].zacceleration += myaccelz







        end # End function







































        
end # Module end 


























