module MyDropletVectorJComponents
include("myparticlepoints.jl")
# include("myquadtreequery.jl")
using Interpolations
using LinearAlgebra
using RegionTrees
using StaticArrays: SVector
import RegionTrees: AbstractRefinery, needs_refinement, refine_data, adaptivesampling!
using .ParticlePoints
using DifferentialEquations



# import RegionTrees: AbstractRefinery, needs_refinement, refine_data


# The exports of this particular module 
# export mydropfun!, myOvitoPrint!
export dropvecj!


# using StaticArrays: SVector

function dropvecj!(ρ, k, L , γ, D, timeend, ICmatrix, boolean_osm)
    ICvector = reshape(ICmatrix, :, 1);
    # display(ICvector)
    # Define some of the model parameters 
    outerdropletcount = size(ICmatrix,1)
    # display(outerdropletcount)
    tspan = (0.0, timeend);
    # The 3 value is the number of dimensions, x, y, z 
    dradii = zeros(outerdropletcount,1)
    dosmolarities = zeros(outerdropletcount,1)
    dvelocities = zeros(3*outerdropletcount,1)
    # maxradius::Vector{Float64} = zeros(1)
    accelerationvec = zeros(3,1)
    xyzaccel::Matrix{Float64} = zeros(3,1)
    flux = zeros(1,1)
    # isdefined(dropletcount)
    
      # Create the vector of droplets that will be fed into the main ODE 
      droplets::Vector{Droplet} = Droplet[] # Change this to a vector that needs to be imported in 
      for i in 1:outerdropletcount
          droplet::Droplet = Droplet(0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)  #< ----- > THE FORMATTING FOR THE DROPLET
          push!(droplets, droplet)
      end # End the for loop
    # typeof(droplets)
    # display(droplets)
    # Create the current residual that makes the osmolarity activate after the droplets settle
    #   normvelocity = zeros(1,1)



    p = [outerdropletcount, k, L, γ, D, boolean_osm, dradii, dosmolarities, dvelocities, ρ, accelerationvec, droplets, xyzaccel, flux] 
    prob = ODEProblem(mydropfun!, ICvector, tspan, p)
    sol = solve(prob, Tsit5(),saveat = 0.1, dt = 0.1)
    # graphics = mySplineConstructor(sol, dropletcount, dimensions, tstepextract)
    W = myOvitoPrint!(sol, outerdropletcount, 3)
    outfile = "E:/Research Scripts and Functions/OVITO_Files/visual_model.ovito"
    f = open(outfile, "w")
    for i in eachindex(W)   
        println(f, W[i])
    end
    close(f)
    end # Function end

function mydropfun!(du, u, p, t)
    # Put in the parameters that are needed, very important line of code: #=========================================================================
    outerdropletcount, k, L, γ, D, boolean_osm, dradii, dosmolarities, dvelocities, ρ, accelerationvec, droplets, xyzaccel, flux = p;
    # ================================================================================================ # 
    # println("NEW FRAME # -------------------------------------------------------------------------------------------------------------- #")
    # display(typeof(u))
    maxradius::Float64 = 0. 
    dropletcount = outerdropletcount 
    # for i in 1:dropletcount
    #     u[i] > maxradius ? maxradius = u[i]::Float64 : nothing 
    # end # End for loop 
    
    maxradius = maximum(view(u,1:dropletcount))

    # display(typeof(maxradius))
    # a = typeof(du)
    # display(a)

    normvelocity::Float64 = 0.
    # vecvelocity = u[5*dropletcount+1:8*dropletcount].^2
    ℷ = abs.(view(u, 5*dropletcount+1:8*dropletcount))
    ℸ = view(u, 5*dropletcount+1:8*dropletcount)
    Ζ = view(u, 5*dropletcount+1:8*dropletcount)
    broadcast(abs, ℸ)
    ℶ = view(u, 5*dropletcount+1:8*dropletcount).^2
    ℵ::Float64 = sum(ℶ)::Float64
    ℵ  = sqrt(ℵ)
    # for i in 1:ℵ
    # for i in 1:ℵ
    #     display(typeof(view(u, 5*dropletcount+i,:)^2))
    #     normvelocity::Float64 += view(u, 5*dropletcount+i).^2
    # end # End for loop
    # printstyled("The type of the normvelocity prior to the square root operator is:"; color =:green)
    # display(typeof(normvelocity))
    # printstyled("The type of the normvelocity after  the square root operator is:"; color =:green)
    # display(typeof(normvelocity))




    # display(mynormvelocity(u[5*dropletcount+1:8*dropletcount]))

    if (t > 0.1 && normvelocity/dropletcount <= 1e-4)  
        boolean_osm[1] = 1.0
    end # End conditional

    for i in 1:dropletcount
        # display(typeof(droplets[i]))
        droplets[i].radius = getindex(u,i,1)::Float64
        droplets[i].osmolarity = getindex(u, i+dropletcount)::Float64
        droplets[i].x = getindex(u, Int(i+2*dropletcount))::Float64
        droplets[i].y = getindex(u, Int(i+3*dropletcount))::Float64
        droplets[i].z = getindex(u, Int(i+4*dropletcount))::Float64
        droplets[i].xvelocity =  getindex(u, Int(i+5*dropletcount))::Float64
        droplets[i].yvelocity =  getindex(u, Int(i+6*dropletcount))::Float64
        droplets[i].zvelocity =  getindex(u, Int(i+7*dropletcount))::Float64
        droplets[i].searchradius =  10.
        # display(typeof(droplets[i]))
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
    r = MyRefinery(0.0)
    asaf = TreeCellData(droplets, dropletcount)
    root::Cell = Cell(SVector(xorigin, yorigin, zorigin), SVector(xwidth, ywidth, zwidth), asaf)
    myquadtree = adaptivesampling!(root, r)

    for i in 1:dropletcount
        
        fill!(xyzaccel, 0.)
        flux[1] = 0.
        jvec = 0.
        query(myquadtree, droplets[i], xyzaccel, flux, L, k, boolean_osm, D, jvec, accelerationvec)


        myvolume = 4/3 * pi * droplets[i].radius^3
        dvelocities[i] = (xyzaccel[1] - γ * droplets[i].xvelocity) / (ρ * myvolume) 
        dvelocities[i+dropletcount] = (xyzaccel[2] - γ * droplets[i].yvelocity) / (ρ * myvolume) 
        dvelocities[i+2*dropletcount] = (xyzaccel[3] - γ * droplets[i].zvelocity) / (ρ * myvolume) 

        mymoles = droplets[i].osmolarity * myvolume
        myvolume = myvolume + flux[1]  
        dradii[i] = (3/(4*pi) * myvolume)^(1/3) - droplets[i].radius
        dosmolarities[i] = mymoles/myvolume - droplets[i].osmolarity
        
    end # End the for loop

    du[1:dropletcount] = dradii
    du[dropletcount+1 : 2*dropletcount] = dosmolarities
    du[2*dropletcount+1 : 5*dropletcount] = view(u, 5*dropletcount+1:8*dropletcount)
    du[5*dropletcount+1 : 8*dropletcount] = dvelocities

    end # End the main function



# =============================================== All the functions that really should go inside their own module ========================# 
#=========================================================================================================================================#
# =============================================== The refinement structure that creates our quad tree ====================================# 
struct MyRefinery <: AbstractRefinery
    tolerance::Int64
end # End struct


function needs_refinement(r::MyRefinery, cell)
    
    # display(cell.data.endindex)
    # test
    # display(typeof(pointvector))
#     display(pointvector)
    myindex::Int64 = 0 
    for i in 1:cell.data.endindex
        # Adress an edge case where the point sits directly in top of the line of division
        edgecorrectionx = 0.
        edgecorrectiony = 0.
        edgecorrectionz = 0. 
        
        edgecasex1::Bool = cell.data.list[i].x == cell.boundary.origin[1] 
        edgecasex2::Bool = cell.data.list[i].x == cell.boundary.origin[1] + cell.boundary.widths[1] 
        edgecasey1::Bool = cell.data.list[i].y == cell.boundary.origin[2] 
        edgecasey2::Bool = cell.data.list[i].y == cell.boundary.origin[2] + cell.boundary.widths[2] 
        edgecasez1::Bool = cell.data.list[i].z == cell.boundary.origin[3] 
        edgecasez2::Bool = cell.data.list[i].z == cell.boundary.origin[3] + cell.boundary.widths[3] 
        
        edgecasex1 ? edgecorrectionx::Float64 += (cell.boundary.widths[1]+10)^(-1/cell.boundary.widths[1]) : nothing 
        edgecasex2 ? edgecorrectionx::Float64 += (cell.boundary.widths[1]+10)^(-1/cell.boundary.widths[1]) : nothing
        edgecasey1 ? edgecorrectiony::Float64 += (cell.boundary.widths[2]+10)^(-1/cell.boundary.widths[2]) : nothing
        edgecasey2 ? edgecorrectiony::Float64 += (cell.boundary.widths[2]+10)^(-1/cell.boundary.widths[2]) : nothing
        edgecasez1 ? edgecorrectionz::Float64 += (cell.boundary.widths[3]+10)^(-1/cell.boundary.widths[3]) : nothing
        edgecasez2 ? edgecorrectionz::Float64 += (cell.boundary.widths[3]+10)^(-1/cell.boundary.widths[3]) : nothing
      
        case1 = cell.boundary.origin[1] < cell.data.list[i].x + edgecorrectionx < cell.boundary.origin[1] + cell.boundary.widths[1]  
        case2 = cell.boundary.origin[2] < cell.data.list[i].y + edgecorrectiony < cell.boundary.origin[2] + cell.boundary.widths[2] 
        case3 = cell.boundary.origin[3] < cell.data.list[i].z + edgecorrectionz < cell.boundary.origin[3] + cell.boundary.widths[3] 
        if case1 && case2 && case3 == true 
            myindex += Int(1)
            cell.data.list[myindex] = cell.data.list[i] 
        end # end conditional
    end # End the for loop
    
    cell.data.endindex = myindex
    myindex == r.tolerance
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



    function query(cell, mydroplet, xyzaccel, flux, L, k, boolean_osm, D, jvec, accelerationvec) 
            if intersects(cell, mydroplet) == false
                return nothing
            else 
                if isleaf(cell)
                    for i in 1 : length(cell.data.list) 
                        contains(cell.data.list[i], mydroplet) ? xyzaccel .+= mygetaccelerationandflux(cell.data.list[i], mydroplet, L, k, boolean_osm, D, jvec, accelerationvec)[1]  : nothing
                        contains(cell.data.list[i], mydroplet) ? flux .+= mygetaccelerationandflux(cell.data.list[i], mydroplet, L, k, boolean_osm, D, jvec, accelerationvec)[2]  : nothing
                    end # End the for loop



                end # End the conditional


                # The recursive part of the equation ------------- # 
                
                if !isleaf(cell) 
                    query(cell[1,1,1], mydroplet, xyzaccel, flux, L, k, boolean_osm, D, jvec, accelerationvec)
                    query(cell[1,1,2], mydroplet, xyzaccel, flux, L, k, boolean_osm, D, jvec, accelerationvec)
                    query(cell[1,2,1], mydroplet, xyzaccel, flux, L, k, boolean_osm, D, jvec, accelerationvec)
                    query(cell[1,2,2], mydroplet, xyzaccel, flux, L, k, boolean_osm, D, jvec, accelerationvec)
                    
                    query(cell[2,2,2], mydroplet, xyzaccel, flux, L, k, boolean_osm, D, jvec, accelerationvec)
                    query(cell[2,2,1], mydroplet, xyzaccel, flux, L, k, boolean_osm, D, jvec, accelerationvec)
                    query(cell[2,1,2], mydroplet, xyzaccel, flux, L, k, boolean_osm, D, jvec, accelerationvec)
                    query(cell[2,1,1], mydroplet, xyzaccel, flux, L, k, boolean_osm, D, jvec, accelerationvec)
                end # End conditional
                
            end # End conditional    
            # Perform the final calculations 




            return xyzaccel
            return flux
        end # End function 

        function mygetaccelerationandflux(fardrop::Droplet, neardrop::Droplet, L, k, boolean_osm, D, jvec, accelerationvec)#::Matrix{Float64}       
            xdist = fardrop.x - neardrop.x
            ydist = fardrop.y - neardrop.y 
            zdist = fardrop.z - neardrop.z
            directionx = xdist/sqrt(xdist^2 + ydist^2 + zdist^2)
            directiony = ydist/sqrt(xdist^2 + ydist^2 + zdist^2)
            directionz = zdist/sqrt(xdist^2 + ydist^2 + zdist^2)
            mybilayerdist = L * (neardrop.radius + fardrop.radius)
            # xaccel = directionx * -k * (abs(mybilayerdist * directionx) - abs(xdist))
            # yaccel = directiony * -k * (abs(mybilayerdist * directiony) - abs(ydist)) 
            # zaccel = directionz * -k * (abs(mybilayerdist * directionz) - abs(zdist))   

            
            # @.accelerationvec .= [xaccel; yaccel; zaccel]
            


            accelerationvec[1] = directionx * -k * (abs(mybilayerdist * directionx) - abs(xdist))
            accelerationvec[2] = directiony * -k * (abs(mybilayerdist * directiony) - abs(ydist)) 
            accelerationvec[3] = directionz * -k * (abs(mybilayerdist * directionz) - abs(zdist))   

            # display(accelerationvec)
            # printstyled("accelerationvec[1] is given to us as $(accelerationvec)")



            # =================== The osmolarity part of the code ================================= # 
            intersectcase1 =  neardrop.radius + fardrop.radius - sqrt(xdist^2 + ydist^2 + zdist^2) >= 0
            intersectcase2 = sqrt(xdist^2 + ydist^2 + zdist^2) + neardrop.radius  - fardrop.radius >= 0
            intersectcase3 = sqrt(xdist^2 + ydist^2 + zdist^2) + fardrop.radius  - neardrop.radius >= 0

            # jvec = 0 
            if boolean_osm[1] >= 1.0 
                if (intersectcase1 && intersectcase2 && intersectcase3)
                    component1 = (sqrt(xdist^2 + ydist^2 + zdist^2)^2 - fardrop.radius^2 + neardrop.radius^2)^2
                    component2 = 4 * sqrt(xdist^2 + ydist^2 + zdist^2)^2 * neardrop.radius^2
                    component3 = component2 - component1
                    component4 = 1/(2*sqrt(xdist^2 + ydist^2 + zdist^2)) * sqrt(component3)
                    mycircular_area = pi * component4^2
                    osmdiff = neardrop.osmolarity - fardrop.osmolarity
                    jvec = D * mycircular_area * osmdiff
                end # End the conditional
            end # End the conditional
            # printstyled("jvec is "; color =:green)
            # display(jvec)
            # display(typeof(jvec))

            return accelerationvec, jvec
        end # End function




        # function mymaxradius()
        #     for i in 1:length(u[5*dropletcount+1:8*dropletcount])
        #         normvelocity[1] += u[5*dropletcount+i]^2
        #     end # End for loop
        #     normvelocity = sqrt(normvelocity)















#================================================================================ OVITO PRINT FUNCTION ================================================================================================#

function myOvitoPrint!(sol, dropletcount, dimensions)
    # x = Array{Float64}(undef, 1,1)
    # y = Array{Float64}(undef, 1,1)
    # z = Array{Float64}(undef, 1,1)
    # radius = Array{Float64}(undef, 1,1)
    # osmolarity = Array{Float64}(undef, 1,1)


    # This function will recast the solution of our ordinary differential equation into a printed file that can be visualized in ovito.
    
    V = convert(Array,sol) # Convert the ODE solution into an array that can be worked on 
    V = dropdims(V, dims = 2) # Dimension 2 is a singleton and should be eliminated
    # display(V)
    totalsteps::Int64 = size(V,2)
    # Create the empty string vector that will be used to create the printout 
    myprintout = Array{String}(undef, (5+dropletcount)*totalsteps, 1)
    # Initialize the timestep 
    timestep::Float64 = 0.
    stepindex::Int64 = 1

    for i in 1:  5+dropletcount:  (5+dropletcount)*totalsteps
        myprintout[i] = "ITEM: TIMESTEP"
        myprintout[i+1] = "$timestep"
        myprintout[i+2] = "ITEM: NUMBER OF ATOMS"
        myprintout[i+3] = "$dropletcount"
        myprintout[i+4] = "ITEM: ATOMS x y z radius osmolarity"
       for j in 1:dropletcount
            # myprintout[i+4+j] = "$x $y $z $radius $osmolarity"
            x = V[2*dropletcount+j, stepindex]
            y = V[3*dropletcount+j, stepindex]
            z = V[4*dropletcount+j, stepindex]
            radius = V[j, stepindex]
            osmolarity = V[dropletcount + j, stepindex]
            myprintout[i+4+j] = "$x $y $z $radius $(osmolarity)"
       end # End sub-loop for droplets
       stepindex += 1 
       timestep += 0.1 
    end # End loop
    return myprintout
    end # Function end 


end # Module end 
























