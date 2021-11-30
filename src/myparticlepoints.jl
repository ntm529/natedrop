module ParticlePoints

mutable struct Droplet
    radius::Float64
    osmolarity::Float64
    x::Float64
    y::Float64
    z::Float64
    xvelocity::Float64
    yvelocity::Float64
    zvelocity::Float64
    xacceleration::Float64
    yacceleration::Float64
    zacceleration::Float64
    searchradius::Float64
end # End the struct


mutable struct TreeCellData
    list::Vector{Droplet}
    endindex::Int64
end # End the struct





export Droplet, TreeCellData
end # End the module

