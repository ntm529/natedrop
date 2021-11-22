function myLoopMatrix(radii ,dropletcount, osmolaritybase, osmolaritytop)

radiivector = radii * ones(dropletcount,1)

dropletnum = Int(dropletcount/2)

ix = range(0, step = 2* radii, stop =  2 * radii * dropletnum - radii)
ix = collect(ix)
iy = zeros(1,dropletnum)
iz = zeros(1,dropletnum)

jx = range(radii, step = 2*radii, stop = 2*radii*dropletnum)
jx = collect(jx)
jy = repeat([sqrt(3*radii^2)], 1, dropletnum)
jz = zeros(1,dropletnum)

velocities = zeros(dropletcount, 3)


loopmatrix = [[[ix  iy'  iz'] ; [jx  jy' jz']]  velocities]


osmolaritybase = osmolaritybase * ones(dropletnum,1)
osmolaritytop = osmolaritytop * ones(dropletnum,1)
osmolarities = [osmolaritybase ; osmolaritytop]

loopmatrix = [radiivector osmolarities loopmatrix ]

end


