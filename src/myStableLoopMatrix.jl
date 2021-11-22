function myStableLoopMatrix(radii , L, dropletcount, osmolaritybase, osmolaritytop)

    radiivector = radii * ones(dropletcount,1)
    
    dropletnum = Int(dropletcount/2)
    
    ix = range(0, step = 2* radii*L, stop =  (2 * radii * dropletnum - radii)*L )
    ix = collect(ix)
    iy = zeros(1,dropletnum)
    iz = zeros(1,dropletnum)
    
    jx = range(radii*L, step = 2*radii*L, stop = 2*radii*dropletnum*L)
    jx = collect(jx)
    jy = repeat([sqrt(3*radii^2)*L], 1, dropletnum)
    jz = zeros(1,dropletnum)
    
    velocities = zeros(dropletcount, 3)
    
    
    loopmatrix = [[[ix  iy'  iz'] ; [jx  jy' jz']]  velocities]
    
    
    osmolaritybase = osmolaritybase * ones(dropletnum,1)
    osmolaritytop = osmolaritytop * ones(dropletnum,1)
    osmolarities = [osmolaritybase ; osmolaritytop]
    
    stableloopmatrix = [radiivector osmolarities loopmatrix ]
    
    end
    
    
    