function mySheetMatrix(r, n)

    ix = collect(0:2*r:2*(n-1)*r)'
    iy = repeat(ix,n,1)
    iy = reshape(iy,:,1)

    ix = repeat(ix',n,1)
    iz = zeros(length(ix),1)
    i = [ix iy iz]
    display(i)
    velocities = zeros(length(ix),3)


    
    jx = collect(r:2*r:r+2*(n-2)*r)'

    jy = repeat(jx,n-1,1)
    jy = reshape(jy,:,1)
    print("jy is ",jy )
    display(jy)


    jx = repeat(jx',n-1,1)
    display(jx)

    testrad = r * ones(length(ix)+length(jx),1)
    display(testrad)

    jz = (2* testrad[length(ix)+1]).^2
    display(jz)

    jz = jz - jx[1].^2 - jy[1].^2
    display(jz)

    jz = sqrt(jz)
    display(jz)
    jz = jz * ones(length(jx),1)
    display(jz)
    j = [jx jy jz]


    testosm = ones(length(ix)+length(jx),1)
    testosm[length(ix)+1:end] .= 0.1 

    velocities = zeros(length(ix)+length(jx),3)


    firstmat = [testrad testosm [i;j] velocities]

    return firstmat


end



function myStableSheetMatrix(r, n, L)

    bottomstart = 0 
    topstart = bottomstart + r*L

    bottomend = 2*(n)*r*L
    ix = zeros(1,n)
  

    for i in 1:length(ix)
        ix[i] = bottomstart
        bottomstart += 2*r*L
        # display(ix)

    end # End the loop

    
    jx = zeros(1,n-1)
    for i in 1:length(jx)
        jx[i] = topstart
        topstart += 2*r*L
        # display(jx)
    end # End the loop
    # display(jx)




    printstyled("ix is given to us as"; color = :green)
    display(ix)
    printstyled("jx is given to us as"; color = :green)
    display(jx)
    bottomstart = 0 

    # ix = collect(bottomstart:          2*r*L:         bottomend)'
    # printstyled("ix is given to us as"; color = :green)
    # display(ix)
    iy = repeat(ix,n,1)
    iy = reshape(iy,:,1)

    ix = repeat(ix',n,1)
    iz = zeros(length(ix),1)
    i = [ix iy iz]
    # display(i)
    velocities = zeros(length(ix),3)


    
    # jx = collect(bottomstart + r*L:                2*r*L:              bottomend - r*L)'
    # printstyled("jx is given to us as"; color = :green)
    # display(jx)

    jy = repeat(jx,n-1,1)
    jy = reshape(jy,:,1)
    # print("jy is ",jy )
    # display(jy)


    jx = repeat(jx',n-1,1)
    display(jx)

    testrad = r * ones(length(ix)+length(jx),1)
    # printstyled("tesrad is given to us as"; color = :green)
    # display(testrad)

    jz = (2* testrad[length(ix)+1]*L).^2
    # printstyled("jz is given to us as"; color = :green)
    # display(jz)

    jz = jz - jx[1].^2 - jy[1].^2
    # printstyled("jz is given to us as"; color = :green)
    # display(jz)

    jz = sqrt(jz)
    # printstyled("jz is given to us as"; color = :green)
    # display(jz)
    jz = jz * ones(length(jx),1)
    # printstyled("jz is given to us as"; color = :green)
    # display(jz)
    j = [jx jy jz]


    testosm = ones(length(ix)+length(jx),1)
    testosm[length(ix)+1:end] .= 0.1 

    velocities = zeros(length(ix)+length(jx),3)


    firstmat = [testrad testosm [i;j] velocities]

    return firstmat


end



