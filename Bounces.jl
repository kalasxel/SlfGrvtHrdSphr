function BounceBoundaries!( prtcl::Particle{T} ) where {T<:Real}

    if prtcl.r.x > xf - prtcl.radius
        #delta = prtcl.r.x + prtcl.radius - xf # distence exeed
        delta = convert(T,0)
        prtcl.r.x = xf - prtcl.radius + delta # move it to boundary and compensate distence exeed
        prtcl.v.x = -prtcl.v.x # change velocity directon 
    end
    if prtcl.r.x < xi + prtcl.radius
        #delta = prtcl.r.x - prtcl.radius - xi
        delta = convert(T,0)
        prtcl.r.x = xi + prtcl.radius + delta
        prtcl.v.x = -prtcl.v.x 
    end

    if prtcl.r.y > yf - prtcl.radius
        #delta = prtcl.r.y + prtcl.radius - zf
        delta = convert(T,0)
        prtcl.r.y = yf - prtcl.radius + delta
        prtcl.v.y = -prtcl.v.y
    end
    if prtcl.r.y < yi + prtcl.radius
        #delta = prtcl.r.y - prtcl.radius - yi
        delta = convert(T,0)
        prtcl.r.y = yi + prtcl.radius + delta
        prtcl.v.y = -prtcl.v.y
    end

    if prtcl.r.z > zf - prtcl.radius
        #delta = prtcl.r.z + prtcl.radius - zf
        delta = convert(T,0)
        prtcl.r.z = zf - prtcl.radius
        prtcl.v.z = -prtcl.v.z
    end
    if prtcl.r.z < zi + prtcl.radius
        #delta = prtcl.r.z - prtcl.radius - zi
        delta = convert(T,0)
        prtcl.r.z = zi + prtcl.radius + delta
        prtcl.v.z = -prtcl.v.z
    end

end


function BounceParticles!( pi::Particle{T}, pj::Particle{T} ) where {T<:Real}

    delta = AbsVec( pi.r - pj.r ) - ( pi.radius + pj.radius ) # distence exeed
    if delta <= 0

        mij = 2*pi.m*pj.m / ( pi.m + pj.m ) # reduced mass
        nij = ( pi.r - pj.r ) / AbsVec( pi.r - pj.r ) # normal to both spheres
        vij = pi.v - pj.v # relative velocity

        dvi = -(mij/pi.m) * ScalProd(nij, vij) * nij 
        dvj =  (mij/pj.m) * ScalProd(nij, vij) * nij 

        pi.v = pi.v + dvi
        pj.v = pj.v + dvj

        pi.r = pi.r - nij * delta / convert(T,2) # shift pi only to half of delta in nij direction
        pj.r = pj.r + nij * delta / convert(T,2) # shift it in another direction

    end
end
