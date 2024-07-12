function SetAccelerations!( prtcls::Vector{Particle{T}} ) where {T<:Real}

    for pi in prtcls
        acceleration = CartesianVector{Float64}( [0e0, 0e0, 0e0] );
        for pj in prtcls
            if pi != pj
                acceleration = acceleration + (G/AbsVec( pj.r - pi.r )^3) * (pj.r - pi.r)
            end
        end
        pi.a = acceleration
    end
    
end


function EulerStep!( prtcls::Vector{Particle{T}} ) where {T<:Real}

    SetAccelerations!(prtcls)

    for pi in prtcls
        pi.v = pi.v + pi.a * dt
        pi.r = pi.r + pi.v * dt
    end

end


function LeapfrogStep!( prtcls::Vector{Particle{T}} ) where {T<:Real}

    prtclsOld = deepcopy(prtcls) :: Vector{Particle{T}}

    for pi in prtcls
        pi.r = pi.r + pi.v*dt + convert(T,2)*pi.a*dt*dt
    end

    SetAccelerations!(prtcls)

    for k in 1:nb
        prtcls[k].v = prtcls[k].v + ( prtcls[k].a + prtclsOld[k].a )*dt / convert(T,2)
    end

end