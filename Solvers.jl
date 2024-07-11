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


function VerletStep!( prtcls::Vector{Particle{T}}, prtclsOld::Vector{Particle{T}} ) where {T<:Real}

    SetAccelerations!(prtcls)

    tmp = deepcopy(prtcls)

    for i in 1:nb
        prtcls[i].r = convert(T,2)*prtcls[i].r - prtclsOld[i].r + prtcls[i].a * dt*dt
        prtcls[i].v = ( prtcls[i].r - prtclsOld[i].r )/(convert(T,2)*dt)
    end

    prtclsOld = deepcopy(tmp)

end