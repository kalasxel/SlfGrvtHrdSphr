function Force( r1::CartesianVector{T}, r2::CartesianVector{T} ) :: CartesianVector{T} where {T<:Real}
    return -(G/AbsVec( r1 - r2 )^3) * (r1 - r2)    
end


function SetAccelerations!( prtcls::Vector{Particle{T}} ) where {T<:Real}

    for pi in prtcls
        acceleration = CartesianVector{Float64}( [0e0, 0e0, 0e0] );
        for pj in prtcls
            if pi != pj
                #acceleration = acceleration + (G/AbsVec( pj.r - pi.r )^3) * (pj.r - pi.r)
                acceleration = acceleration + Force(pi.r, pj.r)
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

    for pi in prtcls
        pi.v = pi.v + pi.a *dt / convert(T,2)       
        pi.r = pi.r + pi.v *dt        
    end

    SetAccelerations!(prtcls)
    for pi in prtcls
        pi.v = pi.v + pi.a *dt / convert(T,2)  
    end

end


const c1 = c4 = 0.67560359597982881702384390448573041
const c2 = c3 = -0.17560359597982881702384390448573041
const d3 = d1 = 1.3512071919596576340476878089714608
const d2 = -1.7024143839193152680953756179429217

function YoshidaStep!( prtcls::Vector{Particle{T}} ) where {T<:Real}

    SetAccelerations!(prtcls)
    for pi in prtcls
        pi.r = pi.r + pi.v *dt *c1  # x1
        pi.v = pi.v + pi.a *dt *d1  # v1
        pi.r = pi.r + pi.v *dt *c2  # x2       
    end

    SetAccelerations!(prtcls)
    for pi in prtcls
        pi.v = pi.v + pi.a *dt *d2  # v2
        pi.r = pi.r + pi.v *dt *c3  # x3         
    end

    SetAccelerations!(prtcls)
    for pi in prtcls
        pi.v = pi.v + pi.a *dt *d3  # v3 = v4
        pi.r = pi.r + pi.v *dt *c4  # x4        
    end

end