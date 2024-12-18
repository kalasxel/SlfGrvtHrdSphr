mutable struct CartesianVector{T<:Real} # to define components of vector
    x :: T
    y :: T
    z :: T
    function CartesianVector{T}( x=convert(T,0)::T, y=convert(T,0)::T, z=convert(T,0)::T ) where {T<:Real}
        new(x,y,z)
    end
    function CartesianVector{T}( r :: Vector{T} ) where {T<:Real}
        new(r[1],r[2],r[3])
    end
    
end

#   Some basic functions for vectors
#   All arguments in these functions must have the same type 

function Base.:+( a::CartesianVector{T}, b::CartesianVector{T} ) :: CartesianVector{T} where {T<:Real} # sum of two vectors
    return CartesianVector{T}( a.x+b.x, a.y+b.y, a.z+b.z )
end

function Base.:-( a::CartesianVector{T}, b::CartesianVector{T} ) :: CartesianVector{T} where {T<:Real} # difference of two vectors
    return CartesianVector{T}( a.x-b.x, a.y-b.y, a.z-b.z )
end

function Base.:*( s::T, b::CartesianVector{T} ) :: CartesianVector{T} where {T<:Real} # left product with scalar
    return CartesianVector{T}( s*b.x, s*b.y, s*b.z )
end

function Base.:*( b::CartesianVector{T}, s::T ) :: CartesianVector{T} where {T<:Real} # right product with scalar
    return CartesianVector{T}( s*b.x, s*b.y, s*b.z )
end 

function Base.:/( b::CartesianVector{T}, s::T ) :: CartesianVector{T} where {T<:Real} # divide vector by scalar
    return CartesianVector{T}( b.x/s, b.y/s, b.z/s )
end 

function ScalProd( a::CartesianVector{T}, b::CartesianVector{T} ) :: T where {T<:Real}  # scalar product of two vectors 
    return a.x*b.x + a.y*b.y + a.z*b.z    
end

function AbsVec( a::CartesianVector{T} ) :: Real where {T<:Real} # absolute value of a vector
    return sqrt( ScalProd(a,a) )
end

function VecProd( a::CartesianVector{T}, b::CartesianVector{T} ) :: CartesianVector{T} where {T<:Real}  # vector product of two vectors
    return CartesianVector{T}( a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x )
end


abstract type Particles end

mutable struct Particle{T<:Real} <: Particles # everything that particle has
    m :: T
    radius :: T
    r :: CartesianVector{T}
    v :: CartesianVector{T}
    a :: CartesianVector{T}
    function Particle{T}( m=convert(T,0)::T, radius=convert(T,0)::T, 
                          r=CartesianVector{T}()::CartesianVector{T},
                          v=CartesianVector{T}()::CartesianVector{T},
                          a=CartesianVector{T}()::CartesianVector{T} ) where {T<:Real}
        new( m, radius, r, v, a )
    end
end

#   Some functions for a particle
#   All arguments in these functions must have the same type 

function Momentum( prtcl::Particle{T} ) :: CartesianVector{T} where {T<:Real} # momentum of particle
    return prtcl.m * prtcl.v
end

function AngularMomentum( r::CartesianVector{T}, prtcl::Particle{T} ) :: CartesianVector{T} where {T<:Real} # angular momentum of particle relative to r
    return VecProd( r, Momentum(prtcl) )
end

function KineticEnergy( prtcl::Particle{T} ) :: T where {T<:Real} # kinetic energy of particle
    return ScalProd( Momentum(prtcl), prtcl.v ) / convert(T,2)
end

function TotKinEnergy( prtcls::Vector{Particle{T}} ) :: T where {T<:Real}
    ϵKin = convert(T,0)
    for p in prtcls
        ϵKin += KineticEnergy(p)
    end
    return ϵKin
end

function Potential( r::CartesianVector{T}, prtcl::Particle{T} ) :: T where {T<:Real}  # potential of a particle "prtcl" in a point "r"
    return -G*prtcl.m / AbsVec( r - prtcl.r )    
end

function TotPotEnergy( prtcls::Vector{Particle{T}} ) :: T where {T<:Real} # calculates full potential energy
    ϵPot = convert(T,0)
    for pi in prtcls
        for pj in prtcls
            if pi != pj
                ϵPot += pi.m*pj.m / AbsVec( pi.r - pj.r )
            end
        end
    end
    return -G*ϵPot/convert(T,2)
end
