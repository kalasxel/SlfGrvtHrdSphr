using DelimitedFiles
import Combinatorics: powerset

additionalConstrain = false

const workingDir = dirname(@__FILE__)
cd(workingDir)


# reading previously obtained data

prms = readdlm( workingDir* "/OUTPUT/PARAMS.dat", Float64 )
const M=prms[1,1]; const R=prms[1,2];
const xi=prms[2,1]; const xf=prms[2,2];
const yi=prms[3,1]; const yf=prms[3,2];
const zi=prms[4,1]; const zf=prms[4,2];
const nb=convert(Int,prms[5,1]); const to=convert(Int,prms[5,2])
const G=prms[6,1]; const σ=prms[6,2];

stc = 500 :: Int # considered time moment
if stc >= to  
    error("Considered step ", stc, " must be less than ", to, ".")
end

ptmp = zeros(4,nb) # [x1, x2, x3, E]

@inbounds for i in 1:nb
    ptmp[:,i] = readdlm( workingDir* "/OUTPUT/particle"* string(i)* ".dat", Float64 )[stc,2:5]
end


# forming array of particles params --- p[1:nb]

struct Prtcl{T<:Real}
    x :: T
    y :: T
    z :: T
    e :: T
    function Prtcl{T}( x=convert(T,0)::T, y=convert(T,0)::T, z=convert(T,0)::T,
            e=convert(T,0)::T ) where {T<:Real}
    new( x, y, z, e )
    end
end

p = [ Prtcl{Float64}(ptmp[1,i],ptmp[2,i],ptmp[3,i],ptmp[4,i]) for i in 1:nb ]

#println("End reading")
#comb = deleteat!(collect(powerset([i for i=1:nb])),1:nb+1) 
#println("The size is ", size(comb))


function AreTheyBound( p::Vector{Prtcl{T}}, composition::Vector{Int64} ) :: Bool where {T<:Real}

  # if each particle has negative energy in the field of all the others in the subset
    if additionalConstrain

        for j in composition

            ϵ = 0e0            
            for i in composition
                if i != j
                    ϵ += - G*M*M / sqrt( (p[j].x-p[i].x)^2 + (p[j].y-p[i].y)^2 + (p[j].z-p[i].z)^2 )
                end
            end 

            if ϵ + p[j].e >= 0e0
                return false
            end

        end

    end


    # the total energy (kinetic + potential) of the subset
    ϵkin = ϵpot = 0e0
    for j in composition

        ϵkin += p[j].e 

        for i in composition
            if i != j
                ϵpot += -(5e-1)*G*M*M/sqrt( (p[j].x-p[i].x)^2 + (p[j].y-p[i].y)^2 + (p[j].z-p[i].z)^2 )
            end
        end

    end

    if ϵkin + ϵpot <= 0
        return true
    else 
        return false
    end    

end


function BoundSubsets( p::Vector{Prtcl{T}} ) :: Vector{T} where {T<:Real}

    # generating index array with all possible subsets exluding empty and one-element sets
    comb = deleteat!(collect(powerset([i for i=1:nb])),1:nb+1) 
    len = length(comb)

    nBound = [0 for i in 1:nb-1] # number of bound coupls [1], threes[2], ..., nb-s [nb-1]
    nBound0 = deepcopy(nBound)

    for i = 1:len
        if AreTheyBound( p, comb[i] )
            nBound[length(comb[i])-1] += 1
        end
        nBound0[length(comb[i])-1] += 1
    end

    return nBound ./ nBound0

end
@time bc = BoundSubsets(p)

open( workingDir* "/OUTPUT/NbBoundClusters"* string(stc) *"_sigma" * string(σ) * ".dat", "w" ) do io writedlm(io, [(i+1, bc[i]) for i in 1:nb-1]) end


