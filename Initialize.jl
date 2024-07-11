function Intersection( pi::Particle{T}, pj::Particle{T} ) :: Bool where {T<:Real} # do pi and pj intersect? 
    if AbsVec( pi.r - pj.r ) - ( pi.radius + pj.radius ) < convert(T,0)
        return true
    else 
        return false
    end
end


function makeInitialSpatialDistribution!( prtcls::Vector{Particle{T}} ) :: Bool where {T<:Real}

    for pi in prtcls
        pi.r = CartesianVector{Float64}( rand(Uniform(xi+R,xf-R),3) )
        for pj in prtcls        
            if pi != pj
                while Intersection(pi,pj)
                    pi.r = CartesianVector{Float64}( rand(Uniform(xi+R,xf-R),3) )
                end
            end
        end
    end

    global doTheyIntersect = false :: Bool

    for pi in prtcls
        for pj in prtcls
            if pi != pj
                tmp = Intersection(pi,pj) :: Bool
                if tmp 
                    global doTheyIntersect = true
                    #println("They intersect")
                end
            end
            if doTheyIntersect break end
        end
        if doTheyIntersect break end
    end  

    return doTheyIntersect

end