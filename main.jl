using Distributions
using DelimitedFiles


include("ParticleStructure.jl")
include("Bounces.jl")
include("Initialize.jl")


println("Start")


const workingDir = dirname(@__FILE__)
cd(workingDir)
foreach( rm, filter(endswith(".dat"), readdir(workingDir* "/OUTPUT/", join=true)) )


const halfLenght = 40e0
const xi, xf = -halfLenght, halfLenght # area. better to keep as cube with equal sides
const yi, yf = -halfLenght, halfLenght
const zi, zf = -halfLenght, halfLenght
const dt = .001 # step
const tf = 3e2 # final time
const nbStps = Int(round( tf/dt, RoundDown )) # number of steps
const stFq = 1000 :: Int # how many steps between each output
const nb = 20:: Int # number of particles
const M = 1e0 # all particles with the same mass
const R = 3e0 # all particles with the same radius
const G = 10 # dimensionless grav const. G -> G* M* t0^2 / L0^3 since all masses are the same
const σ = 1e0 # initial velocity dispersion
const nbTry = 900000 # max number of attempts for creating spatial distribution without intersections


t = 0e0 :: Float64
st = 0 :: Int64

prtcls = [ Particle{Float64}( M, R, CartesianVector{Float64}( [xf+R+1e0, xf+R+1e0, xf+R+1e0] ), 
                                CartesianVector{Float64}( rand(Normal(0,σ), 3) ),
                                CartesianVector{Float64}( [0e0, 0e0, 0e0] ) ) for i in 1:nb ] # generate nb particles outside box with normally distribute velocities

countInit = 1 :: Int64
while makeInitialSpatialDistribution!(prtcls) 
    global countInit += 1
    if countInit >= nbTry
        break
    end
end

if countInit >= nbTry
    error("After ", nbTry, " attemps cannot place all ", nb, " particles. Try again or change parameters.")
end
println(countInit, " from ", nbTry, " attemps were made to place ", nb, " particles without intersections.")                     

outputParticles = zeros( nb, Int(round( tf/(stFq*dt), RoundDown ))+1, 5 ) # output initial
@inbounds for i in 1:nb
    outputParticles[i,1,:] = [t, prtcls[i].r.x, prtcls[i].r.y, prtcls[i].r.z, KineticEnergy(prtcls[i])]
end


println("Total kinetic energy: ", TotKinEnergy(prtcls))
println("Total potential energy: ", TotPotEnergy(prtcls))
const EtotIni =  TotKinEnergy(prtcls) + TotPotEnergy(prtcls)

@time while(t<=tf)

    global t += dt
    global st += 1

    for pi in prtcls
        acceleration = CartesianVector{Float64}( [0e0, 0e0, 0e0] );
        for pj in prtcls
            if pi != pj
                acceleration = acceleration + (G/AbsVec( pj.r - pi.r )^3) * (pj.r - pi.r)
            end
        end
        pi.a = acceleration
    end

    for pi in prtcls
        pi.v = pi.v + pi.a * dt
        pi.r = pi.r + pi.v * dt
    end

    for pi in prtcls
        BounceBoundaries!(pi)
        for pj in prtcls
            if pi != pj
                BounceParticles!(pi,pj)
            end
        end
    end

    if st % stFq == 0
        println( "Step ", st, "/", nbStps, ". ", "Total energy: ", TotKinEnergy(prtcls) + TotPotEnergy(prtcls) )
        @inbounds for i in 1:nb
            outputParticles[i,st÷stFq+1,:] = [t, prtcls[i].r.x, prtcls[i].r.y, prtcls[i].r.z, KineticEnergy(prtcls[i])]
        end
    end

end


paramsOutput = [string(M) string(R)
                string(xi) string(xf)
                string(yi) string(yf)
                string(zi) string(zf)
                string(nb) string(Int(round(tf/dt/stFq+1)))]

open( workingDir* "/OUTPUT/PARAMS.dat", "w" ) do io writedlm(io, paramsOutput) end
@inbounds for i in 1:nb
    open( workingDir* "/OUTPUT/particle"* string(i)* ".dat", "w" ) do io writedlm(io, outputParticles[i,:,:]) end
end


println("Total kinetic energy: ", TotKinEnergy(prtcls))
println("Total potential energy: ", TotPotEnergy(prtcls))

const EtotFin =  TotKinEnergy(prtcls) + TotPotEnergy(prtcls)

println("Full initial energy: ", EtotIni)
println("Full final energy: ", EtotFin)