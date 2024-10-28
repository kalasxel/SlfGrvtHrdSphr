using Distributions
using DelimitedFiles


include("ParticleStructure.jl")
include("Bounces.jl")
include("Initialize.jl")
include("Solvers.jl")


println("Start")


const workingDir = dirname(@__FILE__)
cd(workingDir)
foreach( rm, filter(endswith(".dat"), readdir(workingDir* "/OUTPUT/", join=true)) )


const halfLenght = 12.5992
const xi, xf = -halfLenght, halfLenght # area. better to keep as cube with equal sides
const yi, yf = -halfLenght, halfLenght
const zi, zf = -halfLenght, halfLenght
const dt = 2.e-6 # step
const tf = 1e2 # final time
const nbStps = Int(round( tf/dt, RoundDown )) # number of steps
const stFq = 100000 :: Int # how many steps between each output
const nb = 15 :: Int # number of particles
const M = 1e0 # all particles with the same mass
const R = 1e0 # all particles with the same radius
const G = 10e0 # dimensionless grav const. G -> G* M* t0^2 / L0^3 if all masses are the same
const σ = 1.8 # initial velocity dispersion
const nbTry = 900000 :: Int # max number of attempts for creating spatial distribution without intersections


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

SetAccelerations!(prtcls)


println("Total kinetic energy: ", TotKinEnergy(prtcls))
println("Total potential energy: ", TotPotEnergy(prtcls))
const EtotIni =  TotKinEnergy(prtcls) + TotPotEnergy(prtcls)


@time while(t<=tf)

    global t += dt
    global st += 1

    LeapfrogStep!(prtcls)
    AllBounces!(prtcls)

    if st % stFq == 0
        println( "Step ", st, "/", nbStps, ". ", "Total energy: ", TotKinEnergy(prtcls) + TotPotEnergy(prtcls), "." )
        @inbounds for i in 1:nb
            outputParticles[i,st÷stFq+1,:] = [ t, prtcls[i].r.x, prtcls[i].r.y, prtcls[i].r.z, KineticEnergy(prtcls[i]) ]
        end
    end

end


paramsOutput = [string(M) string(R)
                string(xi) string(xf)
                string(yi) string(yf)
                string(zi) string(zf)
                string(nb) string(Int(round(tf/dt/stFq+1)))
                string(G) string(σ)]

open( workingDir* "/OUTPUT/PARAMS.dat", "w" ) do io writedlm(io, paramsOutput) end
@inbounds for i in 1:nb
    open( workingDir* "/OUTPUT/particle"* string(i)* ".dat", "w" ) do io writedlm(io, outputParticles[i,:,:]) end
end


println("Total kinetic energy: ", TotKinEnergy(prtcls))
println("Total potential energy: ", TotPotEnergy(prtcls))

const EtotFin =  TotKinEnergy(prtcls) + TotPotEnergy(prtcls)

println("Full initial energy: ", EtotIni)
println("Full final energy: ", EtotFin)