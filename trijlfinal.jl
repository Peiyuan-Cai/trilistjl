using ITensors
using Random
Random.seed!(39)
pushfirst!(LOAD_PATH, "~/trans/triangularjl/")
include("struct_qshlistit.jl")


let
  Nx = 4
  Ny = 4

  N = Nx*Ny

  sites = siteinds("S=1/2", N)
  scale = 1.

  #interactions
  Jxy = 1. * scale
  Jz = 2. * scale
  Jc = 2. * scale
  Jcz = 1. * scale

  hl = Hlist("torus", Nx, Ny)
  os = OpSum()
  #xxz terms
  for k in 1:length(hl.hlist)
    os .+= Jxy/2, "S+", hl.hlist[k][1], "S-", hl.hlist[k][2]
    os .+= Jxy/2, "S-", hl.hlist[k][1], "S+", hl.hlist[k][2]
    os .+= Jz, "Sz", hl.hlist[k][1], "Sz", hl.hlist[k][2]
  end
  #compass and in-out
  for ktype in 1:3
    type = ktype-1
    if type == 0
      nx, ny, nz = 1, 0, 0
    elseif type == 1
      nx, ny, nz = 1/2, sqrt(3)/2, 0
    else type == 2
      nx, ny, nz = -1/2, sqrt(3)/2, 0
    end
    for k in 1:length(hl.sortedhlist[ktype])
      #compass
      os .+= 4*Jc*((nx*nx)-0.5), "Sx", hl.sortedhlist[ktype][k][1], "Sx", hl.sortedhlist[ktype][k][2]
      os .+= 4*Jc*nx*ny, "Sx", hl.sortedhlist[ktype][k][1], "Sy", hl.sortedhlist[ktype][k][2]
      os .+= 4*Jc*nx*ny, "Sy", hl.sortedhlist[ktype][k][1], "Sx", hl.sortedhlist[ktype][k][2]
      os .+= 4*Jc*((ny*ny)-0.5), "Sy", hl.sortedhlist[ktype][k][1], "Sy", hl.sortedhlist[ktype][k][2]
      #in-out
      os .+= 2*Jcz*ny, "Sx", hl.sortedhlist[ktype][k][1], "Sz", hl.sortedhlist[ktype][k][2]
      os .+= -2*Jcz*nx, "Sy", hl.sortedhlist[ktype][k][1], "Sz", hl.sortedhlist[ktype][k][2]
      os .+= 2*Jcz*ny, "Sz", hl.sortedhlist[ktype][k][1], "Sx", hl.sortedhlist[ktype][k][2]
      os .+= -2*Jcz*nx, "Sz", hl.sortedhlist[ktype][k][1], "Sy", hl.sortedhlist[ktype][k][2]
    end
  end

  H = MPO(os,sites)

  state = [isodd(n) ? "Up" : "Dn" for n=1:N]
  # Initialize wavefunction to a random MPS
  # of bond-dimension 10 with same quantum
  # numbers as `state`
  psi0 = randomMPS(sites,state,20)

  nsweeps = 20
  maxdim = [20,100,200,400]
  cutoff = [1E-9]

  energy,psi = dmrg(H,psi0; nsweeps, maxdim, cutoff)
  @show energy

  czz = correlation_matrix(psi, "Sz", "Sz")
  println("Scorr(1,2) is ", czz[1,2])
  println("Scorr(2,4) is ", czz[2,4])

  mz = expect(psi, "Sz")
  @show mz

  mx = expect(psi, "Sx")
  @show mx

  my = expect(psi, "Sy")
  @show my

  one = norm(psi)
  @show one

  return energy
end