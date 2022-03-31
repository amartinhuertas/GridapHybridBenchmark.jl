module GridapHybridBenchmark

using Gridap
using GridapHybrid
using GridapDistributed
using GridapPETSc
using PartitionedArrays
using FileIO
using MPI
const PArrays=PartitionedArrays


function petsc_gamg_options()
  """
    -ksp_type cg -ksp_rtol 1.0e-06 -ksp_atol 0.0
    -ksp_monitor -pc_type gamg -pc_gamg_type agg -pc_gamg_est_ksp_type cg
    -mg_levels_esteig_ksp_type cg -mg_coarse_sub_pc_type cholesky
    -mg_coarse_sub_pc_factor_mat_ordering_type nd -pc_gamg_process_eq_limit 50
    -pc_gamg_square_graph 9 pc_gamg_agg_nsmooths 1
  """
end

u2(x) = VectorValue(1+x[1],1+x[2])
Gridap.divergence(::typeof(u2)) = (x) -> 2
p(x) = -3.14
∇p2(x) = VectorValue(0,0,0)
Gridap.∇(::typeof(p)) = ∇p
f2(x) = u2(x) + ∇p2(x)
# Normal component of u(x) on Neumann boundary
function g2(x)
  tol=1.0e-14
  if (abs(x[2])<tol)
    return -x[2] #-x[1]-x[2]
  elseif (abs(x[2]-1.0)<tol)
    return x[2] # x[1]+x[2]
  end
  Gridap.Helpers.@check false
end


#3D problem
u3(x) = VectorValue(1+x[1],1+x[2],1+x[3])
Gridap.divergence(::typeof(u3)) = (x) -> 3
∇p3(x) = VectorValue(0,0,0)
f3(x) = u3(x) + ∇p3(x)
function g3(x) # Normal component of u(x) on Neumann boundary
  @assert false
end

function ufg(D::Int)
   if (D==2)
    u2,f2,g2
   elseif (D==3)
    u3,f3,g3
   end
end

function dirichlet_tags(D::Int)
  if (D==2)
    collect(5:8)
  elseif (D==3)
    collect(21:26)
  end
end

function mytic!(t,comm)
  MPI.Barrier(comm)
  PArrays.tic!(t)
end

function _from_setup_fe_space_to_the_end(t,model,order=1)

  comm = model.models.comm

  D = num_cell_dims(model)
  Ω = Triangulation(ReferenceFE{D},model)
  Γ = Triangulation(ReferenceFE{D-1},model)
  ∂K = GridapHybrid.Skeleton(model)

  dtags=dirichlet_tags(D)
  u,f,_ = ufg(D)


  # FE formulation params
  τ = 1.0 # HDG stab parameter

  degree = 2*order+1
  dΩ     = Measure(Ω,degree)
  n      = get_cell_normal_vector(∂K)
  nₒ     = get_cell_owner_normal_vector(∂K)
  d∂K    = Measure(∂K,degree)

  # FESpaces
  mytic!(t,comm)
  # Reference FEs
  reffeᵤ = ReferenceFE(lagrangian,VectorValue{D,Float64},order;space=:P)
  reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
  reffeₗ = ReferenceFE(lagrangian,Float64,order;space=:P)

  # Define test FESpaces
  V = TestFESpace(Ω  , reffeᵤ; conformity=:L2)
  Q = TestFESpace(Ω  , reffeₚ; conformity=:L2)
  M = TestFESpace(Γ,
                  reffeₗ;
                  conformity=:L2,
                  dirichlet_tags=dtags)
  Y = MultiFieldFESpace([V,Q,M])

  # Define trial FEspaces
  U = TrialFESpace(V)
  P = TrialFESpace(Q)
  L = TrialFESpace(M,p)
  X = MultiFieldFESpace([U, P, L])
  PArrays.toc!(t,"FESpaces")

  # FE Affine Operator
  mytic!(t,comm)
  a((uh,ph,lh),(vh,qh,mh)) = ∫( vh⋅uh - (∇⋅vh)*ph - ∇(qh)⋅uh )dΩ +
                            ∫((vh⋅n)*lh)d∂K +
                            #∫(qh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
                            ∫(qh*(uh⋅n))d∂K +
                            ∫(τ*qh*ph*(n⋅nₒ))d∂K -
                            ∫(τ*qh*lh*(n⋅nₒ))d∂K +
                            #∫(mh*(uh⋅n+τ*(ph-lh)*n⋅no))*d∂K
                            ∫(mh*(uh⋅n))d∂K +
                            ∫(τ*mh*ph*(n⋅nₒ))d∂K -
                            ∫(τ*mh*lh*(n⋅nₒ))d∂K
  l((vh,qh,mh)) = ∫( vh⋅f + qh*(∇⋅u))*dΩ
  op=HybridAffineFEOperator((u,v)->(a(u,v),l(v)), X, Y, [1,2], [3])
  PArrays.toc!(t,"HybridAffineFEOperator")

  # Linear Solver
  mytic!(t,comm)
  solver = PETScLinearSolver()
  xh=solve(solver,op)
  PArrays.toc!(t,"Solve")

  # Error norms and print solution
  mytic!(t,comm)
  uh,_=xh
  e = u -uh
  e_l2=sqrt(sum(∫(e⋅e)dΩ))
  PArrays.toc!(t,"L2Norm")

  map_main(get_part_ids(model.models)) do part
    println("$(e_l2)\n")
  end

  ngdofs  = length(Y.gids)
  ngcells = length(model.gids)
  ngcells, ngdofs, e_l2
end


function generate_model_cartesian(parts,subdomains,partition)
  d = length(subdomains)
  domain = Vector{Float64}(undef, 2*d)
  for i = 1:2:2 * d
    domain[i]=0
    domain[i+1]=1
  end
  domain = Tuple(domain)
  model = CartesianDiscreteModel(parts, domain, partition)
end

function main_cartesian(parts,subdomains,partition,title,ir,order=1)
  t = PArrays.PTimer(parts,verbose=true)
  PArrays.tic!(t)
  model=generate_model_cartesian(parts,subdomains,partition)
  PArrays.toc!(t,"Model")
  ngcells, ngdofs, enorm = _from_setup_fe_space_to_the_end(t,model,order)
  display(t)
  nparts = length(parts)
  map_main(t.data) do data
    out = Dict{String,Any}()
    merge!(out,data)
    out["d"] = length(subdomains)
    out["enorm"] = enorm
    out["nparts"] = nparts
    out["ngdofs"] = ngdofs
    out["ngcells"] = ngcells
    out["ls"] = partition[1]/subdomains[1]
    out["nc"] = partition
    out["np"] = subdomains
    out["ir"] = ir
    save("$title.bson",out)
  end
end

########
function main(;
  np::Tuple,
  nr::Integer,
  title::AbstractString,
  nc::Tuple,
  numrefs::Integer=-1,
  k::Integer=1,
  verbose::Bool=true)

  # Process parameters of mesh
  length(np) == length(nc) || throw(ArgumentError("np and nc must be of same length"))
  all(nc .> 0) || throw(ArgumentError("all values in nc should be larger than 0"))

  options=petsc_gamg_options()

  prun(mpi,np) do parts
    for ir in 1:nr
      GridapPETSc.with(args=split(options)) do
          str_r   = lpad(ir,ceil(Int,log10(nr)),'0')
          title_r = "$(title)_ir$(str_r)"
          main_cartesian(parts,np,nc,title_r,ir,1)
          GridapPETSc.gridap_petsc_gc()
      end
    end
  end
end

#main(;np=(1,1),nr=1,title="test",nc=(10,10))
main(;np=(1,1,1),nr=1,title="test",nc=(10,10,10))

end
