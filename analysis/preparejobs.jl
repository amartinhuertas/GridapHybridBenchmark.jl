using Mustache
using DrWatson

jobname(args...) = replace(savename(args...;connector="_"),"="=>"_")
driverdir(args...) = normpath(projectdir("..",args...))

function convert_nc_np_to_prod(d)
  o=Dict()
  for k in keys(d)
    if k==:nc || k==:np
     o[k]=prod(d[k])
    else
     o[k]=d[k]
    end
  end
  return o
end

function jobdict(params)
  nr = params[:nr]
  np = params[:np]
  nc = params[:nc]
  fparams=convert_nc_np_to_prod(params)
  Dict(
  "q" => "normal",
  "o" => datadir(jobname(fparams,"o.txt")),
  "e" => datadir(jobname(fparams,"e.txt")),
  "walltime" => "00:30:00",
  "ncpus" => prod(np),
  "mem" => "$(prod(np)*4)gb",
  "name" => jobname(fparams),
  "nc" => nc,
  "numrefs" => haskey(params,:numrefs) ? params[:numrefs] : -1,
  "n" => prod(np),
  "np" => np,
  "nr" => nr,
  "projectdir" => driverdir(),
  "modules" => driverdir("modules.sh"),
  "title" => datadir(jobname(fparams)),
  "sysimage" => driverdir("GridaHybridBenchmark.so")
  )
end

function generate_2d_dicts(lst_nodes,lst_ls,nr=10)
   dicts = Dict[]
   d=2
   for node in lst_nodes
     px=6*node
     py=8*node
     for ls in lst_ls
        nx=px*ls
        ny=py*ls
        aux=Dict(:d=>2,
                  :nc=>(nx,ny),
                  :np=>(px,py),
                  :nr=>nr)
        push!(dicts,aux)
     end
   end
   dicts
end

function generate_3d_dicts(lst_nodes,lst_ls,nr=10)
  d=3
  for node in lst_nodes
    px=4*node
    py=4*node
    pz=3*node
    for ls in lst_ls
      nx=px*ls
      ny=py*ls
      nz=pz*ls
      aux=Dict(:d=>3,
                :nc=>(nx,ny,nz),
                :np=>(px,py,pz),
                :nr=>nr)
      push!(dicts,aux)
    end
  end
  dicts
end

dicts=generate_2d_dicts(collect(3:8),[16,32,64,128,256,512])
template = read(projectdir("jobtemplate.sh"),String)
for params in dicts
   fparams=convert_nc_np_to_prod(params)
   jobfile = datadir(jobname(fparams,"sh"))
   open(jobfile,"w") do io
     render(io,template,jobdict(params))
   end
end
