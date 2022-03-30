using PackageCompiler
create_sysimage(:GridapHybridBenchmark,
  sysimage_path=joinpath(@__DIR__,"..","GridapHybridBenchmark.so"),
  precompile_execution_file=joinpath(@__DIR__,"..","src","GridapHybridBenchmark.jl"))
