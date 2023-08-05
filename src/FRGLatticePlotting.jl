module FRGLatticePlotting
# using SpinFRGLattices,Parameters,StaticArrays,LaTeXStrings,Plots
using GLMakie, CairoMakie, SpinFRGLattices, MakieHelpers, LaTeXStrings
include("LatticePlot.jl")
export pairsPlot, plotSystem, plotBonds!, plotDistBonds!, plotDistBonds


using PrecompileTools
include("precompile.jl")

function __init__()
    GLMakie.activate!(inline=false)
end

end # module
