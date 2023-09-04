module FRGLatticePlotting

using GLMakie, CairoMakie, SpinFRGLattices, MakieHelpers, LaTeXStrings
include("LatticePlot.jl")
export pairsPlot,pairsPlot!, plotSystem, plotSystem!, plotBonds!, plotDistBonds!, plotDistBonds, scatterRvec!,linesRvec!, scatterRvec, linesRvec, getStandardFigure

using PrecompileTools
include("precompile.jl")

function __init__()
    GLMakie.activate!(inline=false)
end

end # module
