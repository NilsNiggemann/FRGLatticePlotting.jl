module FRGLatticePlotting

using CairoMakie, SpinFRGLattices, MakieHelpers, LaTeXStrings

HASGSL = false
try
    using GSL
    HASGSL = true
catch
    @info "GSL not found, using fallback"
end

include("LatticePlot.jl")
export pairsPlot,pairsPlot!, plotSystem, plotSystem!, plotBonds!, plotDistBonds!, plotDistBonds, scatterRvec!,linesRvec!, scatterRvec, linesRvec, getStandardFigure

using PrecompileTools
include("precompile.jl")

function __init__()
    if HASGSL
        GLMakie.activate!(inline=false)
    end
end

end # module
