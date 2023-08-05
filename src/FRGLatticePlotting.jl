module FRGLatticePlotting
    # using SpinFRGLattices,Parameters,StaticArrays,LaTeXStrings,Plots
    using GLMakie, SpinFRGLattices, MakieHelpers,LaTeXStrings #,PrecompileTools
    include("LatticePlot.jl")
    export pairsPlot, plotSystem, plotBonds!,plotDistBonds!,plotDistBonds

end # module
