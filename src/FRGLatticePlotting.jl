module FRGLatticePlotting
    # using SpinFRGLattices,Parameters,StaticArrays,LaTeXStrings,Plots
    using GLMakie, SpinFRGLattices, MakieHelpers #,PrecompileTools
    include("LatticePlot.jl")
    export pairsPlot, plotSystem, plotBonds!,plotDistBonds!,plotDistBonds

end # module
