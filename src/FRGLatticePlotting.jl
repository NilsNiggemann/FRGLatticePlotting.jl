module FRGLatticePlotting
    # using SpinFRGLattices,Parameters,StaticArrays,LaTeXStrings,Plots
    using GLMakie, SpinFRGLattices, MakieHelpers,LaTeXStrings 
    include("LatticePlot.jl")
    export pairsPlot, plotSystem, plotBonds!,plotDistBonds!,plotDistBonds

    GLMakie.activate!()
    GLMakie.Makie.inline!(false)
    using PrecompileTools
    include("precompile.jl")

end # module
