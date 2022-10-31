module FRGLatticePlotting
    using SpinFRGLattices,Parameters,StaticArrays,LaTeXStrings,Plots

    include("LatticePlot.jl")
    
    export pairsPlot, plotSystem, plotCouplings!,plotCorrelations!,plotBond!,plotBonds!,plotDistBonds!,plotDistBonds,getVertexR

end # module
