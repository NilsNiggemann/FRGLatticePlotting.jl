module FRGLatticePlotting
    using SpinFRGLattices,Parameters,StaticArrays,LaTeXStrings,Plots

    include("Fourier.jl")
    export AbstractLattice, LatticeInfo,FourierTransform, Fourier2D, equalTimeChiBeta, EnergyBeta, get_e_Chi, Chikplot, getFlow, plotFlow, plotMaxFlow,plotMaxFlow!,plotMaxFlow_fast, pointPath, fetchKPath, plotKpath,plotKpath!, pscatter!,pplot!,getkMax,Fourier3D

    export hhlplane,xyplane,zzerocut,sphereplane
    
    include("LatticePlot.jl")
    
    export pairsPlot, plotSystem, plotCouplings!,plotCorrelations!,plotBond!,plotBonds!,plotDistBonds!,plotDistBonds
    
    
    include("Tests/FourierTest.jl")
    export test_fourier_onsite,test_fourier_pairs

    
end # module
