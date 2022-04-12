module FRGLatticePlotting
    using SpinFRGLattices,Parameters,StaticArrays,LaTeXStrings,Plots

    using SLEEFPirates:cos_fast

    include("Fourier.jl")
    export AbstractLattice, LatticeInfo,FourierTransform, Fourier2D, equalTimeChiBeta, EnergyBeta, get_e_Chi, Chikplot, getFlow, plotFlow, getMaxFlow, plotMaxFlow,plotMaxFlow!,plotMaxFlow_fast, pointPath, fetchKPath, plotKpath,plotKpath!, pscatter!,pplot!,getkMax,Fourier3D

    
    include("LatticePlot.jl")
    
    export pairsPlot, plotSystem, plotCouplings!,plotCorrelations!,plotBond!,plotBonds!,plotDistBonds!,plotDistBonds,getVertexR
    
    include("helpers.jl")
    export hhlplane,xyplane,zzerocut,sphereplane,hhllabels,strd,strd3,getPairs,DimerResponse,allOccurIn,findNames,NameFilter,OnlyIndex,getNumberFromName,stringLatex
    

    include("Tests/FourierTest.jl")
    export test_fourier_onsite,test_fourier_pairs

    
end # module
