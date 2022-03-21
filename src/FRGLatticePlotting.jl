module FRGLatticePlotting
    using SpinFRGLattices,Parameters,StaticArrays,LaTeXStrings,Plots

    include("Fourier.jl")
    include("LatticePlot.jl")
    include("Tests/FourierTest.jl")

    export test_fourier_onsite,test_fourier_pairs

end # module
