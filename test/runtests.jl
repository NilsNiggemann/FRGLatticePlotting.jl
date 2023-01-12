import SpinFRGLattices as SL
import FRGLatticePlotting as FLP
using FRGLatticePlotting
using Test

@testset "SquareLattice" begin
    L1 = LatticeInfo(SL.SquareLattice.getSquareLattice(7),SL.SquareLattice)
    FLP.test_fourier_onsite(L1)
    FLP.test_fourier_pairs(L1)
end
@testset "Pyrochlore" begin
    L1 = LatticeInfo(SL.Pyrochlore.getPyrochlore(7),SL.Pyrochlore)
    FLP.test_fourier_onsite(L1)
    FLP.test_fourier_pairs(L1)
end

@testset "SquareKagome" begin
    L1 = LatticeInfo(SL.SquareKagome.getSquareKagome(7),SL.SquareKagome)
    FLP.test_fourier_onsite(L1)
    FLP.test_fourier_pairs(L1)
end

