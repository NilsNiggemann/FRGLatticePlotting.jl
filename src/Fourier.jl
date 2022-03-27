export AbstractLattice, LatticeInfo,FourierTransform, Fourier2D, equalTimeChiBeta, EnergyBeta, get_e_Chi, Chikplot, getFlow, plotFlow, plotMaxFlow,plotMaxFlow!,plotMaxFlow_fast, pointPath, fetchKPath, plotKpath,plotKpath!, pscatter!,pplot!,getkMax,Fourier3D

struct FourierInfo{Dim}
    pairs::Vector{Int}
    Rij_vec::Vector{SVector{Dim,Float64}}
end
abstract type AbstractLattice end
@with_kw struct LatticeInfo{BasisType,RvecType,FunctionType,Dim} <: AbstractLattice
    System::Geometry
    Basis::BasisType
    NLen::Int = System.NLen
    Npairs::Int = System.Npairs
    NUnique::Int = System.NUnique
    PairList::Vector{RvecType} = System.PairList
    PairTypes::Vector{sitePair} = System.PairTypes
    SiteList::Vector{RvecType} = unique(SpinFRGLattices.sortedPairList(NLen,Basis)[1])
    UnitCell::Vector{RvecType} = [SpinFRGLattices.getRvec(b,Basis) for b in Basis.b]
    pairToInequiv::FunctionType
    FourierInfos::FourierInfo{Dim} = PrecomputeFourier(UnitCell,SiteList,PairList,PairTypes,pairToInequiv,Basis)
end
LatticeInfo(System,Mod::Module,PTI =Mod.pairToInequiv ) = LatticeInfo(System=System,Basis = Mod.Basis,pairToInequiv = PTI)

getDim(B::Basis_Struct_2D) = 2
getDim(B::Basis_Struct_3D) = 3

function PrecomputeFourier(UnitCell::AbstractVector{T},SiteList::AbstractVector{T},PairList::AbstractVector{T},PairTypes,pairToInequiv,Basis) where T <: Rvec
    Rij_vec = SVector{getDim(Basis),Float64}[]
    pairs = Int[]
    for i_site in UnitCell
        Ri = getCartesian(i_site,Basis)
        # println(Ri)
        for j_site in SiteList # site summation
            Rj = getCartesian(j_site,Basis)
            Rij = Ri - Rj
            R_Ref,ij = pairToInequiv(i_site,j_site) #Map j to correct pair so that we may use Chi_0,j'
            xi = getSiteType(R_Ref,Basis)
            pair = MapToPair(xi,ij,PairList,PairTypes)
            if pair !== 0
                push!(Rij_vec,Rij)
                push!(pairs,pair)
            end
        end
    end
    return FourierInfo(pairs,Rij_vec)
end

# function FourierTransform(k::StaticVector,Chi_R, NCell,pairs,Rij_vec)
#     Chi_k = 0. +0im
#     for (p,Rij) in zip(pairs,Rij_vec)
#         Chi_k += exp(1im * k' * Rij) * Chi_R[p]
#     end
#     return 1/NCell * real(Chi_k)
# end

# function FourierTransform(k::StaticVector,Chi_R, NCell,pairs,Rij_vec)
#     Chi_k = 0.
#     for (p,Rij) in zip(pairs,Rij_vec)
#         Chi_k += cos(k' * Rij) * Chi_R[p]
#     end
#     return 1/NCell * Chi_k
# end

@inline function dot(v1::AbstractVector{T},v2::AbstractVector{T}) where T
    res = zero(T)
    @inbounds @simd for i in eachindex(v1,v2)
        res+= v1[i]*v2[i]
    end
    return res
end

@inline function FourierTransform(k::AbstractVector,Chi_R, NCell,pairs,Rij_vec)
    Chi_k = 0.
    # @inbounds @simd for (p,Rij) in zip(pairs,Rij_vec)
    @inbounds @simd for i in eachindex(pairs,Rij_vec)
        p = pairs[i]
        Rij = Rij_vec[i]
        Chi_k += cos(dot(k, Rij)) * Chi_R[p]
    end
    return 1/NCell * Chi_k
end

@inline FourierTransform(k,Chi_R, Lattice::AbstractLattice) = FourierTransform(k,Chi_R, Lattice.Basis.NCell,Lattice.FourierInfos.pairs,Lattice.FourierInfos.Rij_vec) 

"""Returns 2D Fourier trafo in plane as specified by the "regionfunc" function. Eg for a plot in the xy plane we can use plane = (ki,kj) -> SA[ki,kj] """
function Fourier2D(Chi_R::AbstractArray,regionfunc::Function,Lattice::AbstractLattice;res=100,ext = pi,minext = -ext)
    karray = range(minext,stop = ext,length = res)
    Chi_k = Fourier2D(Chi_R::AbstractArray,karray,karray ,regionfunc::Function,Lattice::AbstractLattice)
    return karray,Chi_k
end

function Fourier2D(Chi_R::AbstractArray,x::AbstractVector,y::AbstractVector, regionfunc::Function,Lattice::AbstractLattice)
    Chi_k = zeros(length(x),length(y))

    Threads.@threads for j in eachindex(y)
        kj = x[j]
        for (i,ki) in enumerate(x)
            Chi_k[i,j] = FourierTransform(regionfunc(ki,kj),Chi_R,Lattice)
        end
    end
    return Chi_k
end

"""Returns 3D Fourier trafo"""
function Fourier3D(Chi_R::AbstractArray,Lattice::AbstractLattice;res=50,ext = pi,minext = -ext)
    k = range(minext,stop = ext,length = res)
    k, Fourier3D(Chi_R,Lattice,k,k,k)
end

function Fourier3D(Chi_R::AbstractArray,Lattice::AbstractLattice,kx_vec::AbstractVector,ky_vec::AbstractVector,kz_vec::AbstractVector)
    Chi_k = zeros(length(kx_vec),length(ky_vec),length(kz_vec))

    Threads.@threads for iz in eachindex(kz_vec)
        kz = kz_vec[iz]
        for (iy,ky) in enumerate(ky_vec),(ix,kx) in enumerate(kx_vec)
            Chi_k[ix,iy,iz] = FourierTransform(SA[kx,ky,kz],Chi_R,Lattice)
        end
    end
    return Chi_k
end
##

"""Computes χ_ij(τ=0)/T"""
function equalTimeChiBeta(Chi_RNu,N_nu = size(Chi_RNu,2))
    # Npairs,N_nu = size(Chi_RNu)
    Chi_Tau0 = Chi_RNu[:,begin] # add static nu=0 component, appears only once in sum
    for R in eachindex(Chi_Tau0)
        sum = 0.
        for n in 2:N_nu # dynamic components
            sum += Chi_RNu[R,n]
        end
        Chi_Tau0[R] += 2*sum# Chi(nu=0)
    end
    return Chi_Tau0
end

"""Compute energy divided by temperature from spin correlations"""
function EnergyBeta(Chi_RNu, Lattice,Nnu)
    @unpack PairList,SiteList,PairTypes,Basis,UnitCell,pairToInequiv,Npairs = Lattice
    J_ij = Lattice.System.couplings
    E = 0.
    Chi_Tau0 = equalTimeChiBeta(Chi_RNu,Nnu)

    for i_site in UnitCell
        # println(Ri)
        for j_site in SiteList # site summation
            R_Ref,ij = pairToInequiv(i_site,j_site) #Map j to correct pair so that we may use Chi_0,j'
            xi = getSiteType(R_Ref,Basis)
            pair = MapToPair(xi,ij,PairList,PairTypes)
            if pair !== 0
                E += 3/(2*Basis.NCell) *J_ij[pair] * Chi_Tau0[pair]
            end
            # println(j_site,Chi_R[pair])
        end
    end
    return E
end

function get_e_Chi(Chi_TRnu,Trange,Lattice,Nnu)
    e_Chi = similar(Trange)
    for (iT,T) in enumerate(Trange)
        e_Chi[iT] = @views T*EnergyBeta(Chi_TRnu[iT,:,:],Lattice,Nnu)
    end
    return e_Chi
end


##
pitick(x) = latexstring("$(round(Int,x/pi)) \\pi") 
PiMultipleTicks(ticks) = [pitick(x) for x in ticks  ]


function pitick(start, stop, denom; mode=:text)
    #taken from https://discourse.julialang.org/t/plot-with-x-axis-in-radians-with-ticks-at-multiples-of-pi/65325/3
    a = Int(cld(start, π/denom))
    b = Int(fld(stop, π/denom))
    tick = range(a*π/denom, b*π/denom; step=π/denom)
    ticklabel = piticklabel.((a:b) .// denom, Val(mode))
    tick, ticklabel
end

function piticklabel(x::Rational, ::Val{:text})
    #taken from https://discourse.julialang.org/t/plot-with-x-axis-in-radians-with-ticks-at-multiples-of-pi/65325/3
    iszero(x) && return "0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return S * N * "π"
    S * N * "π/" * repr(d)
end

function piticklabel(x::Rational, ::Val{:latex})
    #taken from https://discourse.julialang.org/t/plot-with-x-axis-in-radians-with-ticks-at-multiples-of-pi/65325/3
    iszero(x) && return L"0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return L"%$S%$N\pi"
    L"%$S\frac{%$N\pi}{%$d}"
end


function Chikplot(kx::AbstractArray,ky::AbstractArray,Chi_k;xlabel = L"k_x",ylabel= L"k_y",colorscheme = :viridis,tickfontsize = 14,labelfontsize = 20,PiStep = 1, kwargs...)
    denom = 1//PiStep
    piticks(ax) = pitick(minimum(ax),maximum(ax),denom)

    pl = heatmap(kx,ky,transpose(Chi_k),size = (570, 600),xticks=piticks(kx),yticks=piticks(ky),linewidths=0.0,xlabel=xlabel ,ylabel= ylabel,c= colorscheme,right_margin = 15 *Plots.px,aspectratio = 1,tickfontsize = tickfontsize,labelfontsize = labelfontsize;kwargs...)
    return pl
end

Chikplot(k,Chi_k; kwargs...) = Chikplot(k,k,Chi_k; kwargs...)
##
function getFlow(k::StaticArray,Chi_LR,Lambdas,Lattice)
    flow = similar(Lambdas)
    FT(k,Chi) = FourierTransform(k,Chi,Lattice)
    for i in eachindex(Lambdas)
        flow[i] =  @views FT(k,Chi_LR[i,:])
    end
    return flow
end

function plotFlow(k::StaticArray,Chi_LR,Lambdas,Lattice,pl = plot();method = plot!,xmax=1.,kwargs...)
    flow = getFlow(k,Chi_LR,Lambdas,Lattice)
    method(pl,Lambdas,flow, xlims = (0.,xmax);kwargs...)
    return pl
end 

function plotMaxFlow(Chi_LR,Lambdas,Lattice,regionfunc::Function,pl = plot();  res = 30,ext = pi,xmax=1.,method = plot!,kwargs...)
    flow = similar(Lambdas)
    FT(Chi) = maximum(Fourier2D(Chi,regionfunc,Lattice,res = res,ext=ext)[2])
    for i in eachindex(Lambdas)
        flow[i] =  @views FT(Chi_LR[i,:])
    end
    method(pl,Lambdas,flow, xlims = (0.,xmax);kwargs...)
    return pl
end

function plotMaxFlow(Chi_LR,Lambdas,Lattice::LatticeInfo{Basis_Struct_3D,Rvec_3D,F,3},pl = plot();  res = 30,ext = pi,xmax=1.,method = plot!,kwargs...) where {F}
    flow = similar(Lambdas)
    FT(Chi) = maximum(Fourier3D(Chi,Lattice,res = res,ext=ext)[2])
    for i in eachindex(Lambdas)
        flow[i] =  @views FT(Chi_LR[i,:])
    end
    method(pl,Lambdas,flow, xlims = (0.,xmax);kwargs...)
    return pl
end

plotMaxFlow!(pl,Chi_LR,Lambdas,Lattice,regionfunc::Function;kwargs...) = plotMaxFlow(Chi_LR,Lambdas,Lattice,regionfunc,pl;kwargs...)

plotMaxFlow!(Chi_LR,Lambdas,Lattice,regionfunc::Function;kwargs...) = plotMaxFlow(Chi_LR,Lambdas,Lattice,regionfunc,current();kwargs...)

function getkMax(k::AbstractVector,Chik::AbstractMatrix)
    ik1,ik2 = Tuple(argmax(Chik))
    return k[ik1],k[ik2]
end

function getkMax(Chi_R,Lattice::LatticeInfo,regionfunc::Function;kwargs...)
    k1,k2 = getkMax(Fourier2D(Chi_R,regionfunc,Lattice;kwargs...)...)
    return regionfunc(k1,k2)
end

function plotMaxFlow_fast(Chi_LR,Lambdas,Lattice,regionfunc::Function,pl = plot();  res = 90,ext = pi,xmax=1.,method = plot!,kwargs...)
    flow = similar(Lambdas)
    flowEnd = argmin(Lambdas)
    chiend = Chi_LR[flowEnd,:]
    getkMax(chiend,Lattice,regionfunc;res=res,ext = ext)
    for i in eachindex(Lambdas)
        flow[i] =  @views FourierTransform(point,Chi_LR[i,:],Lattice)
    end
    method(pl,Lambdas,flow, xlims = (0.,xmax);kwargs...)
    return pl
end  
##
function pointPath(p1::StaticArray,p2::StaticArray,res)
    Path = Vector{typeof(p1)}(undef,res)
    for i in eachindex(Path)
        Path[i] = p1 + i/res*(p2 -p1)
    end
    return Path
end
"""res contains the number of points along -pi,pi"""
function fetchKPath(points,res = 100)
    Path = Vector{typeof(points[begin])}(undef,0)
    # Path = []
    PointIndices = [1]
    for i in eachindex(points[begin:end-1])
        p1 = points[i]
        p2 = points[i+1]
        append!(Path,pointPath(p1,p2,round(Int,norm(p1-p2)/2pi * res)))
        append!(PointIndices,length(Path)) # get indices corresponding to points
    end
    return PointIndices,Path
end
##

function plotKpath!(pl,Chi_R,Lattice,points; PointTicks = [round.(p,digits = 2) for p in points],res = 100,xmax=1.,method = plot!,kwargs...)

    # PointTicks(ticks) = [p for p in points ]
    
    FT(k) = FourierTransform(k,Chi_R,Lattice)
    ticks,Path = fetchKPath(points,res)
    FTPath = FT.(Path)
    p0 = first(Path)
    pend = last(Path)
    lpath = norm(p0-pend)
    xarr = lpath.* LinRange(0,1,length(Path))
    method(pl,xarr,FTPath,xticks = (xarr[ticks],PointTicks);kwargs...)
    return pl
end 

plotKpath!(args...;kwargs...) = plotKpath!(current(),args...;kwargs...)
plotKpath(args...;kwargs...) = plotKpath!(plot(),args...;kwargs...)
##
# pscatter!(pArray;kwargs...) = scatter!(Tuple(p for p in pArray);kwargs...)
# pscatter(pArray;kwargs...) = scatter(Tuple(p for p in pArray);kwargs...)

pscatter!(pArray;kwargs...) = scatter!(Tuple([p[i] for p in pArray] for i in 1:length(pArray[1]));kwargs...)
pscatter(pArray;kwargs...) = scatter(Tuple([p[i] for p in pArray] for i in 1:length(pArray[1]));kwargs...)
pplot!(pArray;kwargs...) = plot!(Tuple(p for p in pArray);kwargs...)
