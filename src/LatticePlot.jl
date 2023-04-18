
"""Plots sites in PairList."""
function pairsPlot(PairList,Basis,pl = plot(size = (700,700),aspectratio = 1);colors = ("blue","red","black","cyan","yellow","green","pink","orange","lime","brown","grey"),color = "",colorBasis = false,kwargs...)
    uniquepairs = unique(PairList)
    
    if color != ""
        colors = Tuple(color for i in 1:Basis.NCell)
    end
    prepdata(x) = Tuple([i] for i in x)
    for R in uniquepairs
        if colorBasis # color each site differently according to basis
            col = colors[R.b]
        else 
            col = colors[Basis.SiteType[R.b]]
        end
        points = prepdata(getCartesian(R,Basis))
        scatter!(pl,points,label = nothing,color = col;kwargs...)
    end
    return pl
end

"""Plot all sites and inequivalent pairs"""
function plotSystem(System,Basis;
    plotAll = true,
    refSite = nothing,
    markersize = 5,
    inequivColor = "green",
    inequivalpha = 0.5,
    plotBonds=true,
    plotCouplings=true,
    CouplingColors = nothing,
    bondlw = 7,
    Bonds = [(minDist = Basis.NNdist-1e-3,maxDist = Basis.NNdist+1e-3,colorRGB = [0,0,0])],
    allpairs = unique!(SpinFRGLattices.sortedPairList(System.NLen,Basis)[1]),
    kwargs...)
    (;PairList,OnsitePairs )= System
    
    indices = copy(OnsitePairs)
    push!(indices,length(PairList)) # get final index
    if refSite === nothing 
        plotpairs = unique(PairList)
    else
        # allpairs = unique!(generatePairSites(System.NLen,Basis,Basis.refSites[refSite]))
        plotpairs = PairList[indices[refSite]:indices[refSite+1]]
    end
    filter!(x-> x in allpairs,plotpairs)

    plotAll || (allpairs = plotpairs)
    pl = pairsPlot(allpairs,Basis,markersize = markersize,aspect_ratio=:equal;kwargs...)
    if plotBonds
        for b in Bonds
            plotDistBonds!(allpairs,Basis,minDist = b.minDist, maxDist = b.maxDist,lw = bondlw,color = Plots.Colors.RGB((b.colorRGB./255)...))
        end
    end
    # plotBonds && plotDistBonds!(allpairs,Basis;color = Bondcolor,lw = bondlw, minDist = bondDist-1e-3, maxDist = bondDist+1e-3)

    plotAll && pairsPlot(plotpairs,Basis,pl,color = inequivColor,alpha = inequivalpha,markersize = 2*markersize)
    pairsPlot([Basis.refSites[refSite]],Basis,pl,color = "darkred",markershape = :cross,markersize = 1.5*markersize)

    plotCouplings && plotCouplings!(System,Basis;refSite = refSite,colors = CouplingColors)
    return pl
end
function plotCorrelations!(System::Geometry,Basis::Basis_Struct,couplings::AbstractVector,pl::Plots.Plot=current();refSite = nothing,colors = nothing,kwargs...)

    @unpack PairTypes,PairList = System
    inds = findall(x-> abs(x)>1E-14,couplings)

    for (i,Ind) in enumerate(inds)
    type,R_pair,J = PairTypes[Ind],PairList[Ind],couplings[Ind]

        refSite in (nothing,type.xi) || continue

        R_Ref = Basis.refSites[type.xi]
        plotCoup!(col::Nothing) = plotBond!(R_Ref,R_pair,Basis,pl,label = round(J,digits=3),lw = 1+4*abs(J);kwargs...)
        plotCoup!(col::Function) = plotBond!(R_Ref,R_pair,Basis,pl,label = round(J,digits=3),lw = 1+4*abs(J);kwargs...,col(J)...)
        plotCoup!(col) = plotBond!(R_Ref,R_pair,Basis,pl,label = round(J,digits=3),lw = 1+4*abs(J);kwargs...,color = colors[i])

        plotCoup!(colors)
    end
    return pl
end

@inline dim(B::Basis_Struct_2D) = 2
@inline dim(B::Basis_Struct_3D) = 3

"""Returns list of sites within bond length"""
function getBonds(R::RType,minDist::Real,maxDist::Real,Basis::Basis_Struct) where RType <: Rvec
    NN = RType[]
    getLConst(x) = getLatticeConst(x,Basis)
    maxN = div(maxDist,minimum(getLConst,1:dim(Basis)))+1

    AllSites = generateLUnitCells(maxN,Basis,R)
    for Rnew in AllSites
        d = dist(R,Rnew,Basis)
        if (isapprox(d,maxDist,atol = 1E-10) || d < maxDist) && d > minDist
            push!(NN,Rnew)
        end
    end
    return NN
end

getBonds(R::Rvec,maxDist::Real,Basis::Basis_Struct) = getBonds(R,0.,maxDist,Basis)

function getLatticeConst(direction,Basis::Basis_Struct_3D)
    if direction == 1
        R2 = Rvec(1,0,0,1)
    elseif direction == 2
        R2 = Rvec(0,1,0,1)
    elseif direction == 3
        R2 = Rvec(0,0,1,1)
    else 
        error("direction must be <3")
    end
    dist(Rvec(0,0,0,1),R2,Basis)
end

function getLatticeConst(direction,Basis::Basis_Struct_2D)
    if direction == 1
        R2 = Rvec(1,0,1)
    elseif direction == 2
        R2 = Rvec(0,1,1)
    else 
        error("direction must be <2")
    end
    dist(Rvec(0,0,1),R2,Basis)
end


function plotCouplings!(System,Basis,pl=current();kwargs...)
    plotCorrelations!(System,Basis,System.couplings,pl;kwargs...)
end

function plotBond!(R1::Rvec,R2::Rvec,Basis::Basis_Struct,pl=current();kwargs...)
    r1 = getCartesian(R1,Basis)
    r2 = getCartesian(R2,Basis)
    points = Tuple(SA[r1[i],r2[i]] for i in eachindex(r1))
    plot!(pl,points...;kwargs...)
end

function plotBonds!(Site::Rvec,Sites::AbstractVector{T},Basis::Basis_Struct,pl = current();kwargs...) where T <: Rvec
    for R2 in Sites
        plotBond!(Site,R2,Basis,pl;kwargs...)
    end
    return pl
    
end

"""

    Plots bonds of specified length between sites in siteList   
    SiteList::AbstractVector,
    Basis,
    pl = current();
    minDist::Real=0.,
    maxDist::Real=Basis.NNdist,
    kwargs...

"""
function plotDistBonds!(SiteList::AbstractVector,Basis,pl = current();minDist::Real=0., maxDist::Real=Basis.NNdist,label = "",kwargs...)
    for (i,Site) in enumerate(SiteList)
        Bonds = getBonds(Site,minDist,maxDist,Basis)
        filter!(x-> x in SiteList && x != Site,Bonds)
        plotBonds!(Site,Bonds,Basis,pl;color = :black,lw = 0.1,kwargs...,label = "")
    end    
    plotBond!(first(SiteList),first(SiteList),Basis;color = :black,lw = 0.1,label=label,kwargs...)# dummy to write label only once
    return pl
end

plotDistBonds!(Site::Rvec,Basis,pl = current();kwargs...) = plotBonds!(Site,getBonds(Site,minDist,maxDist,Basis),Basis,pl;kwargs...)

function plotDistBonds!(System::Geometry,Basis::Basis_Struct,sites,pl = current();kwargs...)
    println(sites)
    plotDistBonds!(sites,Basis,pl = current();kwargs...)
end

plotDistBonds(System::Geometry,Basis::Basis_Struct;kwargs...) = plotDistBonds!(generatePairSites(System.NLen,Basis),Basis,pl = plot();kwargs...)


function getVertexR(Vertex::AbstractVector,Lattice)
    @unpack Basis,PairList,PairTypes = Lattice
    @unpack refSites = Basis
    norm(i) = dist(refSites[PairTypes[i].xi],PairList[i],Basis)
    norms = norm.(eachindex(PairList))
    return norms,Vertex
end
