
getPoint(R::Rvec,Basis) = Point(getCartesian(R,Basis)) 
Makie.scatter!(ax,Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = scatter!(ax,getPoint.(Rs,(Basis,)),args...;kwargs...)
Makie.scatter!(Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = scatter!(current_axis(),getPoint.(Rs,(Basis,)),args...;kwargs...)

Makie.lines!(ax,Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = lines!(ax,getPoint.(Rs,(Basis,)),args...;kwargs...)
Makie.lines!(Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = lines!(current_axis(),getPoint.(Rs,(Basis,)),args...;kwargs...)

Makie.scatter(Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = scatter(getPoint.(Rs,(Basis,)),args...;kwargs...)
Makie.lines(Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = lines(getPoint.(Rs,(Basis,)),args...;kwargs...)

function pairsPlot!(ax,PairList,Basis,args...;colors = (:blue,:red,:black,:cyan,:yellow,:green,:pink,:orange,:lime,:brown,:grey),colorBasis = false,kwargs...)
    uniquepairs = unique(PairList)
    
    colors = colors[colors[Basis.SiteType[R.b]] for R in uniquepairs] 

    scatter!(ax,uniquepairs,color = colors,args...;kwargs...)
end

function getStandardFigure(::Rvec_2D)
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1],aspect = 1,xticks = SimpleTicks(),yticks = SimpleTicks(),xlabel = L"x", ylabel = L"y")
    fig,ax
end

function getStandardFigure(::Rvec_3D)
    fig = Figure(resolution = (800, 600))
    ax = Axis3(fig[1, 1],aspect = (1,1,1),xticks = SimpleTicks(),yticks = SimpleTicks(),zticks = SimpleTicks(),xlabel = L"x", ylabel = L"y",zlabel = L"z")
    fig,ax
end

function pairsPlot(PairList::AbstractVector{R},Basis,args...;kwargs...) where R <: Rvec
    fig,ax = getStandardFigure(R)
    plotSites!(ax,sites,Basis)
    fig
end

function connectWithinDist(points::AbstractVector{P},dmin,dmax;tol = 1e-8) where {P <: Point}

    missingPoint = Point((NaN for _ in first(points))...)
    isMissingPoint(p) = any(isnan.(p))

    connectedPoints = P[]
    for p1 in points
        for p2 in points
            if dmin-tol < abs(norm(p1-p2) ) < dmax+tol
                push!(connectedPoints,p1)
                push!(connectedPoints,p2)
            elseif isempty(connectedPoints) || !isMissingPoint(last(connectedPoints))
                push!(connectedPoints,missingPoint)
            end
        end
    end
    return connectedPoints
end

connectWithinDist(points::AbstractVector{Rvec},Basis,dmin,dmax) = connectWithinDist(getPoint.(points,Basis),dmin,dmax)
"""

    Plots bonds of specified length between sites in siteList   
    ax,
    SiteList::AbstractVector{<:Rvec},
    Basis,
    minDist::Real=0.,
    maxDist::Real=Basis.NNdist,
    kwargs...

"""
function plotDistBonds!(ax,SiteList::AbstractVector{<:Rvec},Basis,args...;minDist::Real=0., maxDist::Real=Basis.NNdist,kwargs...)
    connectedPoints = connectWithinDist(SiteList,Basis,minDist,maxDist)
    lines!(connectedPoints,color = (:black,1),linestyle = :solid;kwargs...)
    return ax
end

function plotDistBonds(SiteList::AbstractVector{<:Rvec},Basis,args...;kwargs...)
    fig,ax = getStandardFigure(R)
    plotDistBonds!(ax,SiteList,Basis;minDist=minDist,maxDist=maxDist,args...,kwargs...)
    fig
end

"""Plot all sites and inequivalent pairs"""
function plotSystem(System,Basis;
    plotAll = true,
    refSite = nothing,
    markersize = 5,
    inequivColor = :green,
    inequivalpha = 0.5,
    plotBonds=true,
    plotCouplings=true,
    CouplingColors = nothing,
    Bonds = [(minDist = Basis.NNdist-1e-3,maxDist = Basis.NNdist+1e-3,colorRGB = [0,0,0])],
    bondlw = 7,
    allpairs = unique!(SpinFRGLattices.sortedPairList(System.NLen,Basis)[1]),
    args...,
    kwargs...)
    (;PairList,OnsitePairs )= System
    
    indices = copy(OnsitePairs)
    push!(indices,length(PairList)) # get final index
    if refSite === nothing 
        plotpairs = unique(PairList)
    else
        # allpairs = unique!(generatePairSites(System.NLen,Basis,Basis.refSites[refSite]))
        plotpairs = PairList[indices[refSite]:indices[refSite+1]-1]
    end
    filter!(x-> x in allpairs,plotpairs)

    plotAll || (allpairs = plotpairs)
    fig,ax = getStandardFigure(first(plotpairs))

    if plotBonds
        if bondlw isa Real
            bondlw = [bondlw for b in Bonds]
        end
        for (lw,b) in zip(bondlw,Bonds)
            plotDistBonds!(ax,allpairs,Basis,minDist = b.minDist, maxDist = b.maxDist;lw,color = Plots.Colors.RGB((b.colorRGB./255)...))
        end
    end
    refSite !== nothing && pairsPlot!(ax,[Basis.refSites[refSite]], Basis,color = "black",markersize = 1.5*markersize,markershape = '×')

    plotAll && pairsPlot!(ax,plotpairs,Basis,color = inequivColor,alpha = inequivalpha,markersize = 2*markersize)
    refSite !== nothing && pairsPlot!(ax,[Basis.refSites[refSite]],Basis,color = "darkred",markershape = '×',markersize = 1.5*markersize)

    plotCouplings && plotCouplings!(System,Basis;refSite = refSite,colors = CouplingColors)

    pairsPlot!(ax,allpairs,Basis,markersize = markersize;kwargs...)

    return pl
end