
getPoint(R::Rvec,Basis) = Point(getCartesian(R,Basis)) 
scatterRvec!(ax::Makie.Block,Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = scatter!(ax::Makie.Block,getPoint.(Rs,Ref(Basis)),args...;kwargs...)
scatterRvec!(Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = scatter!(current_axis(),getPoint.(Rs,Ref(Basis)),args...;kwargs...)

linesRvec!(ax::Makie.Block,Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = lines!(ax,getPoint.(Rs,Ref(Basis)),args...;kwargs...)
linesRvec!(Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = lines!(current_axis(),getPoint.(Rs,Ref(Basis)),args...;kwargs...)

scatterRvec(Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = scatter(getPoint.(Rs,Ref(Basis)),args...;kwargs...)
lineRvecs(Rs::AbstractVector{<:Rvec},Basis,args...;kwargs...) = lines(getPoint.(Rs,Ref(Basis)),args...;kwargs...)

function pairsPlot!(ax,PairList,Basis,args...;colors = [:blue,:red,:black,:cyan,:yellow,:green,:pink,:orange,:lime,:brown,:grey],colorBasis = false,kwargs...)
    uniquepairs = unique(PairList)
    
    inds = [Basis.SiteType[R.b] for R in uniquepairs]
    col = colors[inds]

    scatterRvec!(ax,uniquepairs,Basis,color = col,args...;kwargs...)
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

function connectWithinDist(p1::P,points::AbstractVector{P},dmin,dmax;tol = 1e-8) where {P <: Point}

    missingPoint = Point((NaN for _ in first(points))...)
    isMissingPoint(p) = any(isnan.(p))

    connectedPoints = P[]
    for p2 in points
        if dmin-tol < abs(norm(p1-p2) ) < dmax+tol
            push!(connectedPoints,p1)
            push!(connectedPoints,p2)
        elseif isempty(connectedPoints) || !isMissingPoint(last(connectedPoints))
            push!(connectedPoints,missingPoint)
        end
    end
    return connectedPoints
end

connectWithinDist(points::AbstractVector{<:Rvec},Basis,dmin,dmax) = connectWithinDist(getPoint.(points,Ref(Basis)),dmin,dmax)
connectWithinDist(p1::Rvec,points::AbstractVector{<:Rvec},Basis,dmin,dmax) = connectWithinDist(getPoint(p1,Basis),getPoint.(points,Ref(Basis)),dmin,dmax)


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
function plotSystem(System,Basis,args...;
    plotAll = true,
    refSite = nothing,
    markersize = 15,
    inequivColor = :green,
    inequivalpha = 0.7,
    plotBonds=true,
    plotCouplings=true,
    CouplingColors = nothing,
    Bonds = [(minDist = Basis.NNdist-1e-3,maxDist = Basis.NNdist+1e-3,colorRGB = [0,0,0])],
    bondlw = 7,
    allpairs = unique!(SpinFRGLattices.sortedPairList(System.NLen,Basis)[1]),
    kwargs...)
    (;PairList,OnsitePairs )= System
    
    indices = copy(OnsitePairs)
    push!(indices,length(PairList)) # get final index
    plotpairs = if refSite === nothing 
        unique(PairList)
    else
        PairList[indices[refSite]:indices[refSite+1]-1]
    end
    filter!(x-> x in allpairs,plotpairs)

    plotAll || (allpairs = plotpairs)
    fig,ax = getStandardFigure(first(plotpairs))

    if plotBonds
        if bondlw isa Real
            bondlw = [bondlw for b in Bonds]
        end
        for (lw,b) in zip(bondlw,Bonds)
            plotDistBonds!(ax,allpairs,Basis,minDist = b.minDist, maxDist = b.maxDist;lw,color = Makie.Colors.RGB((b.colorRGB./255)...))
        end
    end

    
    plotAll && pairsPlot!(ax,plotpairs,Basis,color = inequivColor,alpha = inequivalpha,markersize = 2*markersize)

    if refSite !== nothing
        plotCouplings && plotCorrelations!(ax,System,Basis,System.couplings;allpairs,refSite = Basis.refSites[refSite],colors = CouplingColors)
    end
    
    pairsPlot!(ax,allpairs,Basis,markersize = markersize;kwargs...)
    refSite !== nothing && scatterRvec!(ax,[Basis.refSites[refSite]],Basis,color = :darkred,marker = '×',markersize = 3.5*markersize)

    # refSite !== nothing && scatterRvec!(ax,[Basis.refSites[refSite]], Basis)

    return fig
end

function plotCorrelations!(ax,System,Basis,Correlations::AbstractVector,args...; allpairs = unique!(SpinFRGLattices.sortedPairList(System.NLen,Basis)[1]),minCorr = 1e-14,refSite = nothing,colors = nothing,kwargs...)

    (;PairTypes,PairList) = System

    inds = findall(x-> abs(x)>minCorr,Correlations)
    Correlations = Correlations[inds]
    PairList = PairList[inds]
    label(J) = round(J,digits=3)
    lw(J) = 0.1+3*abs(J)

    for i in eachindex(Correlations,PairList)
        J = Correlations[i]
        color = J > 0 ? :red : :blue
        linesRvec!(ax,[refSite,PairList[i]],Basis;color,linewidth = lw(J),args...,kwargs...)
    end

    return ax
    
end