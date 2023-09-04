getPoint(R::Rvec, Basis) = Point(getCartesian(R, Basis))
scatterRvec!(ax::Makie.Block, Rs::AbstractVector{<:Rvec}, Basis, args...; kwargs...) = scatter!(ax::Makie.Block, getPoint.(Rs, Ref(Basis)), args...; kwargs...)
scatterRvec!(Rs::AbstractVector{<:Rvec}, Basis, args...; kwargs...) = scatter!(current_axis(), getPoint.(Rs, Ref(Basis)), args...; kwargs...)

linesRvec!(ax::Makie.Block, Rs::AbstractVector{<:Rvec}, Basis, args...; kwargs...) = lines!(ax, getPoint.(Rs, Ref(Basis)), args...; kwargs...)
linesRvec!(Rs::AbstractVector{<:Rvec}, Basis, args...; kwargs...) = lines!(current_axis(), getPoint.(Rs, Ref(Basis)), args...; kwargs...)

scatterRvec(Rs::AbstractVector{<:Rvec}, Basis, args...; kwargs...) = scatter(getPoint.(Rs, Ref(Basis)), args...; kwargs...)
lineRvecs(Rs::AbstractVector{<:Rvec}, Basis, args...; kwargs...) = lines(getPoint.(Rs, Ref(Basis)), args...; kwargs...)

function pairsPlot!(ax, PairList, Basis, args...; colors=[:black, :red, :blue, :grey, :darkred, :pink, :orange, :lime, :cyan, :yellow, :green, :brown], colorBasis=false, kwargs...)
    uniquepairs = unique(PairList)

    inds =
        if colorBasis
            [R.b for R in uniquepairs]
        else
            [Basis.SiteType[R.b] for R in uniquepairs]
        end
    col = colors[inds]
    scatterRvec!(ax, uniquepairs, Basis, color=col, args...; inspector_label=getInspector(PairList, Basis), kwargs...)
end

function getStandardFigure(::Type{Rvec_2D};inspect=true,kwargs...)
    fig = Figure(resolution=(800, 600))
    ax = Axis(fig[1, 1], aspect=1, xticks=SimpleTicks(), yticks=SimpleTicks(), xlabel=L"x", ylabel=L"y";kwargs...)
    if inspect
        DataInspector(fig)
    end
    fig, ax
end


function getStandardFigure(::Type{Rvec_3D};inspect=true,kwargs...)
    fig = Figure(resolution=(800, 600))
    ax = Axis3(fig[1, 1], aspect=(1, 1, 1), xticks=SimpleTicks(), yticks=SimpleTicks(), zticks=SimpleTicks(), xlabel=L"x", ylabel=L"y", zlabel=L"z";kwargs...)
    if inspect
        DataInspector(fig)
    end
    fig, ax
end

getStandardFigure(R::Rvec;kwargs...) = getStandardFigure(typeof(R); kwargs...)

function pairsPlot(PairList::AbstractVector{R}, Basis, args...; kwargs...) where {R<:Rvec}
    fig, ax = getStandardFigure(R)
    scatterRvec!(ax, PairList, Basis) 
    fig
end

function connectWithinDist(points::AbstractVector{P}, dmin, dmax; tol=1e-8) where {P<:Point}

    missingPoint = Point((NaN for _ in first(points))...)
    isMissingPoint(p) = any(isnan.(p))

    connectedPoints = P[]
    for p1 in points
        for p2 in points
            if dmin - tol < abs(norm(p1 - p2)) < dmax + tol
                push!(connectedPoints, p1)
                push!(connectedPoints, p2)
            elseif isempty(connectedPoints) || !isMissingPoint(last(connectedPoints))
                push!(connectedPoints, missingPoint)
            end
        end
    end
    return connectedPoints
end

function connectWithinDist(p1::P, points::AbstractVector{P}, dmin, dmax; tol=1e-8) where {P<:Point}

    missingPoint = Point((NaN for _ in first(points))...)
    isMissingPoint(p) = any(isnan.(p))

    connectedPoints = P[]
    for p2 in points
        if dmin - tol < abs(norm(p1 - p2)) < dmax + tol
            push!(connectedPoints, p1)
            push!(connectedPoints, p2)
        elseif isempty(connectedPoints) || !isMissingPoint(last(connectedPoints))
            push!(connectedPoints, missingPoint)
        end
    end
    return connectedPoints
end

connectWithinDist(points::AbstractVector{<:Rvec}, Basis, dmin, dmax) = connectWithinDist(getPoint.(points, Ref(Basis)), dmin, dmax)
connectWithinDist(p1::Rvec, points::AbstractVector{<:Rvec}, Basis, dmin, dmax) = connectWithinDist(getPoint(p1, Basis), getPoint.(points, Ref(Basis)), dmin, dmax)


"""

    Plots bonds of specified length between sites in siteList   
    ax,
    SiteList::AbstractVector{<:Rvec},
    Basis,
    minDist::Real=0.,
    maxDist::Real=Basis.NNdist,
    kwargs...

"""
function plotDistBonds!(ax, SiteList::AbstractVector{<:Rvec}, Basis, args...; minDist::Real=0.0, maxDist::Real=Basis.NNdist, kwargs...)
    connectedPoints = connectWithinDist(SiteList, Basis, minDist, maxDist)
    lines!(connectedPoints, color=(:black, 1), linestyle=:solid; kwargs...)
    return ax
end

plotDistBonds!(SiteList::AbstractVector{<:Rvec}, Basis, args...; kwargs...) = plotDistBonds!(current_axis(), SiteList, Basis, args...; kwargs...)

function plotDistBonds(SiteList::AbstractVector{R}, Basis, args...; kwargs...) where {R<:Rvec}
    fig, ax = getStandardFigure(R)
    plotDistBonds!(ax, SiteList, Basis; minDist=minDist, maxDist=maxDist, args..., kwargs...)
    fig
end

strd(x) = string(round(x, digits=3))
getstring(R::Rvec_3D) = string("(", R.n1, ", ", R.n2, ", ", R.n3, ", ", R.b, ")")
getstring(R::Rvec_2D) = string("(", R.n1, ", ", R.n2, ", ", R.b, ")")
getstring(r::AbstractArray) = string("[", join(strd.(r),", "),"]")

function getInspector(RvecList::AbstractVector{RT}, Basis) where {RT<:Rvec}
    fig, ax = getStandardFigure(RT)
    function inspector(self, i, p)
        R = RvecList[i]
        r = getPoint(R, Basis)
        r = getstring(r)
        latexstring(getstring(R), " →\n", r)
    end
    return inspector
end


function getPairNumberInspector(PairList::AbstractVector{<:Rvec}, Basis,offset=0)
    function inspector(self, i, p)
        R = PairList[i]
        latexstring(getstring(R), " → ", i+offset)
    end
    return inspector
end

function getCorrelationInspector(J, R2)
    R2str = getstring(R2)
    function inspector(self, i, p)
        latexstring("J_{ij}$(R2str)", " = ", J)
    end
    return inspector
end

"""
Plot all sites and inequivalent pairs

    function plotSystem!(ax, System,Basis,args...;
        plotAll=true,
        refSite=1,
        markersize=20,
        inequivColor=:green,
        inequivalpha=0.6,
        plotBonds=true,
        plotCouplings=true,
        colorBasis=false,
        CouplingColors=nothing,
        bondDist = Basis.NNdist,
        Bonds=[(minDist=bondDist- 1e-3, maxDist=bondDist + 1e-3, colorRGB=[0, 0, 0])],
        bondlw=4,
        bondstyle = [:solid for _ in Bonds],
        inequivScale=2.5,
        allpairs=unique!(SpinFRGLattices.sortedPairList(System.NLen, Basis)[1]),
        linewidthscaling=J -> 2 * abs(J),
        kwargs...
    )
"""
function plotSystem!(ax, System, Basis, args...;
    plotAll=true,
    refSite=1,
    markersize=20,
    inequivColor=:green,
    inequivalpha=0.6,
    plotBonds=true,
    plotCouplings=true,
    colorBasis=false,
    CouplingColors=nothing,
    bondDist = Basis.NNdist,
    Bonds=[(minDist=bondDist- 1e-3, maxDist=bondDist + 1e-3, colorRGB=[0, 0, 0])],
    bondlw=4,
    
    bondstyle = [:solid for _ in Bonds],
    inequivScale=2.5,
    allpairs=unique!(SpinFRGLattices.sortedPairList(System.NLen, Basis)[1]),
    linewidthscaling=J -> 2 * abs(J),
    kwargs...
)

    (; PairList, OnsitePairs) = System

    if plotBonds
        if bondlw isa Real
            bondlw = [bondlw for b in Bonds]
        end
        for (lw, b,ls) in zip(bondlw, Bonds, bondstyle)
            plotDistBonds!(ax, allpairs, Basis, minDist=b.minDist, maxDist=b.maxDist; linestyle = ls, linewidth = lw, color=Makie.Colors.RGB((b.colorRGB ./ 255)...), inspectable=false)
        end
    end

    plotCouplings && plotCorrelations!(ax, System, Basis, System.couplings; allpairs, refSite, colors=CouplingColors, linewidthscaling)

    indices = copy(OnsitePairs)
    push!(indices, length(PairList)) # get final index
    plotpairs = if refSite === nothing
        unique(PairList)
    else
        PairList[indices[refSite]:indices[refSite+1]]
    end
    filter!(x -> x in allpairs, plotpairs)

    plotAll || (allpairs = plotpairs)


    plotAll && pairsPlot!(ax, plotpairs, Basis, color=inequivColor, alpha=inequivalpha, markersize=inequivScale * markersize; inspector_label=getPairNumberInspector(plotpairs, Basis,indices[refSite]-1))

    R0 = Basis.refSites[refSite]
    scatterRvec!(ax, [Basis.refSites[refSite]], Basis, color=:darkred, marker='×', markersize=5.0 * markersize, inspector_label=(self, i, p) -> latexstring("R_0 = ", round.(getPoint(R0, Basis), digits=3)))

    pairsPlot!(ax, allpairs, Basis, markersize=markersize; colorBasis, inspectable=true, kwargs...)

end

"""
Plot all sites and inequivalent pairs

    plotSystem(ax, System,Basis,args...;
        plotAll=true,
        refSite=1,
        markersize=20,
        inequivColor=:green,
        inequivalpha=0.6,
        plotBonds=true,
        plotCouplings=true,
        colorBasis=false,
        CouplingColors=nothing,
        Bonds=[(minDist=Basis.NNdist - 1e-3, maxDist=Basis.NNdist + 1e-3, colorRGB=[0, 0, 0])],
        bondlw=4,
        bondstyle = [:solid for _ in Bonds],
        inequivScale=2.5,
        allpairs=unique!(SpinFRGLattices.sortedPairList(System.NLen, Basis)[1]),
        inspect=true,
        linewidthscaling=J -> 2 * abs(J),
        kwargs...
    )

"""
function plotSystem(System, Basis, args...;
    Axis = (;),
    inspect=true,
    kwargs...
)
    fig, ax = getStandardFigure(first(System.PairList);inspect,Axis...)

    plotSystem!(ax, System, Basis, args...; kwargs...)

    return fig
end

function plotCorrelations!(ax, System, Basis, Correlations::AbstractVector, args...; allpairs=unique!(SpinFRGLattices.sortedPairList(System.NLen, Basis)[1]), minCorr=1e-14, refSite=1, colors=nothing,
    linewidthscaling=J -> 3 * abs(J),
    kwargs...)

    (; PairTypes, PairList) = System
    inds1 = findall(x -> x.xi == refSite, PairTypes)

    inds = findall(x -> abs(x) > minCorr, Correlations)
    inds = intersect(inds, inds1)

    Correlations = Correlations[inds]
    PairList = PairList[inds]
    label(J) = latexstring(round(J, digits=3))
    R0 = Basis.refSites[refSite]
    for i in eachindex(Correlations, PairList)
        J = Correlations[i]
        color = J > 0 ? :red : :blue
        R = PairList[i]
        linesRvec!(ax, [R0, R], Basis; color, label=label(J), linewidth=linewidthscaling(J), inspector_label=getCorrelationInspector(J, R), args..., kwargs...)
    end
    
    if length(Correlations) > 0
        axislegend(L"J_{ij}")
    end

    return ax

end
