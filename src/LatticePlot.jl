export pairsPlot, plotSystem, plotCouplings!,plotCorrelations!

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
function plotSystem(System,Basis;plotAll = true,refSite = 0,markersize = 5,inequivColor = "green",inequivalpha = 0.5,kwargs...)
    @unpack PairList,OnsitePairs = System
    indices = copy(OnsitePairs)
    push!(indices,length(PairList)) # get final index
    allpairs = unique!(SpinFRGLattices.sortedPairList(System.NLen,Basis)[1])
    if refSite == 0
        plotpairs = unique(PairList)
    else
        # allpairs = unique!(generatePairSites(System.NLen,Basis,Basis.refSites[refSite]))
        plotpairs = PairList[indices[refSite]:indices[refSite+1]]
    end

    plotAll || (allpairs = plotpairs)
    pl = pairsPlot(allpairs,Basis,markersize = markersize;kwargs...)
    plotAll && pairsPlot(plotpairs,Basis,pl,color = inequivColor,alpha = inequivalpha,markersize = 2*markersize)
    return pl
end

function plotCorrelations!(System,Basis,couplings,pl=current();kwargs...)
    prepdata(r1,r2) = Tuple(SA[r1[i],r2[i]] for i in eachindex(r1))

    @unpack PairTypes,PairList = System
    for (type,R_pair,J) in zip(PairTypes,PairList,couplings)
        x = type.xi
        r_ref = getCartesian(Basis.refSites[x],Basis)
        if abs(J)> 1E-10
            r_pair = getCartesian(R_pair,Basis)
            points = prepdata(r_ref,r_pair)
            plot!(pl,points...,label = round(J,digits=3),lw = 1+4*abs(J);kwargs...)
        end
    end
    return pl
end

function plotCouplings!(System,Basis,pl=current();kwargs...)
    plotCorrelations(System,Basis,System.couplings,pl,kwargs...)
end
