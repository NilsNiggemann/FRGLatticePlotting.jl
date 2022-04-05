strd(x,n=2) = string(round(x,digits = n))
strd3(x) = string(round(x,digits = 3))

allOccurIn(name,args...) = all((occursin(arg,name) for arg in args))
findNames(names,args...) = findall(x->allOccurIn(x,args...),names)
NameFilter(names,args...) = filter(x->allOccurIn(x,args...),names)
OnlyIndex(names,args...) = only(findNames(names,args...))

function getNumberFromName(Name,subName)
    res_string = split(Name,subName*"=")[end]
    for i in length(res_string):-1:1
        N = tryparse(Int,res_string[1:i])
        if N !== nothing
            return N
        end
    end
    error("Could not get ", subName, "from string ",Name)
end


# planes for Fourier space
@inline hhlplane(x,z) = SA[x,x,z]
@inline xyplane(x,y) = SA[x,y]
@inline zzerocut(x,y) = SA[x,y,0]
@inline function sphereplane(origin,radius)
    sphere(θ,ϕ) = radius*SA[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ),cos(θ)] +origin
end
hhllabels() = Dict([:xlabel => L"[hh0]",:ylabel => L"[00l]"])


getPair(R1::Rvec,R2::Rvec,Lattice) = SpinFRGLattices.pairNumber(R1,R2,Lattice.System.PairList,Lattice.System.PairTypes,Lattice.Basis,Lattice.pairToInequiv)

function getPairs(R1,R2,UnitCellTranslation::StaticVector,Lattice::LatticeInfo)
    @unpack System,Basis,pairToInequiv =Lattice
    Pair1 = getPair(R1,R2,Lattice)
    R1_pr = translation(R1,UnitCellTranslation,Basis)
    R2_pr = translation(R2,UnitCellTranslation,Basis)
    Pair2 = getPair(R1_pr,R2_pr,Lattice)
    return Pair1,Pair2
end

function DimerResponse(Chi1,Chi2,delta)
    return @. (Chi1 - Chi2) / (Chi1 + Chi2) / delta
end
function DimerResponse(R1::Rvec,R2::Rvec,UnitCellTranslation::StaticVector,Chi_LR,delta,Lattice::LatticeInfo)
    Pair1,Pair2 = getPairs(R1,R2,UnitCellTranslation,Lattice)
    @views Chi1 = Chi_LR[:,Pair1]
    @views Chi2 = Chi_LR[:,Pair2]
    DimerResponse(Chi1,Chi2,delta)
end

function DimerResponse(R1::Rvec,R2::Rvec,Chi_LR::AbstractMatrix,Chi_LRNew::AbstractMatrix,delta,Lattice::LatticeInfo)
    Chiold = DimerFlow(R1,R2,Chi_LR,Lattice)
    ChiNew = DimerFlow(R1,R2,Chi_LRNew,Lattice)
    DimerResponse(ChiNew,Chiold,delta)
end

function DimerFlow(R1::Rvec,R2::Rvec,Chi_LR::AbstractMatrix,Lattice::LatticeInfo)
    Pair = getPair(R1,R2,Lattice)
    @views Chi1 = Chi_LR[:,Pair]
end

function stringLatex(args...)
    res = ""
    for arg in args
        if typeof(arg) == LaTeXString
            res = string(res,"\\ ",arg.s[2:end-1])
        elseif typeof(arg) == String
            res = string(res,"\\ ","\\textrm{$arg}")
        end
    end
    return string("\$",res,"\$")
end
