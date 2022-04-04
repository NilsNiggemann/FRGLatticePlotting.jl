strd(x,n=2) = string(round(x,digits = n))
strd3(x) = string(round(x,digits = 3))

# planes for Fourier space
@inline hhlplane(x,z) = SA[x,x,z]
@inline xyplane(x,y) = SA[x,y]
@inline zzerocut(x,y) = SA[x,y,0]
@inline function sphereplane(origin,radius)
    sphere(θ,ϕ) = radius*SA[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ),cos(θ)] +origin
end
hhllabels() = Dict([:xlabel => L"[hh0]",:ylabel => L"[00l]"])