using FRGLatticePlotting, SpinFRGLattices
S = Pyrochlore.getPyrochlore(3, [5.,-3,0,1])
fig = plotSystem(S,Pyrochlore.Basis,refSite = 1,colorBasis=false,inspect = true)
##
S = SquareKagome.getMirrorSquareKagome(6, 5.,-3)
Basis = SquareKagome.ShurikenBasis() #chose any you want for visualization
fig = plotSystem(S,Basis,
    refSite = 2,
    inspect = true,
    inequivScale = 2.5,
    Bonds = [
        (minDist = Basis.NNdist,maxDist = Basis.NNdist,colorRGB = [0,0,0]),
        (minDist = dist(Rvec(0,0,2),Rvec(0,0,4),Basis),maxDist = dist(Rvec(0,0,2),Rvec(0,0,4),Basis),colorRGB = [200,50,50]),
        (minDist = dist(Rvec(0,0,2),Rvec(-1,1,4),Basis),maxDist = dist(Rvec(0,0,2),Rvec(-1,1,4),Basis),colorRGB = [10,70,70])
    ],
    bondlw = [4,3,2],
    bondstyle = [:solid,:dot,:dash]
)
