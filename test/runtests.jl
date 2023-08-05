using FRGLatticePlotting, SpinFRGLattices

S = Pyrochlore.getPyrochlore(3, [5.,-3,0,1]./2)
plotSystem(S,Pyrochlore.Basis,refSite = 1)
