using FRGLatticePlotting, SpinFRGLattices

S = Pyrochlore.getPyrochlore(4, [1.,0.7,.5,0.1])
plotSystem(S,Pyrochlore.Basis)
