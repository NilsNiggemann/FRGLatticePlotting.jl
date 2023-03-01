using PMFRG, FRGLatticePlotting, SpinFRGLattices

S = SquareKagome.getSquareKagome(4, (1.,0.7,.5,0.1))
plotSystem(S,SquareKagome.Basis)
