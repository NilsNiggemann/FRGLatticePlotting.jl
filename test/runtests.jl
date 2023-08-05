using FRGLatticePlotting, SpinFRGLattices
S = Pyrochlore.getPyrochlore(3, [5.,-3,0,1])
fig = plotSystem(S,Pyrochlore.Basis,refSite = 1,colorBasis=false,inspect = true)
# x = DataInspector(fig)
# fig
##
S = SquareKagome.getMirrorSquareKagome(3, 5.,-3)
fig = plotSystem(S,SquareKagome.Basis,refSite = 1,inspect = true)