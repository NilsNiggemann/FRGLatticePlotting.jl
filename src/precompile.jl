
@setup_workload begin

    @compile_workload begin
        GLMakie.activate!()
        GLMakie.Makie.inline!(false)
        using FRGLatticePlotting, SpinFRGLattices
        S = Pyrochlore.getPyrochlore(3, [5.,-3,0,1])
        fig = plotSystem(S,Pyrochlore.Basis,refSite = 1,colorBasis=false,inspect = true)
        
        S = SquareKagome.getMirrorSquareKagome(6, 5.,-3)
        fig = plotSystem(S,SquareKagome.ShurikenBasis(),refSite = 2,inspect = true,inequivScale = 2.5)
    end
end

