
@setup_workload begin

    @compile_workload begin
        if HASGSL
            GLMakie.activate!()
            GLMakie.Makie.inline!(false)
        end
        using FRGLatticePlotting, SpinFRGLattices


        S1 = Pyrochlore.getPyrochlore(3, [5.,-3,0,1])
        S2 = SquareKagome.getMirrorSquareKagome(6, 5.,-3)
        
        for (Basis,S) in zip(
            (Pyrochlore.Basis,SquareKagome.ShurikenBasis()),
            (S1,S2)
            )
            fig = plotSystem(S,Basis,refSite = 1,inspect = true,inequivScale = 2.5)
            
            FRGLatticePlotting.getInspector(S.PairList,Basis)
            FRGLatticePlotting.getPairNumberInspector(S.PairList,Basis)
            FRGLatticePlotting.getCorrelationInspector(1.0,first(S.PairList))
        end

    end
end

