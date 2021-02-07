using AtomEnergyLevels
using PyPlot
pygui(true)

function main(z, rgrid; config = atom[z], xcfun = xc_lda(LDA_X, LDA_C_VWN), n = 500)
    x, en, HOMO = log.(rgrid), similar(rgrid), similar(rgrid)
    for (i, xmax) in enumerate(x)        
        results = lda(z, range(-35, xmax, length = n), conf = config, xc_func! = xcfun)
        en[i] = results.energy.total
        HOMO[i] = maximum(map(λ -> λ[1], values(results.orbitals)))
    end
    return en, HOMO
end

begin
    r = range(0.90, 6, length = 50)
    E, HOMO = main(2, r, config = c"[He]", xcfun = xc_lda(LDA_X, ρ -> 0))

    title("Total energy of the helium atom for several confinements")
    xlabel("Rc (a.u.)")
    ylabel("Total energy (a.u.)")
    grid("on")

    plot(r, E)

    # TABLE III. Jorge Garza, Rubicelia Vargas, and Alberto Vela. 
    # Numerical self-consistent-field method to solve the Kohn-Sham equations
    # in confined many-electron atoms // Phys. Rev. E. 58 (1998) 3949-3954
    # DOI:10.1103/PhysRevE.58.3949

    plot(1:6, [1.354, -2.384, -2.682, -2.718, -2.723, -2.724], linestyle="None", marker="o")
end
