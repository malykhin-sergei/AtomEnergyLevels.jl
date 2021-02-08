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

#= 
Jorge Garza, Rubicelia Vargas, and Alberto Vela. 
Numerical self-consistent-field method to solve the Kohn-Sham equations
in confined many-electron atoms // Phys. Rev. E. 58 (1998) 3949-3954
DOI:10.1103/PhysRevE.58.3949
=#

r = range(0.90, 6, length = 50)
E, HOMO = main(2, r, config = c"[He]", xcfun = xc_lda(LDA_X, ρ -> 0))

begin
    title("Total energy of the helium atom for several confinements")
    xlabel("Rc (a.u.)")
    ylabel("Total energy (a.u.)")
    grid("on")

    plot(r, E)

    # TABLE III. 
    plot(1:6, [1.354, -2.384, -2.682, -2.718, -2.723, -2.724], linestyle="None", marker="o", 
         label="Jorge Garza et al // Phys. Rev. E. 58 (1998) 3949")
    legend(loc = "upper right", fancybox = "true")
end

#=
Jan C. A. Boeyens. Ionization Radii of Compressed Atoms. 
J. CHEM. SOC. FARADAY TRANS., 1994, 90(22), 3377-3381

NB. Hartree-Fock-Slater method is simply Xα functional, where α = 1

Data taken from Fig. 4 One-electron compression curve for the valence level of oxygen.
=#

r = [2.2059, 2.2752, 2.3446, 2.4002, 2.4696, 2.5526, 
     2.6494, 2.7598, 2.8839, 3.0218, 3.1733, 3.3522, 
     3.5722, 3.8195, 4.0805, 4.3688, 4.6980, 5.0272, 
     5.3564, 5.6992, 6.0283, 6.3711, 6.7139, 7.0429, 
     7.3582, 7.7010, 8.0437, 8.3728, 8.7019, 9.0446, 
     9.3874]

ϵₕ = [0.95604,  0.81394,  0.66592,  0.51786,  0.37576, 
      0.23370,  0.09166, -0.03848, -0.17452, -0.31052, 
     -0.44649, -0.57054, -0.68264, -0.78281, -0.87110, 
     -0.94154, -0.97038, -0.99330, -1.01622, -1.02725, 
     -1.03239, -1.03749, -1.04259, -1.03587, -1.04104, 
     -1.04021, -1.03939, -1.04452, -1.04373, -1.04291, 
     -1.04208] # Ry

E, HOMO = main(8, r, config = c"[He] 2s2 2p4", xcfun = xc_lda(ρ -> LDA_X(ρ, α = 1), ρ -> 0))

begin
    ax1 = PyPlot.axes()
    ax1.set_xlim([ 2.0, 10])
    ax1.set_ylim([-1.5,  1])
    plot(r, ϵₕ)
    plot(r, 2HOMO)
    grid("on")
end
