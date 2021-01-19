using AtomEnergyLevels
using PyPlot
pygui(true)

b = 1/2 * (3π/4)^(2/3)

x = -35.0:0.1:5.0; r = exp.(x);

Z = 92; ρ = lda(Z, x).density;

begin
    title("Radial electronic density of the U and Thomas-Fermi atoms")
    xlabel(L"\frac{r Z^{1/3}}{b},\,a.u.")
    ylabel(L"\rho Z^{-\frac{4}{3}}")

    grid("on")

    ax1 = PyPlot.axes()
    ax1.set_xlim([0.0001,100])
    PyPlot.xscale("log")
    
    plot(r/b, 4π * r .^ 2 .* TF.(r, 1), label = "TF")
    plot(r/b * Z^(1/3), 4π * r .* ρ / Z^(4/3), label = "U (LDA)")
    legend(loc = "upper right", fancybox = "true")
end
