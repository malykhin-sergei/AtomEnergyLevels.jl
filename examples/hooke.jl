"""
# Hooke's atom

Find the total energy and density of the 
[Hooke's atom](https://en.wikipedia.org/wiki/Hooke's_atom) using LDA DFT.
"""

using AtomEnergyLevels
using SpecialFunctions, PyPlot
pygui(true);

x = -30.0:0.1:20.0;
r  = exp.(x);
N² = 1 / (π^(3/2) * (8 + 5 * sqrt(π)));
ρₑ = @. 2N² * exp(-1/2 * r^2) * (sqrt(π/2) * (7/4 + 1/4 * r^2 + (r + 1/r) * erf(r/sqrt(2))) + exp(-1/2 * r^2));

α = 0.798
xalpha = lda(2, x, conf = 2, Vex = r -> 1/8 * r^2,
             xc! = (ρ, vxc, exc) -> Xα!(ρ, vxc, exc, α = α), β = 0.7);

svwn5 = lda(2, x, conf = 2, Vex = r -> 1/8 * r^2,  β = 0.7);

begin
    title("Hooke's atom density")
    xlabel(L"\frac{r Z^{1/3}}{b},\,a.u.")
    ylabel(L"\rho Z^{-\frac{4}{3}}")
    
    ax1 = PyPlot.axes()
    ax1.set_xlim([0,5])

    plot(r, ρₑ, label = "Exact")
    plot(r, xalpha.density ./ r, label = "Xα, α=$α")
    plot(r, svwn5.density ./ r, label = "LDA")

    legend(loc = "upper right", fancybox = "true")
    grid("on")
end
