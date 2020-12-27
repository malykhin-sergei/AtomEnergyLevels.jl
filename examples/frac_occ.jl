using AtomEnergyLevels
using PyPlot
pygui(true);

function cobalt()
    E = []
    for q = 7:0.1:9
        push!(E, lda(27, conf = conf_enc("[Ar] 3d$q 4s$(9 - q)")).energy.total)
    end
    return E
end

E = cobalt()

begin
    title(L"Energy of the cobalt atom with configuration $3d^q\,4s^{9-q}$")
    xlabel(L"q")
    ylabel(L"E - E_{min}")

    ax1 = PyPlot.axes()
    ax1.set_xlim([7,9])
    
    plot(7:0.1:9, E .- minimum(E))
    grid("on")
end