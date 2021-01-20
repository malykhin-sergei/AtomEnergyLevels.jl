using AtomEnergyLevels, Printf

const eV = 27.211384

function slater_ts(Z)
    IPs = []
    for (l, subshell) in enumerate(atom[Z])
        for nᵣ in eachindex(subshell)
            config_ts = deepcopy(atom[Z])
            config_ts[l][nᵣ] -= 0.5
            results = lda(Z, conf = config_ts).orbitals
            push!(IPs, (results[(nᵣ = nᵣ - 1, l = l - 1)].ϵᵢ, l - 1, nᵣ - 1))
        end
    end
    @printf("\n== Electron binding energies of the %s atom, eV ==\n\n", keys(atom)[Z])
    for level in sort(IPs, by = x -> first(x))
        ϵᵢ, l, nᵣ = level
        n = nᵣ + l + 1
        @printf("\t%i%s\t%12.2f\n", n, shells[l + 1], ϵᵢ*eV)
    end
end

