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

slater_ts(36)
#=
== Electron binding energies of the Kr atom, eV ==

        1S         -14097.00
        2S          -1838.84
        2P          -1676.35
        3S           -266.07
        3P           -205.35
        3D            -96.10
        4S            -27.95
        4P            -14.40
=#