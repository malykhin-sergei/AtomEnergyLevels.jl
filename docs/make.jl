using Documenter
using AtomEnergyLevels

makedocs(
    sitename = "AtomEnergyLevels",
    format = Documenter.HTML(),
    modules = [AtomEnergyLevels]
)

deploydocs(
    repo = "github.com/malykhin-sergei/AtomEnergyLevels.jl.git",
)
