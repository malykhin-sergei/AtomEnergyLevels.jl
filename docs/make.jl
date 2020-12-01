using Documenter
using AtomEnergyLevels

DocMeta.setdocmeta!(AtomEnergyLevels, :DocTestSetup, :(using AtomEnergyLevels); recursive=true)

makedocs(
    sitename = "AtomEnergyLevels",
    format = Documenter.HTML(),
    modules = [AtomEnergyLevels]
)

deploydocs(
    repo = "github.com/malykhin-sergei/AtomEnergyLevels.jl.git",
)
