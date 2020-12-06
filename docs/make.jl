using Documenter
using AtomEnergyLevels

DocMeta.setdocmeta!(AtomEnergyLevels, :DocTestSetup, :(using AtomEnergyLevels); recursive=true)

makedocs(
    sitename = "AtomEnergyLevels",
    format = Documenter.HTML(),
    modules = [AtomEnergyLevels]
)

if "deploy" in ARGS
  include("../../faketravis.jl")
end

deploydocs(
    repo = "github.com/malykhin-sergei/AtomEnergyLevels.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
)
