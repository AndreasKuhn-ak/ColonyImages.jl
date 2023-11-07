using ColonyImages
using Documenter

DocMeta.setdocmeta!(ColonyImages, :DocTestSetup, :(using ColonyImages); recursive=true)

makedocs(;
    modules=[ColonyImages],
    authors="AndreasKuhn-ak <andreaskuhn92@gmx.net> and contributors",
    repo="https://github.com/AndreasKuhn-ak/ColonyImages.jl/blob/{commit}{path}#{line}",
    sitename="ColonyImages.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://AndreasKuhn-ak.github.io/ColonyImages.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AndreasKuhn-ak/ColonyImages.jl",
    devbranch="master",
)
