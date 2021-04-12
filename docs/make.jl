using MicMods
using Documenter

DocMeta.setdocmeta!(MicMods, :DocTestSetup, :(using MicMods); recursive=true)

makedocs(;
    modules=[MicMods],
    authors="Thomas Wutzler <progtw@arcor.de> and contributors",
    repo="https://github.com/bgctw/MicMods.jl/blob/{commit}{path}#{line}",
    sitename="MicMods.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bgctw.github.io/MicMods.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bgctw/MicMods.jl",
    devbranch="main",
)
