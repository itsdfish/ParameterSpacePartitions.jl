using Documenter
using ParameterSpacePartitions


makedocs(
    sitename = "ParameterSpacePartitions",
    format = Documenter.HTML(
        assets = [
            asset(
                "https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap",
                class = :css,
            ),
        ],
        collapselevel = 1,
    ),

    modules = [ParameterSpacePartitions],
    pages = ["home" => "index.md",
            "examples" => ["example 1" => "example1.md"],
            "api" => "api.md"]
)

deploydocs(
    repo = "github.com/itsdfish/ParameterSpacePartitions.jl.git",
)