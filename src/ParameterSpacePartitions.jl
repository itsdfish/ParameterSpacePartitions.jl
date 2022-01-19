module ParameterSpacePartitions
    using Distributions, ConcreteStructs, LinearAlgebra
    using ThreadedIterables, SpecialFunctions 

    export find_partitions,
        adapt!,
        no_adaption!

    export Options,
        Results

    include("structs.jl")
    include("functions.jl")
    include("TestModels.jl")

end
