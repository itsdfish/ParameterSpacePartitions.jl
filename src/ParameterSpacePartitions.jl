module ParameterSpacePartitions
    using Distributions, ConcreteStructs, LinearAlgebra
    using ThreadedIterables 

    export find_partitions

    export Options,
        Results

    include("structs.jl")
    include("functions.jl")

end
