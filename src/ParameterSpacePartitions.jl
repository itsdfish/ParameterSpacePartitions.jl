module ParameterSpacePartitions
    using Distributions, ConcreteStructs, LinearAlgebra
    using ThreadedIterables 

    export random_position,
        generate_proposal,
        find_partitions,
        eval_patterns,
        initialize,
        update_position!


    export Chain, 
        Options,
        Results

    include("structs.jl")
    include("functions.jl")

end
