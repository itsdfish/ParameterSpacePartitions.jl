module ParameterSpacePartitions
    using Distributions, ConcreteStructs, LinearAlgebra
    using ThreadedIterables, SpecialFunctions 

    export find_partitions,
        adapt!,
        no_adaption!,
        estimate_volume,
        bias_correction

    export Options,
        Results

    include("structs.jl")
    include("sampler.jl")
    include("volume.jl")
    include("TestModels.jl")

end
