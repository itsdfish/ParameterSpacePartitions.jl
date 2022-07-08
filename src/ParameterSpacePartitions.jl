module ParameterSpacePartitions
    using Requires, ProgressMeter, Distributions, ConcreteStructs, LinearAlgebra
    using ThreadedIterables, SpecialFunctions, ComponentArrays

    export find_partitions,
        adapt!,
        no_adaption!,
        estimate_volume,
        bias_correction,
        psp_slices

    export Options,
        Results

    include("init.jl")
    include("structs.jl")
    include("sampler.jl")
    include("volume.jl")
    include("TestModels.jl")

end
