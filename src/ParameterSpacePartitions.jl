module ParameterSpacePartitions
    using Requires, Distributions, ConcreteStructs, LinearAlgebra
    using ThreadedIterables, SpecialFunctions, ComponentArrays
    #using ProgressMeter 

    export find_partitions,
        adapt!,
        no_adaption!,
        estimate_volume,
        bias_correction,
        psp_slices,
        intersects

    export Options,
        Results

    include("init.jl")
    include("structs.jl")
    include("sampler.jl")
    include("volume.jl")
    include("intersection_test.jl")
    include("TestModels.jl")

end
