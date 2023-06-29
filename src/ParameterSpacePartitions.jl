module ParameterSpacePartitions
    using ComponentArrays
    using ConcreteStructs
    using Distributions
    using LinearAlgebra
    using ProgressMeter
    using Requires
    using SpecialFunctions
    using ThreadedIterables

    export find_partitions
    export adapt!
    export no_adaption!
    export estimate_volume
    export psp_slices

    export Options
    export Results

    include("init.jl")
    include("structs.jl")
    include("sampler.jl")
    include("volume.jl")
    include("TestModels.jl")

end
