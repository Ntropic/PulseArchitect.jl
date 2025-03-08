using PulseArchitect
using Test

@testset "PulseArchitect.jl" begin
    # Define a simple melody.
    parallel_melody = [
        ("C4", 0.5, 1.0),
        ("E4", 0.5, 1.0),
        ("G4", 0.5, 1.0),
        ("C5", 1.0, 1.0),
        ("--", 3.0, 1.0)
    ]

    # Generate the melody.
    res_parallel = generate_melody(parallel_melody, 100;
        waveform = soft_square,
        portamento_time = 0.15,
        fs = 96000,
        decay_func = t -> exp(-2*t)
    )

    # Build two effect pipelines.
    branch1 = AudioChain(vibrato; strength = 0.04, rate = 7.0)
    branch2 = AudioChain(delay_reverb; delay_time = 0.25, feedback = 0.5, mix = 0.7)

    # Create a parallel branch using the Splitter.
    # Passing a vector of AudioChains creates a Splitter that processes the input in parallel.
    chain = AudioChain(Splitter([0.5, 0.75], [branch1, branch2]))
    push!(chain, bandpass; low_cutoff = 200, high_cutoff = 2500)
    push!(chain, flanger; depth = 0.02, rate = 0.5, mix = 0.5, feedback = 0.5)

    # Process the generated melody through the effect chain.
    res = chain(res_parallel)

    # (Optionally) Test that the output has the expected structure.
    @test isa(res, Tuple)
    @test length(res) == 2  # Expecting (signal, fs)
    @test isa(res[1], AbstractVector{Float64})
    @test isa(res[2], Integer)

end
