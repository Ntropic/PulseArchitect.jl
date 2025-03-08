# PulseArchitect.jl

**PulseArchitect.jl** is a Julia-based synthesizer package allowing users to easily generate and modify audio signals. Melodies can be constructed from notes and transformed through customizable effect chains, including parallel processing. Audio playback is performed asynchronously, ensuring smooth integration with other code execution.

## Installation

Install the package using:

```julia
using Pkg
Pkg.add("PulseArchitect")
```

## Quick Example

The following example illustrates how to define a melody, generate audio from it using a waveform generator, process it through an effect chain (including parallel effects via a Splitter), and finally play and save the resulting audio:

```julia
# Define a melody as (note, duration, amplitude).
parallel_melody = [
    ("C4", 0.5, 1.0),
    ("E4", 0.5, 1.0),
    ("G4", 0.5, 1.0),
    ("C5", 1.0, 1.0),
    ("--", 3.0, 1.0)
]

# Generate audio from the melody using a waveform.
res_parallel = generate_melody(parallel_melody, 100;
    waveform = soft_square,
    portamento_time = 0.15,
    fs = 96000,
    decay_func = t -> exp(-2*t)
)

# Build two separate effect pipelines.
branch1 = AudioChain(vibrato; strength = 0.04, rate = 7.0)
branch2 = AudioChain(delay_reverb; delay_time = 0.25, feedback = 0.5, mix = 0.7)

# Combine branches into a parallel splitter.
chain = AudioChain(Splitter([0.5, 0.75], [branch1, branch2]))
push!(chain, bandpass; low_cutoff = 200, high_cutoff = 2500)
push!(chain, flanger; depth = 0.02, rate = 0.5, mix = 0.5, feedback = 0.5)

# Process audio through the chain.
res = chain(res_parallel)

# Play audio asynchronously and save to disk.
play(res, volume = 0.125)
println("Parallel branch melody played with Splitter.")
save(res, "parallel_melody")
```

## Authors

- [Michael Schilling](https://github.com/Ntropic)