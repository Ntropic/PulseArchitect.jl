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
using PulseArchitect

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

# Build two effect pipelines:
branch1 = AudioChain(vibrato; strength = 0.04, rate = 7.0)
branch2 = AudioChain(delay_reverb; delay_time = 0.25, feedback = 0.5, mix = 0.7)

# Create a parallel branch using the Splitter.
# Passing a vector of AudioChains creates a Splitter that processes the input in parallel.
chain = AudioChain(Splitter([0.5, 0.5], [branch1, branch2]))
push!(chain, bandpass; low_cutoff=200, high_cutoff=2500)

res = chain(res_parallel)
# Process and play the sound asynchronously.
play(res, volume = 0.5)
println("Parallel branch melody played with Splitter.")
save(res, "parallel_melody")
```

## Authors
- [Michael Schilling](https://github.com/Ntropic)