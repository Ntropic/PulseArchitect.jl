# PulseArchitect.jl

[![Build Status](https://github.com/Ntropic/PulseArchitect.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Ntropic/PulseArchitect.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://Ntropic.github.io/PulseArchitect.jl/)

**PulseArchitect.jl** uis a Julia based Synthesizer, alolowing the generation and modification of audio signals from melodies and by defining process Chains. 
The audio is played via `play` without blocking the execution of further code. 

## Installation
Install it using 
```julia
using Pkg
Pkg.add("PulseArchitect")
```

## Usage 
A quick example, first defining a melody, the wave_generator for it and then an effect Chain including the parallel execution of two effects, 
and finally playing the sound and saving it to disk.

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
