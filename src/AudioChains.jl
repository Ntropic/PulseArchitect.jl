module AudioChains

export AudioChain, Splitter

import Base: push!, append!

# --- Data Structures ---


struct ChainStep
    f::Function
    kwargs::NamedTuple
end

# Helper to create a ChainStep.
chain_step(f::Function; kwargs...) = ChainStep(f, (; kwargs...))

"""
    AudioChain

Represents an audio processing pipeline composed of multiple steps. Each step can be either:
  - a `ChainStep` (i.e. a single sequential processing operation), or
  - a `Vector` of `AudioChain` (representing a parallel branch where the input is split among subchains).

The `AudioChain` type supports various input options via its constructors:
  
1. **Empty Chain:**  
   Create an empty audio chain with no processing steps.
   ```julia
   chain = AudioChain()
   ```

2. **Iterable of Steps:**  
   Build an audio chain from an iterable of arguments where each element can be:
     - a `ChainStep`,
     - an `AudioChain` (in which case its steps are flattened into the parent chain),
     - or a `Function` (which is automatically wrapped into a `ChainStep` using `chain_step`).
   ```julia
   chain = AudioChain(f, some_subchain, g)
   ```

3. **Single Function with Keyword Arguments:**  
   Shortcut constructor that accepts a single function and optional keyword arguments.
   ```julia
   chain = AudioChain(f; multiplier=2)
   ```

The `AudioChain` is also callable as a function. When invoked on an input (wave_vector::Vector{Real}, sampling_rate::Int), it processes the input
through its sequence of steps.
"""
struct AudioChain
    steps::Vector{ChainStep}
    # Default constructor for an empty AudioChain.
    function AudioChain()
        new(ChainStep[])
    end
    # Construct an AudioChain from an iterable of arguments.
    # Each argument can be:
    #   - a ChainStep,
    #   - an AudioChain (its steps will be flattened),
    #   - or a Function (which is wrapped into a ChainStep).
    function AudioChain(args...)
        vec = ChainStep[]
        for arg in args
            if arg isa ChainStep
                push!(vec, arg)
            elseif arg isa AudioChain
                # Flatten the steps of the subchain into this chain.
                append!(vec, arg.steps)
            elseif arg isa Function
                push!(vec, chain_step(arg))
            elseif arg isa Splitter # Special case for Splitter objects.
                push!(vec, chain_step(x -> arg(x)))
            else
                error("Unsupported argument type: $(typeof(arg))")
            end
        end
        new(vec)
    end
    # Shortcut constructor for a single function with optional kwargs.
    function AudioChain(f::Function; kwargs...)
        return new([chain_step(f; kwargs...)])
    end
end

# Make an AudioChain callable.
function (chain::AudioChain)(input)
    result = input
    for step in chain.steps
        if step isa ChainStep
            result = step.f(result; step.kwargs...)
        else
            error("Unsupported chain step type: $(typeof(step))")
        end
    end
    return result
end

# --- Overloaded push! Methods ---

# 1. When pushing a function (with or without keyword arguments).
function push!(chain::AudioChain, f::Function; kwargs...)
    push!(chain.steps, ChainStep(f, (; kwargs...)))
    return chain
end

# 2. When pushing a tuple: (function, kwargs::NamedTuple).
function push!(chain::AudioChain, tpl::Tuple)
    if length(tpl) == 1
        return push!(chain, tpl[1])
    elseif length(tpl) == 2 && tpl[2] isa NamedTuple
        return push!(chain, tpl[1]; tpl[2]...)
    else
        error("Tuple must be (function) or (function, NamedTuple)")
    end
end

# 3. When pushing an AudioChain, flatten its steps into the parent.
function push!(chain::AudioChain, subchain::AudioChain)
    append!(chain.steps, subchain.steps)
    return chain
end

# 4. When pushing a vector of AudioChain objects, add them as a parallel branch.
function push!(chain::AudioChain, chains::Vector{AudioChain})
    push!(chain.steps, chains)
    return chain
end

# --- Overloaded append! Methods ---
# For convenience, append! behaves exactly like push! for our supported types.
function append!(chain::AudioChain, x)
    push!(chain, x)
    return chain
end


"""
    Splitter(weights::Vector{<:Real}, pipelines::Vector)

A splitter that processes an input signal through multiple pipelines concurrently, and then merges the results in a weighted manner.

# Fields
- `weights`: A vector of amplitude weights for each pipeline.
- `pipelines`: A vector of pipelines, each being either an `AudioChain` or a function that accepts an input signal and returns a vector of samples.

# Constructor
The constructor accepts a vector of weights and a vector of pipelines (which can be functions or `AudioChain` objects). It verifies that both vectors have the same length.

# Functor
Calling an instance of `Splitter` with an input signal will:
1. Execute all pipelines concurrently (in parallel) on the input signal.
2. Scale each pipeline's output by its corresponding weight.
3. Sum the weighted outputs element‐wise to produce the final merged output.

Assumes that all pipelines return arrays of samples of the same length.
"""
struct Splitter
    weights::Vector{<:Real}
    pipelines::Vector{Union{AudioChain, Function}}
    function Splitter(weights::Vector{<:Real}, pipelines::Vector)
        @assert length(weights) == length(pipelines) "The number of weights must match the number of pipelines"
        pipelines_converted = Union{AudioChain, Function}[p for p in pipelines]
        new(weights, pipelines_converted)
    end
end



# Define the functor to process the input signal.
function (s::Splitter)(input)
    # Run each pipeline concurrently, collecting full tuples.
    results = Vector{Tuple{Vector{Float64}, Int}}(undef, length(s.pipelines))
    @sync for i in eachindex(s.pipelines)
        @async begin
            res = s.pipelines[i](input)
            results[i] = res
        end
    end
    # Assume fs is the same for each branch; get fs from the first branch.
    fs = results[1][2]
    # Extract the signal from the first branch.
    merged_signal = similar(results[1][1])
    weights = s.weights
    sig0 = results[1][1]
    @inbounds @simd for j in eachindex(sig0)
        merged_signal[j] = sig0[j] * weights[1]
    end
    # Add weighted outputs of the remaining branches.
    for i in 2:length(results)
        branch_signal = results[i][1]
        @inbounds @simd for j in eachindex(branch_signal)
            merged_signal[j] += branch_signal[j] * weights[i]
        end
    end
    return (merged_signal, fs)
end



function Base.push!(s::Splitter, weight::Real, pipeline::Union{AudioChain, Function})
    push!(s.weights, weight)
    push!(s.pipelines, pipeline)
    return s
end

end  # module AudioChains