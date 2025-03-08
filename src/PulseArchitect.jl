module PulseArchitect

export AudioChain, chain_step, Splitter
export tremolo, delay_reverb, compressor, vibrato, chorus, flanger,
       phaser, distortion, bitcrusher, pitch_shift, lowpass, highpass, bandpass
export generate_melody
export square_wave, saw_wave, triangle_wave, pulse_wave, soft_square, fm_wave, wavy_wave
export play, save#, load

include("AudioChains.jl")
include("Effects.jl")
include("Wave_Generators.jl")

using .AudioChains
using .Effects
using .Wave_Generators
using WAV

"""
generate_melody(notes, bpm; waveform, portamento_time, fs, decay_func, end_pause_mult)

Generates a raw audio signal (melody) that transitions smoothly between notes, inserting silence as specified.

# Arguments
- `notes`: A vector where each element is either:
    - A `String` representing the note (e.g., "C5")—for which default values are assumed, or
    - A tuple `(note, duration, loudness)`:
         - `note`: String (e.g., "C5"); special strings like `"--"`, `"-"`, or `" "` denote a pause.
         - `duration`: Note length multiplier (in beats).
         - `loudness`: Amplitude scaling factor.
- `bpm`: Beats per minute.
  
# Keyword Arguments
- `waveform`: Either a function `phase -> sample` (default is `sine_wave`) or an array of weights used by `default_waveform`.
- `portamento_time`: Crossfade duration in seconds between adjacent notes (if no pause is inserted). If too large, it will be capped so that crossfade indices remain within bounds.
- `fs`: Sample rate (default 44100 Hz).
- `decay_func`: A function `t -> amplitude` applied to each note (default returns 1.0).
- `end_pause_mult`: Multiplier (in beats) for a final appended pause if none was explicitly provided (default 0.2).

# Behavior
- For each note, computes its duration in seconds (using the beat multiplier and bpm) and generates a signal.
- For non-pause notes, the note's frequency is determined via `note_to_freq(note)`, and the waveform is generated accordingly.
- A fade-in and fade-out is applied over a crossfade length (portamento), which is automatically limited so that indices never fall outside the signal boundaries.
- When consecutive non-pause notes are adjacent, they are crossfaded over an effective crossfade window.
- A final pause is appended at the end if none was selected.

Returns a tuple `(overall, fs)` where `overall` is the generated audio signal.
"""
function generate_melody(notes::Vector, bpm::Real;
    waveform = sine_wave,
    portamento_time::Real = 0.05,
    fs::Int = 44100,
    decay_func = t -> 1.0,
    end_pause_mult::Real = 0.2)

    overall = Float64[]
    # Compute the number of samples for portamento and ensure it's at least 1.
    portamento_samples = max(1, round(Int, portamento_time * fs))
    prev_is_pause = false

    for (i, note_item) in enumerate(notes)
        # Expand note_item: if it's a String, use defaults.
        note_str = ""
        duration_mult = 1.0
        loudness = 1.0
        if isa(note_item, String)
            note_str = note_item
        elseif isa(note_item, Tuple)
            note_str = note_item[1]
            duration_mult = length(note_item) >= 2 ? note_item[2] : 1.0
            loudness = length(note_item) >= 3 ? note_item[3] : 1.0
        else
            error("Unsupported note type: $note_item")
        end

        is_pause = note_str in ["--", "-", " "]
        note_duration = duration_mult * (60 / bpm)
        L = max(1, round(Int, note_duration * fs))
        t = range(0, step=1/fs, length=L)
        note_signal = zeros(Float64, L)

        if !is_pause
            freq = note_to_freq(note_str)
            for j in 1:L
                phase = mod(2π * freq * t[j], 2π)
                note_signal[j] = (typeof(waveform) <: Function ? waveform(phase) : default_waveform(phase, waveform))
            end
            # Apply decay to each sample.
            for j in 1:L
                note_signal[j] *= decay_func(t[j])
            end
            note_signal .*= loudness
            # Ensure the fade duration doesn't exceed half the note length.
            fade_samples = portamento_time == 0 ? min(0, div(L,2)) : min(max(2, portamento_samples), div(L,2))

            # Create fade envelopes.
            fade_in = collect(range(0, stop=1, length=fade_samples))
            fade_out = collect(range(1, stop=0, length=fade_samples))

            note_signal[1:fade_samples] .*= fade_in
            note_signal[end-fade_samples+1:end] .*= fade_out
            note_signal[1:fade_samples] .*= fade_in
            note_signal[end-fade_samples+1:end] .*= fade_out
        end

        if i == 1
            append!(overall, note_signal)
        else
            if !prev_is_pause && !is_pause
                # Determine an effective crossfade length that is safe.
                effective_portamento = min(portamento_samples, length(note_signal), length(overall))
                if effective_portamento > 1
                    crossfade = zeros(effective_portamento)
                    for j in 1:effective_portamento
                        weight = (j - 1) / (effective_portamento - 1)
                        overall_index = length(overall) - effective_portamento + j
                        crossfade[j] = (1 - weight) * overall[overall_index] + weight * note_signal[j]
                    end
                    overall[end-effective_portamento+1:end] = crossfade
                    append!(overall, note_signal[effective_portamento+1:end])
                else
                    append!(overall, note_signal)
                end
            else
                append!(overall, note_signal)
            end
        end

        prev_is_pause = is_pause
    end

    # Append final pause: if no explicit pause was provided, add silence for end_pause_mult beats.
    end_pause_samples = round(Int, end_pause_mult * (60 / bpm) * fs)
    append!(overall, zeros(end_pause_samples))

    return (overall, fs)
end

function note_to_freq(note::String)
    m = match(r"^([A-Ga-g])([#b]?)(\d+)$", note)
    if m === nothing
        error("Invalid note format: $note")
    end
    letter = uppercase(m.captures[1])
    accidental = m.captures[2]
    octave = parse(Int, m.captures[3])
    offsets = Dict("C"=> -9, "D"=> -7, "E"=> -5, "F"=> -4, "G"=> -2, "A"=> 0, "B"=> 2)
    n = offsets[letter]
    if accidental == "#"
        n += 1
    elseif accidental == "b"
        n -= 1
    end
    semitone_diff = n + 12*(octave - 4)
    return 440.0 * 2^(semitone_diff/12)
end

# Default waveform: weighted sum of sine waves (overtones).
function default_waveform(phase::Real, weights::AbstractVector{<:Real}=[1.0, 0.6, 0.3, 0.1])
    total = 0.0
    for (i, w) in enumerate(weights)
        total += w * sin(i * phase)
    end
    return total / sum(weights)
end

"""
    play(audio_tuple::Tuple{Vector{Float64}, Int}; volume=1.0)

Plays the given audio signal asynchronously without blocking the main Thread.
The audio samples are scaled by the given volume factor before playback.

- `audio_tuple`: A tuple `(samples, fs)` where:
    - `samples` is a vector of Float64 audio samples (normalized, e.g., in [-1, 1]).
    - `fs` is the sample rate (in Hz).
- `volume`: A multiplier for the audio amplitude (default is 1.0).
"""
function play(audio_tuple::Tuple{Vector{Float64}, Int}; volume=1.0)
    samples, fs = audio_tuple
    # Scale the samples by the volume factor.
    scaled_samples = volume .* samples
    # Launch playback in a separate thread.
    Threads.@spawn wavplay(scaled_samples, fs)
end

"""
    save(audio_tuple::Tuple{Vector{Float64}, Int}, name::AbstractString)

Saves the audio tuple (samples, fs) as a WAV file.
If the provided name does not end with ".wav", it appends the extension.
Returns the full filename.

# Arguments
- `audio_tuple`: A tuple `(samples, fs)` where `samples` is a vector of Float64 and `fs` is the sample rate.
- `name`: The desired base name for the file.
"""
function save(audio_tuple::Tuple{Vector{Float64}, Int}, name::AbstractString)
    filename = endswith(name, ".wav") ? name : name * ".wav"
    WAV.wavwrite(audio_tuple[1], filename, Fs=audio_tuple[2])
    println("Audio saved to $filename")
end

end
