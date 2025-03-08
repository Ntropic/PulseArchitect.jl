module Effects
using FFTW

using Base.Threads
export tremolo, delay_reverb, compressor, vibrato, chorus, flanger,
       phaser, distortion, bitcrusher, pitch_shift, lowpass, highpass, bandpass



"""
tremolo(signal; rate=5.0, depth=0.5)

Modulates the input signal’s amplitude using a low-frequency oscillator.
- `signal`: Input audio array.
- `rate`: (Optional) Modulation frequency in Hz.
- `strength`: (Optional) Modulation strength (0 to 1), where 0 is no modulation.
Returns the modulated signal.
"""
function tremolo(input::Tuple{Vector{Float64},Int}; strength::Real=0.3, rate::Real=4.0)
    signal, fs = input
    L = length(signal)
    t = (0:(L-1)) ./ fs
    modulator = 1 .- strength .* ((1 .- sin.(2π .* rate .* t)) ./ 2)
    return (signal .* modulator, fs)
end

"""
delay_reverb(signal; delay=0.3, feedback=0.4, mix=0.5)

Applies a delay-based reverb by mixing delayed copies with the original signal.
- `signal`: Input audio array.
- `delay`: (Optional) Delay time in seconds.
- `feedback`: (Optional) Proportion of delayed signal fed back (0 to 1).
- `mix`: (Optional) Dry/wet mix ratio (0 to 1).
Returns the reverberated signal.
"""
function delay_reverb(input::Tuple{Vector{Float64},Int};
                      delay_time::Real=0.3, feedback::Real=0.4, mix::Real=0.3)
    signal, fs = input
    L = length(signal)
    delay_samples = round(Int, delay_time * fs)
    output = copy(signal)
    @inbounds for n in (delay_samples+1):L
        output[n] += feedback * output[n - delay_samples]
    end
    output = (1 - mix) .* signal .+ mix .* output
    return (output, fs)
end

"""
compressor(signal; threshold=-20.0, ratio=4.0, attack=0.01, release=0.1)

Compresses the dynamic range by reducing levels above the threshold.
- `signal`: Input audio array.
- `threshold`: (Optional) Level in dB above which compression occurs.
- `ratio`: (Optional) Compression ratio (e.g., 4 means 4:1 compression).
- `attack`: (Optional) Time in seconds to start compression.
- `release`: (Optional) Time in seconds to cease compression.
Returns the compressed signal.
"""
function compressor(input::Tuple{Vector{Float64},Int};
                    threshold::Real=0.5, ratio::Real=4.0, attack::Real=0.01, release::Real=0.1)
    signal, fs = input
    L = length(signal)
    output = similar(signal)
    env = 0.0
    attack_coeff = exp(-1/(attack * fs))
    release_coeff = exp(-1/(release * fs))
    Threads.@threads for n in 1:L
        abs_val = abs(signal[n])
        if abs_val > env
            env = attack_coeff * (env - abs_val) + abs_val
        else
            env = release_coeff * (env - abs_val) + abs_val
        end
        if env > threshold
            gain = threshold + (env - threshold) / ratio
            output[n] = signal[n] * (gain / env)
        else
            output[n] = signal[n]
        end
    end
    return (output, fs)
end


function interp(signal::Vector{Float64}, idx::Real)
    L = length(signal)
    if idx < 1
        return signal[1]
    elseif idx >= L
        return signal[end]
    else
        i0 = floor(Int, idx)
        i1 = min(i0 + 1, L)
        frac = idx - i0
        return (1 - frac)*signal[i0] + frac*signal[i1]
    end
end

"""
vibrato(input::Tuple{Vector{Float64}, Int};
        strength::Real=0.003,
        rate::Real=5.0,
        mod_func::Function = x -> sin(2π * rate * x))

Applies a vibrato effect by modulating the pitch of the input signal.
- `input`: A tuple where the first element is a vector of audio samples and the second is the sample rate.
- `strength`: Maximum modulation deviation (as a fraction of the sample rate).
- `rate`: Modulation frequency in Hz.
- `mod_func`: Function mapping time to modulation factor (default is a sine wave).
Returns the modulated audio as a tuple.
"""
function vibrato(input::Tuple{Vector{Float64}, Int};
                 strength::Real=0.003,
                 rate::Real=5.0,
                 mod_func::Function = x -> sin(2π * rate * x))
    signal, fs = input
    L = length(signal)
    output = zeros(Float64, L)
    t = (0:(L-1)) ./ fs
    max_delay = strength * fs  # maximum delay in samples
    @inbounds Threads.@threads for n in 1:L
        delay = max_delay * mod_func(t[n])
        idx = n - delay
        output[n] = interp(signal, idx)
    end
    return (output, fs)
end

"""
chorus(input::Tuple{Vector{Float64}, Int};
       depth::Real=0.003,
       rate::Real=1.5,
       mix::Real=0.5)

Creates a chorus effect by mixing the original signal with delayed, modulated copies.
- `input`: A tuple containing the audio samples and the sample rate.
- `depth`: Maximum delay modulation (in seconds).
- `rate`: Modulation frequency in Hz.
- `mix`: Dry/wet mix ratio (0 to 1).
Returns the processed audio as a tuple.
"""
function chorus(input::Tuple{Vector{Float64}, Int};
                depth::Real=0.003,
                rate::Real=1.5,
                mix::Real=0.5)
    signal, fs = input
    L = length(signal)
    output = zeros(Float64, L)
    t = (0:(L-1)) ./ fs
    max_delay = depth * fs
    @inbounds Threads.@threads for n in 1:L
        delay = max_delay * sin(2π * rate * t[n])
        delayed = interp(signal, n - delay)
        output[n] = (1 - mix)*signal[n] + mix*delayed
    end
    return (output, fs)
end

"""
flanger(input::Tuple{Vector{Float64}, Int};
        depth::Real=0.002,
        rate::Real=0.25,
        mix::Real=0.5,
        feedback::Real=0.5)

Applies a flanger effect by mixing the signal with a time-delayed version of itself.
- `input`: A tuple with audio samples and sample rate.
- `depth`: Maximum modulation delay in seconds.
- `rate`: Modulation frequency in Hz.
- `mix`: Blend ratio between dry and flanged signals.
- `feedback`: Proportion of the delayed signal fed back into the input.
Returns the flanged audio as a tuple.
"""
function flanger(input::Tuple{Vector{Float64}, Int};
                 depth::Real=0.002,
                 rate::Real=0.25,
                 mix::Real=0.5,
                 feedback::Real=0.5)
    signal, fs = input
    L = length(signal)
    output = zeros(Float64, L)
    t = (0:(L-1)) ./ fs
    max_delay = depth * fs
    @inbounds for n in 1:L
        current_delay = max_delay * (0.5 * (1 + sin(2π * rate * t[n])))
        d = interp(output, n - current_delay)  # using output for feedback
        fb = n > 1 ? output[n-1] : 0.0
        output[n] = (1 - mix)*signal[n] + mix*(d + feedback*fb)
    end
    return (output, fs)
end

"""
phaser(input::Tuple{Vector{Float64}, Int};
       depth::Real=0.7,
       rate::Real=0.5,
       stages::Int=4)

Creates a phaser effect by applying multiple all-pass filters to shift the phase of the signal.
- `input`: A tuple containing the audio samples and sample rate.
- `depth`: Intensity of the phase modulation.
- `rate`: Modulation frequency in Hz.
- `stages`: Number of all-pass filter stages.
Returns the phase-shifted audio as a tuple.
"""
function phaser(input::Tuple{Vector{Float64}, Int};
                depth::Real=0.7,
                rate::Real=0.5,
                stages::Int=4)
    signal, fs = input
    L = length(signal)
    output = copy(signal)
    t = (0:(L-1)) ./ fs
    for stage in 1:stages
        y = zeros(Float64, L)
        @inbounds for n in 1:L
            g = depth * sin(2π * rate * t[n] + 2π*(stage-1)/stages)
            if n == 1
                y[n] = output[n]
            else
                y[n] = -g * output[n] + output[n-1] + g * y[n-1]
            end
        end
        output = y
    end
    return (output, fs)
end

"""
distortion(signal; gain=20.0, threshold=0.8)

Applies distortion by amplifying the signal and clipping its peaks.
- `signal`: Input audio array.
- `gain`: Amplification factor before clipping.
- `threshold`: Clipping threshold (normalized between 0 and 1).
Returns the distorted signal.
"""
function distortion(input::Tuple{Vector{Float64}, Int};
                     gain::Real=2.0,
                     threshold::Real=0.6)
    signal, fs = input
    output = gain .* signal
    output = clamp.(output, -threshold, threshold)
    return (output, fs)
end

"""
bitcrusher(input::Tuple{Vector{Float64}, Int};
           bit_depth::Int=8,
           sample_rate_reduction::Int=1)

Reduces audio fidelity by quantizing the input signal to a lower bit depth and by reducing its effective sample rate.
Assumes that the input samples are normalized in the range [-1, 1].

- `input`: A tuple where the first element is a vector of audio samples and the second is the sample rate.
- `bit_depth`: The target bit depth (e.g., 8 for 8-bit resolution), which determines the number of quantization levels (2^bit_depth).
- `sample_rate_reduction`: The factor by which to reduce the sample rate. For instance, if set to 2, the function processes the signal in blocks of 2 samples,
  replacing each block with a single quantized value.
  
The function first computes the quantization step based on the number of levels, then processes the signal in blocks defined by `sample_rate_reduction`.
For each block, it quantizes a representative sample (the first sample of the block) and replicates this value over the entire block.

Returns a tuple `(crushed_samples, sample_rate)`, where `crushed_samples` is the processed audio signal.
"""
function bitcrusher(input::Tuple{Vector{Float64}, Int}; bit_depth::Int=8, sample_rate_reduction::Int=1)
    samples, fs = input
    levels = 2^bit_depth
    step = 2.0 / (levels - 1)  # Assuming sample range is [-1, 1]
    crushed = similar(samples)
    N = length(samples)
    i = 1
    while i <= N
        # Use the first sample in the block as the representative value.
        value = samples[i]
        # Quantize the value to the nearest level.
        quantized = round(value / step) * step
        # Fill the block with the quantized value.
        j_end = min(i + sample_rate_reduction - 1, N)
        @inbounds @simd for j in i:j_end
            crushed[j] = quantized
        end
        i += sample_rate_reduction
    end
    return (crushed, fs)
end


"""
pitch_shift(input::Tuple{Vector{Float64}, Int}, semitones; window_size=1024, hop_size=256)

Shifts the pitch of the input audio signal by the specified number of semitones using a simplified phase vocoder.
The algorithm performs the following steps:
  1. Computes the short-time Fourier transform (STFT) of the input using a Hann window.
  2. Time-stretches the signal by a factor of 1/factor (where factor = 2^(semitones/12)) by linearly interpolating between STFT frames.
  3. Reconstructs the time-domain signal via the inverse STFT (using overlap-add).
  4. Resamples the resulting signal to match the original length.

- `input`: A tuple `(samples, fs)` where `samples` is a vector of normalized audio samples (Float64, range [-1,1])
  and `fs` is the sample rate.
- `semitones`: The number of semitones to shift the pitch (positive to raise, negative to lower).
- `window_size`: Size of the FFT window (default 1024).
- `hop_size`: Hop size between frames (default 256).

Returns a tuple `(shifted_signal, fs)`.
"""
function pitch_shift(input::Tuple{Vector{Float64}, Int}, semitones; window_size=1024, hop_size=256)


    samples, fs = input
    factor = 2^(semitones/12)      # Frequency scaling factor.
    ratio = 1 / factor           # Time-stretch ratio.

    # Create Hann window.
    win = 0.5 .- 0.5*cos.(2π*(0:window_size-1)/(window_size-1))
    N = length(samples)
    n_frames = div(N - window_size, hop_size)

    # Compute STFT (store as a vector of FFT frames).
    stft = [fft(win .* samples[(i-1)*hop_size+1 : (i-1)*hop_size+window_size]) for i in 1:n_frames]

    # Time-stretch STFT: determine new number of frames.
    n_new = max(1, Int(floor(n_frames * ratio)))
    stft_new = Vector{Vector{ComplexF64}}(undef, n_new)
    for i in 1:n_new
        pos = (i-1) / ratio + 1  # fractional frame index in original STFT
        i0 = clamp(floor(Int, pos), 1, n_frames)
        i1 = clamp(i0 + 1, 1, n_frames)
        frac = pos - i0
        stft_new[i] = (1 - frac) .* stft[i0] .+ frac .* stft[i1]
    end

    # Inverse STFT via overlap-add.
    y_len = (n_new-1)*hop_size + window_size
    y = zeros(Float64, y_len)
    win_sum = zeros(Float64, y_len)
    for i in 1:n_new
        frame = real(ifft(stft_new[i]))
        start_idx = (i-1)*hop_size + 1
        y[start_idx:start_idx+window_size-1] .+= frame .* win
        win_sum[start_idx:start_idx+window_size-1] .+= win.^2
    end
    # Avoid division by zero.
    nz = win_sum .> 1e-6
    y[nz] ./= win_sum[nz]

    # Resample y to original length using linear interpolation.
    L = N
    resampled = similar(y, L)
    for i in 1:L
        t = (i-1) / (L-1) * (length(y)-1) + 1
        i0 = clamp(floor(Int, t), 1, length(y))
        i1 = clamp(i0+1, 1, length(y))
        frac = t - i0
        resampled[i] = (1-frac)*y[i0] + frac*y[i1]
    end

    return (resampled, fs)
end

"""
lowpass(input::Tuple{Vector{Float64}, Int}; cutoff::Real)

Applies a first‐order lowpass filter to the input signal using an exponential moving average.
- `input`: Tuple `(samples, fs)` where `samples` is a vector of normalized audio samples (range [-1, 1])
  and `fs` is the sample rate in Hz.
- `cutoff`: Cutoff frequency in Hz. Frequencies above this are increasingly attenuated.
Returns a tuple `(filtered_samples, fs)`.
"""
function lowpass(input::Tuple{Vector{Float64}, Int}; cutoff::Real)
    samples, fs = input
    dt = 1 / fs
    RC = 1 / (2π * cutoff)
    α = dt / (RC + dt)
    y = similar(samples)
    y[1] = samples[1]
    for i in 2:length(samples)
        y[i] = y[i-1] + α * (samples[i] - y[i-1])
    end
    return (y, fs)
end

"""
highpass(input::Tuple{Vector{Float64}, Int}; cutoff::Real)

Applies a first‐order highpass filter to the input signal using a simple recursive formula.
- `input`: Tuple `(samples, fs)` where `samples` is a vector of normalized audio samples (range [-1, 1])
  and `fs` is the sample rate in Hz.
- `cutoff`: Cutoff frequency in Hz. Frequencies below this are increasingly attenuated.
Returns a tuple `(filtered_samples, fs)`.
"""
function highpass(input::Tuple{Vector{Float64}, Int}; cutoff::Real)
    samples, fs = input
    dt = 1 / fs
    RC = 1 / (2π * cutoff)
    α = RC / (RC + dt)
    y = similar(samples)
    y[1] = samples[1]
    for i in 2:length(samples)
        y[i] = α * (y[i-1] + samples[i] - samples[i-1])
    end
    return (y, fs)
end

"""
bandpass(input::Tuple{Vector{Float64}, Int}; low_cutoff::Real, high_cutoff::Real)

Creates a bandpass filter by cascading a highpass and a lowpass filter.
- `input`: Tuple `(samples, fs)` where `samples` is a vector of normalized audio samples (range [-1, 1])
  and `fs` is the sample rate in Hz.
- `low_cutoff`: Lower cutoff frequency in Hz (highpass threshold).
- `high_cutoff`: Upper cutoff frequency in Hz (lowpass threshold).
Returns a tuple `(filtered_samples, fs)` where only frequencies between `low_cutoff` and `high_cutoff` are preserved.
"""
function bandpass(input::Tuple{Vector{Float64}, Int}; low_cutoff::Real, high_cutoff::Real)
    # First, remove frequencies below low_cutoff.
    hp_output = highpass(input; cutoff=low_cutoff)
    # Then, remove frequencies above high_cutoff.
    bp_output = lowpass(hp_output; cutoff=high_cutoff)
    return bp_output
end


end # module Effects
