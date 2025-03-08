var documenterSearchIndex = {"docs":
[{"location":"api/pulse_architect/#Pulse-Architect","page":"Pulse Architect","title":"Pulse Architect","text":"","category":"section"},{"location":"api/pulse_architect/","page":"Pulse Architect","title":"Pulse Architect","text":"In order to apply effects or chains of effects to an audio signal, we first have to create an audio signal. One way to do so is the generate_melody function.  To use it, specify a melody as a Vector of Tuples. Each Tuple contains three elements: the note name (e.g., \"C4\"), its duration in beats, and its amplitude. Generate melody can use this to generate an initial audio signal, using optional parameters, such as a waveform generator, the sample rate, portamento time (smoothly transitioning between notes), a decay function (for plucked note sound) and more.","category":"page"},{"location":"api/pulse_architect/","page":"Pulse Architect","title":"Pulse Architect","text":"generate_melody","category":"page"},{"location":"api/pulse_architect/#PulseArchitect.generate_melody","page":"Pulse Architect","title":"PulseArchitect.generate_melody","text":"generatemelody(notes, bpm; waveform, portamentotime, fs, decayfunc, endpause_mult)\n\nGenerates a raw audio signal (melody) that transitions smoothly between notes, inserting silence as specified.\n\nArguments\n\nnotes: A vector where each element is either:\nA String representing the note (e.g., \"C5\")—for which default values are assumed, or\nA tuple (note, duration, loudness):\nnote: String (e.g., \"C5\"); special strings like \"--\", \"-\", or \" \" denote a pause.\nduration: Note length multiplier (in beats).\nloudness: Amplitude scaling factor.\nbpm: Beats per minute.\n\nKeyword Arguments\n\nwaveform: Either a function phase -> sample (default is sine_wave) or an array of weights used by default_waveform.\nportamento_time: Crossfade duration in seconds between adjacent notes (if no pause is inserted). If too large, it will be capped so that crossfade indices remain within bounds.\nfs: Sample rate (default 44100 Hz).\ndecay_func: A function t -> amplitude applied to each note (default returns 1.0).\nend_pause_mult: Multiplier (in beats) for a final appended pause if none was explicitly provided (default 0.2).\n\nBehavior\n\nFor each note, computes its duration in seconds (using the beat multiplier and bpm) and generates a signal.\nFor non-pause notes, the note's frequency is determined via note_to_freq(note), and the waveform is generated accordingly.\nA fade-in and fade-out is applied over a crossfade length (portamento), which is automatically limited so that indices never fall outside the signal boundaries.\nWhen consecutive non-pause notes are adjacent, they are crossfaded over an effective crossfade window.\nA final pause is appended at the end if none was selected.\n\nReturns a tuple (overall, fs) where overall is the generated audio signal.\n\n\n\n\n\n","category":"function"},{"location":"api/pulse_architect/","page":"Pulse Architect","title":"Pulse Architect","text":"Every audio object is a tuple of the form (audio_signal, fs) where audio_signal is an array of floats representing the audio waveform and fs is the sample rate in Hz. This allows us to manipulate the audio signal using various functions provided by PulseArchitect.","category":"page"},{"location":"api/pulse_architect/","page":"Pulse Architect","title":"Pulse Architect","text":"We can play an audio Tuple via the play function. This does not block the execution of other commands, as it is started in a separate thread.  Audio Signals can be saved via the save function.","category":"page"},{"location":"api/pulse_architect/","page":"Pulse Architect","title":"Pulse Architect","text":"play\nsave","category":"page"},{"location":"api/pulse_architect/#PulseArchitect.play","page":"Pulse Architect","title":"PulseArchitect.play","text":"play(audio_tuple::Tuple{Vector{Float64}, Int}; volume=1.0)\n\nPlays the given audio signal asynchronously without blocking the main Thread. The audio samples are scaled by the given volume factor before playback.\n\naudio_tuple: A tuple (samples, fs) where:\nsamples is a vector of Float64 audio samples (normalized, e.g., in [-1, 1]).\nfs is the sample rate (in Hz).\nvolume: A multiplier for the audio amplitude (default is 1.0).\n\n\n\n\n\n","category":"function"},{"location":"api/pulse_architect/#PulseArchitect.save","page":"Pulse Architect","title":"PulseArchitect.save","text":"save(audio_tuple::Tuple{Vector{Float64}, Int}, name::AbstractString)\n\nSaves the audio tuple (samples, fs) as a WAV file. If the provided name does not end with \".wav\", it appends the extension. Returns the full filename.\n\nArguments\n\naudio_tuple: A tuple (samples, fs) where samples is a vector of Float64 and fs is the sample rate.\nname: The desired base name for the file.\n\n\n\n\n\n","category":"function"},{"location":"api/wave_generators/#Wave-Generators","page":"Wave Generators","title":"Wave Generators","text":"","category":"section"},{"location":"api/wave_generators/","page":"Wave Generators","title":"Wave Generators","text":"For your convenience we provide a few wave generators, which you can use with generate_melody. These are simple functions, feel free to define your own. Here is a list of available wave generators:","category":"page"},{"location":"api/wave_generators/","page":"Wave Generators","title":"Wave Generators","text":"square_wave\nsaw_wave\ntriangle_wave\npulse_wave\nsoft_square\nfm_wave\nwavy_wave","category":"page"},{"location":"api/wave_generators/#PulseArchitect.Wave_Generators.square_wave","page":"Wave Generators","title":"PulseArchitect.Wave_Generators.square_wave","text":"square_wave(phase::Real)\n\nReturns 1.0 for phase < π, else -1.0.\n\n\n\n\n\n","category":"function"},{"location":"api/wave_generators/#PulseArchitect.Wave_Generators.saw_wave","page":"Wave Generators","title":"PulseArchitect.Wave_Generators.saw_wave","text":"saw_wave(phase::Real)\n\nReturns a sawtooth wave in [-1, 1] from phase.\n\n\n\n\n\n","category":"function"},{"location":"api/wave_generators/#PulseArchitect.Wave_Generators.triangle_wave","page":"Wave Generators","title":"PulseArchitect.Wave_Generators.triangle_wave","text":"triangle_wave(phase::Real)\n\nReturns a triangle wave in [-1, 1] from phase.\n\n\n\n\n\n","category":"function"},{"location":"api/wave_generators/#PulseArchitect.Wave_Generators.pulse_wave","page":"Wave Generators","title":"PulseArchitect.Wave_Generators.pulse_wave","text":"pulse_wave(phase::Real; duty=0.3)\n\nReturns a pulse wave with given duty cycle.\n\n\n\n\n\n","category":"function"},{"location":"api/wave_generators/#PulseArchitect.Wave_Generators.soft_square","page":"Wave Generators","title":"PulseArchitect.Wave_Generators.soft_square","text":"soft_square(phase::Real)\n\nReturns a soft square wave via tanh(5*sin(phase)).\n\n\n\n\n\n","category":"function"},{"location":"api/wave_generators/#PulseArchitect.Wave_Generators.fm_wave","page":"Wave Generators","title":"PulseArchitect.Wave_Generators.fm_wave","text":"fm_wave(phase::Real)\n\nReturns a frequency-modulated wave: sin(phase + 0.5sin(3phase)).\n\n\n\n\n\n","category":"function"},{"location":"api/wave_generators/#PulseArchitect.Wave_Generators.wavy_wave","page":"Wave Generators","title":"PulseArchitect.Wave_Generators.wavy_wave","text":"wavy_wave(phase::Real)\n\nReturns a composite wave: sin(phase) + 0.5sin(2phase + 0.5).\n\n\n\n\n\n","category":"function"},{"location":"#PulseArchitect.jl","page":"Home","title":"PulseArchitect.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PulseArchitect.jl is a Julia-based synthesizer package allowing users to easily generate and modify audio signals. Melodies can be constructed from notes and transformed through customizable effect chains, including parallel processing. Audio playback is performed asynchronously, ensuring smooth integration with other code execution.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Install the package using:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"PulseArchitect\")","category":"page"},{"location":"#Quick-Example","page":"Home","title":"Quick Example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The following example illustrates how to define a melody, generate audio from it using a waveform generator, process it through an effect chain (including parallel effects via a Splitter), and finally play and save the resulting audio:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using PulseArchitect\n\n# Define a simple melody.\nparallel_melody = [\n    (\"C4\", 0.5, 1.0),\n    (\"E4\", 0.5, 1.0),\n    (\"G4\", 0.5, 1.0),\n    (\"C5\", 1.0, 1.0),\n    (\"--\", 3.0, 1.0)\n]\n\n# Generate the melody.\nres_parallel = generate_melody(parallel_melody, 100;\n    waveform = soft_square,\n    portamento_time = 0.15,\n    fs = 96000,\n    decay_func = t -> exp(-2*t)\n)\n\n# Build two effect pipelines:\nbranch1 = AudioChain(vibrato; strength = 0.04, rate = 7.0)\nbranch2 = AudioChain(delay_reverb; delay_time = 0.25, feedback = 0.5, mix = 0.7)\n\n# Create a parallel branch using the Splitter.\n# Passing a vector of AudioChains creates a Splitter that processes the input in parallel.\nchain = AudioChain(Splitter([0.5, 0.5], [branch1, branch2]))\npush!(chain, bandpass; low_cutoff=200, high_cutoff=2500)\n\nres = chain(res_parallel)\n# Process and play the sound asynchronously.\nplay(res, volume = 0.5)\nprintln(\"Parallel branch melody played with Splitter.\")\nsave(res, \"parallel_melody\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"This code snippet creates the audio file ","category":"page"},{"location":"","page":"Home","title":"Home","text":"<audio controls>   <source src=\"assets/parallel_melody.wav\" type=\"audio/wav\">   Your browser does not support the audio element. </audio>","category":"page"},{"location":"#Functions","page":"Home","title":"Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The provided functions are split into ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pulse Architect\nEffects\nAudio Chains\nWave Generators","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Michael Schilling","category":"page"},{"location":"api/audio_chains/#Audio-Chains-and-Splitters","page":"Audio Chains","title":"Audio Chains and Splitters","text":"","category":"section"},{"location":"api/audio_chains/","page":"Audio Chains","title":"Audio Chains","text":"AudioChains and Splitters allow us to define sequential and parallel effect processing pipelines for audio signals. They are particularly useful in scenarios where multiple effects need to be applied or when different parts of the signal need to be processed independently.","category":"page"},{"location":"api/audio_chains/#Audio-Chain","page":"Audio Chains","title":"Audio Chain","text":"","category":"section"},{"location":"api/audio_chains/","page":"Audio Chains","title":"Audio Chains","text":"AudioChains are defined as","category":"page"},{"location":"api/audio_chains/","page":"Audio Chains","title":"Audio Chains","text":"AudioChain","category":"page"},{"location":"api/audio_chains/#PulseArchitect.AudioChains.AudioChain","page":"Audio Chains","title":"PulseArchitect.AudioChains.AudioChain","text":"AudioChain\n\nRepresents an audio processing pipeline composed of multiple steps. Each step can be either:\n\na ChainStep (i.e. a single sequential processing operation), or\na Vector of AudioChain (representing a parallel branch where the input is split among subchains).\n\nThe AudioChain type supports various input options via its constructors:\n\nEmpty Chain:   Create an empty audio chain with no processing steps.\nchain = AudioChain()\nIterable of Steps:   Build an audio chain from an iterable of arguments where each element can be:\na ChainStep,\nan AudioChain (in which case its steps are flattened into the parent chain),\nor a Function (which is automatically wrapped into a ChainStep using chain_step).\nchain = AudioChain(f, some_subchain, g)\nSingle Function with Keyword Arguments:   Shortcut constructor that accepts a single function and optional keyword arguments.\nchain = AudioChain(f; multiplier=2)\n\nThe AudioChain is also callable as a function. When invoked on an input (wavevector::Vector{Real}, samplingrate::Int), it processes the input through its sequence of steps.\n\n\n\n\n\n","category":"type"},{"location":"api/audio_chains/#Splitter","page":"Audio Chains","title":"Splitter","text":"","category":"section"},{"location":"api/audio_chains/","page":"Audio Chains","title":"Audio Chains","text":"Splitter allows us to split an audio signal into multiple channels and process each channel separately. It is useful when we need to apply different effects to different parts of the signal.","category":"page"},{"location":"api/audio_chains/","page":"Audio Chains","title":"Audio Chains","text":"Splitter","category":"page"},{"location":"api/audio_chains/#PulseArchitect.AudioChains.Splitter","page":"Audio Chains","title":"PulseArchitect.AudioChains.Splitter","text":"Splitter(weights::Vector{<:Real}, pipelines::Vector)\n\nA splitter that processes an input signal through multiple pipelines concurrently, and then merges the results in a weighted manner.\n\nFields\n\nweights: A vector of amplitude weights for each pipeline.\npipelines: A vector of pipelines, each being either an AudioChain or a function that accepts an input signal and returns a vector of samples.\n\nConstructor\n\nThe constructor accepts a vector of weights and a vector of pipelines (which can be functions or AudioChain objects). It verifies that both vectors have the same length.\n\nFunctor\n\nCalling an instance of Splitter with an input signal will:\n\nExecute all pipelines concurrently (in parallel) on the input signal.\nScale each pipeline's output by its corresponding weight.\nSum the weighted outputs element‐wise to produce the final merged output.\n\nAssumes that all pipelines return arrays of samples of the same length.\n\n\n\n\n\n","category":"type"},{"location":"api/effects/#Effects","page":"Effects","title":"Effects","text":"","category":"section"},{"location":"api/effects/","page":"Effects","title":"Effects","text":"The pulse Architect provides you with a wide range of audio processing effects that can be applied to your audio signals. As input they take a Tuple of the wave Vector and the sample rate, as output they return a Tuple of the modified wave vector and the new sample rate. ","category":"page"},{"location":"api/effects/","page":"Effects","title":"Effects","text":"A audio effect can be directly called for example via tremolo((wave, fs); strength=0.3, rate=4.0) or indirectly, by adding it to an AudioChain via ","category":"page"},{"location":"api/effects/","page":"Effects","title":"Effects","text":"chain = AudioChain(tremolo; strength=0.3, rate=4.0)\n# or by pushing it into an existing chain\npush!(chain, tremolo; strength=0.5, rate=2.0)","category":"page"},{"location":"api/effects/","page":"Effects","title":"Effects","text":"Below is a list of the available effects:","category":"page"},{"location":"api/effects/","page":"Effects","title":"Effects","text":"tremolo\ndelay_reverb\ncompressor\nvibrato\nchorus\nflanger\nphaser\ndistortion\nbitcrusher\npitch_shift\nlowpass\nhighpass\nbandpass","category":"page"},{"location":"api/effects/#PulseArchitect.Effects.tremolo","page":"Effects","title":"PulseArchitect.Effects.tremolo","text":"tremolo(signal; rate=5.0, depth=0.5)\n\nModulates the input signal’s amplitude using a low-frequency oscillator.\n\nsignal: Input audio array.\nrate: (Optional) Modulation frequency in Hz.\nstrength: (Optional) Modulation strength (0 to 1), where 0 is no modulation.\n\nReturns the modulated signal.\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.delay_reverb","page":"Effects","title":"PulseArchitect.Effects.delay_reverb","text":"delay_reverb(signal; delay=0.3, feedback=0.4, mix=0.5)\n\nApplies a delay-based reverb by mixing delayed copies with the original signal.\n\nsignal: Input audio array.\ndelay: (Optional) Delay time in seconds.\nfeedback: (Optional) Proportion of delayed signal fed back (0 to 1).\nmix: (Optional) Dry/wet mix ratio (0 to 1).\n\nReturns the reverberated signal.\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.compressor","page":"Effects","title":"PulseArchitect.Effects.compressor","text":"compressor(signal; threshold=-20.0, ratio=4.0, attack=0.01, release=0.1)\n\nCompresses the dynamic range by reducing levels above the threshold.\n\nsignal: Input audio array.\nthreshold: (Optional) Level in dB above which compression occurs.\nratio: (Optional) Compression ratio (e.g., 4 means 4:1 compression).\nattack: (Optional) Time in seconds to start compression.\nrelease: (Optional) Time in seconds to cease compression.\n\nReturns the compressed signal.\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.vibrato","page":"Effects","title":"PulseArchitect.Effects.vibrato","text":"vibrato(input::Tuple{Vector{Float64}, Int};         strength::Real=0.003,         rate::Real=5.0,         mod_func::Function = x -> sin(2π * rate * x))\n\nApplies a vibrato effect by modulating the pitch of the input signal.\n\ninput: A tuple where the first element is a vector of audio samples and the second is the sample rate.\nstrength: Maximum modulation deviation (as a fraction of the sample rate).\nrate: Modulation frequency in Hz.\nmod_func: Function mapping time to modulation factor (default is a sine wave).\n\nReturns the modulated audio as a tuple.\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.chorus","page":"Effects","title":"PulseArchitect.Effects.chorus","text":"chorus(input::Tuple{Vector{Float64}, Int};        depth::Real=0.003,        rate::Real=1.5,        mix::Real=0.5)\n\nCreates a chorus effect by mixing the original signal with delayed, modulated copies.\n\ninput: A tuple containing the audio samples and the sample rate.\ndepth: Maximum delay modulation (in seconds).\nrate: Modulation frequency in Hz.\nmix: Dry/wet mix ratio (0 to 1).\n\nReturns the processed audio as a tuple.\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.flanger","page":"Effects","title":"PulseArchitect.Effects.flanger","text":"flanger(input::Tuple{Vector{Float64}, Int};         depth::Real=0.002,         rate::Real=0.25,         mix::Real=0.5,         feedback::Real=0.5)\n\nApplies a flanger effect by mixing the signal with a time-delayed version of itself.\n\ninput: A tuple with audio samples and sample rate.\ndepth: Maximum modulation delay in seconds.\nrate: Modulation frequency in Hz.\nmix: Blend ratio between dry and flanged signals.\nfeedback: Proportion of the delayed signal fed back into the input.\n\nReturns the flanged audio as a tuple.\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.phaser","page":"Effects","title":"PulseArchitect.Effects.phaser","text":"phaser(input::Tuple{Vector{Float64}, Int};        depth::Real=0.7,        rate::Real=0.5,        stages::Int=4)\n\nCreates a phaser effect by applying multiple all-pass filters to shift the phase of the signal.\n\ninput: A tuple containing the audio samples and sample rate.\ndepth: Intensity of the phase modulation.\nrate: Modulation frequency in Hz.\nstages: Number of all-pass filter stages.\n\nReturns the phase-shifted audio as a tuple.\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.distortion","page":"Effects","title":"PulseArchitect.Effects.distortion","text":"distortion(signal; gain=20.0, threshold=0.8)\n\nApplies distortion by amplifying the signal and clipping its peaks.\n\nsignal: Input audio array.\ngain: Amplification factor before clipping.\nthreshold: Clipping threshold (normalized between 0 and 1).\n\nReturns the distorted signal.\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.bitcrusher","page":"Effects","title":"PulseArchitect.Effects.bitcrusher","text":"bitcrusher(input::Tuple{Vector{Float64}, Int};            bitdepth::Int=8,            samplerate_reduction::Int=1)\n\nReduces audio fidelity by quantizing the input signal to a lower bit depth and by reducing its effective sample rate. Assumes that the input samples are normalized in the range [-1, 1].\n\ninput: A tuple where the first element is a vector of audio samples and the second is the sample rate.\nbit_depth: The target bit depth (e.g., 8 for 8-bit resolution), which determines the number of quantization levels (2^bit_depth).\nsample_rate_reduction: The factor by which to reduce the sample rate. For instance, if set to 2, the function processes the signal in blocks of 2 samples, replacing each block with a single quantized value.\n\nThe function first computes the quantization step based on the number of levels, then processes the signal in blocks defined by sample_rate_reduction. For each block, it quantizes a representative sample (the first sample of the block) and replicates this value over the entire block.\n\nReturns a tuple (crushed_samples, sample_rate), where crushed_samples is the processed audio signal.\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.pitch_shift","page":"Effects","title":"PulseArchitect.Effects.pitch_shift","text":"pitchshift(input::Tuple{Vector{Float64}, Int}, semitones; windowsize=1024, hop_size=256)\n\nShifts the pitch of the input audio signal by the specified number of semitones using a simplified phase vocoder. The algorithm performs the following steps:\n\nComputes the short-time Fourier transform (STFT) of the input using a Hann window.\nTime-stretches the signal by a factor of 1/factor (where factor = 2^(semitones/12)) by linearly interpolating between STFT frames.\nReconstructs the time-domain signal via the inverse STFT (using overlap-add).\nResamples the resulting signal to match the original length.\n\ninput: A tuple (samples, fs) where samples is a vector of normalized audio samples (Float64, range [-1,1]) and fs is the sample rate.\nsemitones: The number of semitones to shift the pitch (positive to raise, negative to lower).\nwindow_size: Size of the FFT window (default 1024).\nhop_size: Hop size between frames (default 256).\n\nReturns a tuple (shifted_signal, fs).\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.lowpass","page":"Effects","title":"PulseArchitect.Effects.lowpass","text":"lowpass(input::Tuple{Vector{Float64}, Int}; cutoff::Real)\n\nApplies a first‐order lowpass filter to the input signal using an exponential moving average.\n\ninput: Tuple (samples, fs) where samples is a vector of normalized audio samples (range [-1, 1]) and fs is the sample rate in Hz.\ncutoff: Cutoff frequency in Hz. Frequencies above this are increasingly attenuated.\n\nReturns a tuple (filtered_samples, fs).\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.highpass","page":"Effects","title":"PulseArchitect.Effects.highpass","text":"highpass(input::Tuple{Vector{Float64}, Int}; cutoff::Real)\n\nApplies a first‐order highpass filter to the input signal using a simple recursive formula.\n\ninput: Tuple (samples, fs) where samples is a vector of normalized audio samples (range [-1, 1]) and fs is the sample rate in Hz.\ncutoff: Cutoff frequency in Hz. Frequencies below this are increasingly attenuated.\n\nReturns a tuple (filtered_samples, fs).\n\n\n\n\n\n","category":"function"},{"location":"api/effects/#PulseArchitect.Effects.bandpass","page":"Effects","title":"PulseArchitect.Effects.bandpass","text":"bandpass(input::Tuple{Vector{Float64}, Int}; lowcutoff::Real, highcutoff::Real)\n\nCreates a bandpass filter by cascading a highpass and a lowpass filter.\n\ninput: Tuple (samples, fs) where samples is a vector of normalized audio samples (range [-1, 1]) and fs is the sample rate in Hz.\nlow_cutoff: Lower cutoff frequency in Hz (highpass threshold).\nhigh_cutoff: Upper cutoff frequency in Hz (lowpass threshold).\n\nReturns a tuple (filtered_samples, fs) where only frequencies between low_cutoff and high_cutoff are preserved.\n\n\n\n\n\n","category":"function"}]
}
