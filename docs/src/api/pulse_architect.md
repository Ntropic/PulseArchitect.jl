## Pulse Architect

In order to apply effects or chains of effects to an audio signal, we first have to create an audio signal.
One way to do so is the generate_melody function. 
To use it, specify a melody as a Vector of Tuples. Each Tuple contains three elements: the note name (e.g., "C4"), its duration in beats, and its amplitude. Generate melody can use this to generate an initial audio signal, using optional parameters, such as a waveform generator, the sample rate, portamento time (smoothly transitioning between notes), a decay function (for plucked note sound) and more.
```@docs
generate_melody
```
Every audio object is a tuple of the form `(audio_signal, fs)` where `audio_signal` is an array of floats representing the audio waveform and `fs` is the sample rate in Hz. This allows us to manipulate the audio signal using various functions provided by PulseArchitect.

We can play an audio Tuple via the `play` function. This does not block the execution of other commands, as it is started in a separate thread. 
Audio Signals can be saved via the `save` function.
```@docs
play
save
```








 
