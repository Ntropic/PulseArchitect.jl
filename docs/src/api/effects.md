## Effects 

The pulse Architect provides you with a wide range of audio processing effects that can be applied to your audio signals. As input they take a Tuple of the wave Vector and the sample rate, as output they return a Tuple of the modified wave vector and the new sample rate. 

A audio effect can be directly called for example via `tremolo((wave, fs); strength=0.3, rate=4.0)` or indirectly, by adding it to an `AudioChain` via 
```julia 
chain = AudioChain(tremolo; strength=0.3, rate=4.0)
# or by pushing it into an existing chain
push!(chain, tremolo; strength=0.5, rate=2.0)
```
Below is a list of the available effects:
```@docs
tremolo
delay_reverb
compressor
vibrato
chorus
flanger
phaser
distortion
bitcrusher
pitch_shift
lowpass
highpass
bandpass
```
