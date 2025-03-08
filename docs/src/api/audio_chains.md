## Audio Chains and Splitters

`AudioChain`s and `Splitter`s allow us to define sequential and parallel effect processing pipelines for audio signals. They are particularly useful in scenarios where multiple effects need to be applied or when different parts of the signal need to be processed independently.

### Audio Chain
`AudioChain`s are defined as
```@docs
AudioChain
```

### Splitter
`Splitter` allows us to split an audio signal into multiple channels and process each channel separately. It is useful when we need to apply different effects to different parts of the signal.
```@docs
Splitter
```

