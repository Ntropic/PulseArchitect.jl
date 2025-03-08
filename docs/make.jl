using PulseArchitect
using Documenter

makedocs(
    modules = [PulseArchitect],
    sitename = "PulseArchitect.jl",
    pages = [
        "Home" => "index.md",
        "API Reference" => [
            "Pulse Architect" => "api/pulse_architect.md",
            "Effects" => "api/effects.md",        
            "Audio Chains" => "api/audio_chains.md",      
            "Wave Generators" => "api/wave_generators.md",  
        ]
    ]
)
