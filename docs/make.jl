using PulseArchitect
using Documenter

makedocs(
    sitename = "PulseArchitect.jl",
    modules = [PulseArchitect],
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
