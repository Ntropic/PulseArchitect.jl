module Wave_Generators
    
export square_wave, saw_wave, triangle_wave, pulse_wave, soft_square, fm_wave, wavy_wave


"""
    square_wave(phase::Real)

Returns 1.0 for phase < π, else -1.0.
"""
square_wave(phase::Real)   = ifelse(phase < π, 1.0, -1.0)

"""
    saw_wave(phase::Real)

Returns a sawtooth wave in [-1, 1] from phase.
"""
saw_wave(phase::Real)      = 2*(phase/(2π)) - 1.0

"""
    triangle_wave(phase::Real)

Returns a triangle wave in [-1, 1] from phase.
"""
triangle_wave(phase::Real) = phase < π ? 2*phase/π - 1.0 : 1.0 - 2*(phase-π)/π

"""
    pulse_wave(phase::Real; duty=0.3)

Returns a pulse wave with given duty cycle.
"""
pulse_wave(phase::Real; duty=0.3) = phase < 2π*duty ? 1.0 : -1.0

"""
    soft_square(phase::Real)

Returns a soft square wave via tanh(5*sin(phase)).
"""
soft_square(phase::Real)   = tanh(5*sin(phase))

"""
    fm_wave(phase::Real)

Returns a frequency-modulated wave: sin(phase + 0.5*sin(3*phase)).
"""
fm_wave(phase::Real)       = sin(phase + 0.5*sin(3*phase))

"""
    wavy_wave(phase::Real)

Returns a composite wave: sin(phase) + 0.5*sin(2*phase + 0.5).
"""
wavy_wave(phase::Real)     = sin(phase) + 0.5*sin(2*phase + 0.5)


end