module Tsys

export model_flux
export apparent_flux
export apparent_temperature
export tcal_onoff
export tsys_onoff

"""
Boltzmann constant (Joules/Kelvin)
"""
Kb = 1.380469e-23

"""
Janskys per Joule
"""
JY_PER_J = 1e26

"""
Temperature of the *cosmic microwave background* (CMB) in Kelvin
"""
Tcmb = 2.728

"""
Default calibrator flux model frequency for `ν1` in `flux`.
"""
ν1 = Ref{Float64}(1e6)

"""
    model_flux(model, hz=ν1[], ν1=ν1[])

Compute the spectral flux density at frequency `hz` (in Hz) given flux model
coefficients `model` based on frequency units of `ν1` Hz.  `hz` and `ν1` default
to `Tsys.ν1[]` if not given.
"""
function model_flux(model, hz=ν1[], ν1=ν1[])
    10^evalpoly(log10(hz/ν1), model)
end

"""
    apparent_flux(k, diameter; η=1.0, τ=0.0, airmass=1.0)

Return the apparent spectral flux density in Janskys corresponding to black body
temperature `k` in Janskys and antenna diameter `diameter` in meters with
aperture efficiency of `η`, which defaults to 1.0.  Atmospheric opacity `τ`
defaults to 0.0 (totally transparent) and `airmass` defaults to 1.0 (airmass at
zenith).
"""
function apparent_flux(k, diameter; η=1.0, τ=0.0, airmass=1.0)
	a = η * π * diameter^2 / 4
    (2Kb * JY_PER_J * k) / (a * exp(-τ*airmass))
end

"""
    apparent_temperature(s, diameter; η=1.0, τ=0.0, airmass=1.0)

Return the apparent black body temperature in Kelvin corresponding to spectral
flux density `s` in Janskys and antenna diameter `diameter` in meters with
aperture efficiency of `η`, which defaults to 1.0.  Atmospheric opacity `τ`
defaults to 0.0 (totally transparent) and `airmass` defaults to 1.0 (airmass at
zenith).
"""
function apparent_temperature(s, diameter; η=1.0, τ=0.0, airmass=1.0)
	a = η * π * diameter^2 / 4
	s * a * exp(-τ*airmass) / (2Kb * JY_PER_J)
end

"""
    tcal_onoff(pon, poff, tsys; tsky=Tcmb, clip=true)

Calculates calibrator apparent temperature `Tcal`.  The `pon` and `poff`
parameters are the power measurements when pointing ON and OFF the calibrator
using a radio telescope with system temperature `tsys` (in Kelvin).  `tsky` is
the apparent temperature in Kelvin of the "blank sky" at the OFF position.  If
`clip` is true, any computed value less than 0.0 will be clamped to 0.0.

The calculation is derived by observing that:

    pon  = G(Tsky + Tsys + Tcal)
    poff = G(Tsky + Tsys)

Then:

     pon - poff        Tcal
    ------------ = -------------
         poff       Tsky + Tsys

Leads to:

            pon - poff
    Tcal = ------------ * (Tsky + Tsys)
               poff     
"""
function tcal_onoff(pon, poff, tsys; tsky=Tcmb, clip=true)
    lim = clip ? 0.0 : -Inf
    y = (pon - poff) / poff
    clamp(y * (tsky + tsys), lim, Inf)
end

"""
    tsys_onoff(pon, poff, tcal; tsky=Tcmb, clip=true)

Calculates system temperature `Tsys`.  The `pon` and `poff` parameters are the
power measurements when pointing ON and OFF the calibrator with apparent
temperature `tcal` (in Kelvin).  `tsky` is the apparent temperature in Kelvin of
the "blank sky" at the OFF position.  If `clip` is true, any computed value less
than 0.0 will be clamped to 0.0.

The calculation is derived by observing that:

    pon  = G(Tsky + Tsys + Tcal)
    poff = G(Tsky + Tsys)

Then:

     pon - poff        Tcal
    ------------ = -------------
         poff       Tsky + Tsys

Leads to:

                     poff
    Tsys = Tcal * ------------ - Tsky
                   pon - poff
"""
function tsys_onoff(pon, poff, tcal; tsky=Tcmb, clip=true)
    lim = clip ? 0.0 : -Inf
    y = (pon - poff) / poff
    clamp(tcal / y - tsky, lim, Inf)
end

"""
Returns Tsys by computing Tcal from `scal` and passing it to `tsys_onoff`.
Keyword arguments are the same as for `apparent_temperature` and the three
positional parameter method of `tsys_onoff`.
"""
function tsys_onoff(pon, poff, scal, diameter;
              η=1.0, τ=0.0, airmass=1.0, tsky=Tcmb, clip=true)
    tcal = apparent_temperature(scal, diameter; η=η, τ=τ, airmass=airmass)
    tsys_onoff(pon, poff, tcal; tsky=tsky, clip=clip)
end

end # module Tsys
