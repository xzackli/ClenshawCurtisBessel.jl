
# define operations on ArbFloat
HypergeometricFunctions.logabsgamma(x::ArbNumerics.ArbFloat) =
	ArbNumerics.lgamma(abs(x)), sign(gamma(x))

lgamma(x::ArbNumerics.ArbFloat) = ArbNumerics.lgamma(x)
lgamma(x::DoubleFloats.DoubleFloat) = DoubleFloats.lgamma(x)
lgamma(x::Number) = SpecialFunctions.lgamma(x)

J(ν::ArbNumerics.ArbFloat, x) = ArbNumerics.besselj(ν, x)
J(ν, x::ArbNumerics.ArbFloat) = ArbNumerics.besselj(ν, x)
J(ν::ArbNumerics.ArbFloat, x::ArbNumerics.ArbFloat) = ArbNumerics.besselj(ν, x)
J(ν::Float64, x::Float64) = SpecialFunctions.besselj(ν, x)

# doublefloat is missing half-integer bessel
J(ν::DoubleFloats.DoubleFloat, x::DoubleFloats.DoubleFloat) =
    DoubleFloats.Double64(
        (ArbNumerics.besselj(ArbNumerics.ArbFloat(ν), ArbNumerics.ArbFloat(x))))
J(ν::DoubleFloats.DoubleFloat, x::Complex{DoubleFloats.DoubleFloat}) =
    DoubleFloats.Double64(
        (ArbNumerics.besselj(ArbNumerics.ArbFloat(ν), ArbNumerics.ArbFloat(x))))
