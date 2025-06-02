############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    amu2kg(amu)

Converts atomic mass unit `amu` into kilo gram.

# Arguments
- `amu::Real`: Atomic mass unit value.
"""
amu2kg(amu::T) where {T<:AbstractFloat} = T(amu / AVOGADRO_CONSTANT / 1000)
amu2kg(amu::Integer) = amu / AVOGADRO_CONSTANT / 1000


"""
    eV2J(eV)

Converts energy in electron volt `eV` to joule.

# Arguments
- `eV::Real`: Energy in electron volt.
"""
eV2J(eV::T) where {T<:AbstractFloat} = T(eV * ELEMENTARY_CHARGE)
eV2J(eV::Integer) = eV * ELEMENTARY_CHARGE


"""
    eV2kJpmol(eV)

Converts energy in electron volt `eV` to kilo joule per mole.

# Arguments
- `eV::Real`: Energy in electron volt.
"""
eV2kJpmol(eV::T) where {T<:AbstractFloat} = T(eV2J(eV) * AVOGADRO_CONSTANT / 1000)
eV2kJpmol(eV::Integer) = eV2J(eV) * AVOGADRO_CONSTANT / 1000


"""
    J2eV(J)

Converts energy in joule `J` to electron volt.

# Arguments
- `J::Real`: Energy in joule.
"""
J2eV(J::T) where {T<:AbstractFloat} = T(J / ELEMENTARY_CHARGE)
J2eV(J::Integer) = J / ELEMENTARY_CHARGE


"""
    kJpmol2eV(kJpmol)

Converts energy in kilo joule per mole `kJpmol` to electron volt.

# Arguments
- `kJpmol::Real`: Energy in kilo joule per mole.
"""
kJpmol2eV(kJpmol::T) where {T<:AbstractFloat} = T(J2eV(kJpmol * 1000) / AVOGADRO_CONSTANT)
kJpmol2eV(kJpmol::Integer) = J2eV(kJpmol * 1000) / AVOGADRO_CONSTANT


"""
    limit_acos(x)

Extends the `acos` function for input outside of `[-1,1]` through clipping.

# Arguments
- `x::Real`: Input value for which to compute the arccosine.
"""
limit_acos(x::Real) = (abs(x) > 1) ? acos(sign(x)) : acos(x)


"""
    lon2LT(lon)

Converts a subsolar longitude `lon` (in rad) into local time `LT`.

# Arguments
- `lon::Real`: Subsolar longitude in radians.
"""
lon2LT(lng::T) where {T<:AbstractFloat} = T((lng + pi) / pi * 12)
lon2LT(lng::Real) = lon2LT(Float64(lng))
lon2LT(lng::BigInt) = lon2LT(BigFloat(lng))


"""
    LT2lon(LT)

Converts a local time `LT` into a subsolar longitude `lon` (in rad).

# Arguments
- `LT::Real`: Local time in hours.
"""
LT2lon(LT::T) where {T<:AbstractFloat} = T(LT / 12 * pi - pi)
LT2lon(LT::Real) = LT2lon(Float64(LT))
LT2lon(LT::BigInt) = LT2lon(BigFloat(LT))


"""
    pclamp(x, l, u)

Performs a periodic clamp of `x` between `l` and `u`, i.e., it returns `x` if it is within
the bounds, otherwise it returns `l + mod(x - l, u - l)`.

# Arguments
- `x::Real`: Input value to be clamped.
- `l::Real`: Lower bound of the clamp.
- `u::Real`: Upper bound of the clamp.
"""
function pclamp(x::T, l::T, u::T) where {T<:Real}
    if l <= x <= u; return x; end
    return T(l + mod(x - l, u - l))
end
pclamp(x::Real, l::Real, u::Real) = pclamp(promote(x, l, u)...)
# pclamp(x::Integer, l::Integer, u::Integer) = pclamp(float(x), l, u)


"""
    sgn(x)

Returns the sign of `x` with the custom definition `sgn(0) = 1`.

# Arguments
- `x::Real`: Input value for which to compute the sign.
"""
sgn(x::Real) = x != 0 ? sign(x) : typeof(x)(1)


############################################################################################
#::. EXTENSIONS
############################################################################################
"""
Extension of `Base.isapprox` function to work with tuples. Mostly internal use inside of
the test suite of the `ExESS` package.
"""
Base.isapprox(t1::Tuple, t2::Tuple; kwargs...) = isapprox([t1...], [t2...]; kwargs...)


"""
Promotion rules for `Vector{BigInt}` and `Vector{BigFloat}`. Mostly internal use inside of
the test suite of the `ExESS` package.
"""
Base.promote_rule(::Type{Vector{BigInt}}, ::Type{Vector{T}}) where {T<:Integer} = Vector{BigInt}
Base.promote_rule(::Type{Vector{BigInt}}, ::Type{Vector{T}}) where {T<:Real} = Vector{BigFloat}
Base.promote_rule(::Type{Vector{T}}, ::Type{Vector{BigInt}}) where {T<:Integer} = Vector{BigInt}
Base.promote_rule(::Type{Vector{T}}, ::Type{Vector{BigInt}}) where {T<:Real} = Vector{BigFloat}


"""
Custom extension of the `erfinv` function for `Float16` type inputs.
"""
SpecialFunctions.erfinv(x::Float16) = Float16(erfinv(Float32(x)))


############################################################################################
#::. EXPORTS
############################################################################################
export amu2kg, eV2J, eV2kJpmol, J2eV, kJpmol2eV, limit_acos, lon2LT, LT2lon, pclamp, sgn
