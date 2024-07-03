############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] amu2kg(amu::Real)

Converts atomic mass unit `amu` into kilo gram.
"""
amu2kg(amu::T) where {T<:AbstractFloat} = T(amu / AVOGADRO_CONSTANT / 1000)
amu2kg(amu::Integer) = amu / AVOGADRO_CONSTANT / 1000


"""
    [1] eV2J(eV::Real)

Converts energy in electron volt `eV` to joule.
"""
eV2J(eV::T) where {T<:AbstractFloat} = T(eV * ELEMENTARY_CHARGE)
eV2J(eV::Integer) = eV * ELEMENTARY_CHARGE


"""
    [1] eV2kJpmol(eV::Real)

Converts energy in electron volt `eV` to kilo joule per mole.
"""
eV2kJpmol(eV::T) where {T<:AbstractFloat} = T(eV2J(eV) * AVOGADRO_CONSTANT / 1000)
eV2kJpmol(eV::Integer) = eV2J(eV) * AVOGADRO_CONSTANT / 1000


"""
    [1] J2eV(J::Real)

Converts energy in joule `J` to electron volt.
"""
J2eV(J::T) where {T<:AbstractFloat} = T(J / ELEMENTARY_CHARGE)
J2eV(J::Integer) = J / ELEMENTARY_CHARGE


"""
    [1] kJpmol2eV(kJpmol::Real)

Converts energy in kilo joule per mole `kJpmol` to electron volt.
"""
kJpmol2eV(kJpmol::T) where {T<:AbstractFloat} = T(J2eV(kJpmol * 1000) / AVOGADRO_CONSTANT)
kJpmol2eV(kJpmol::Integer) = J2eV(kJpmol * 1000) / AVOGADRO_CONSTANT


"""
    [1] limit_acos(x::Real)

Extends the `acos` function for input outside of `[-1,1]` through clipping.
"""
limit_acos(x::Real) = (abs(x) > 1) ? acos(sign(x)) : acos(x)


"""
    [1] lng2LT(lng::Real)

Converts a subsolar longitude `lng` into local time `LT`.
"""
lng2LT(lng::T) where {T<:AbstractFloat} = T((lng + pi) / pi * 12)
lng2LT(lng::Real) = lng2LT(Float64(lng))
lng2LT(lng::BigInt) = lng2LT(BigFloat(lng))


"""
    [1] LT2lng(LT::Real)

Converts a local time `LT` into a subsolar longitude `lng`.
"""
LT2lng(LT::T) where {T<:AbstractFloat} = T(LT / 12 * pi - pi)
LT2lng(LT::Real) = LT2lng(Float64(LT))
LT2lng(LT::BigInt) = LT2lng(BigFloat(LT))


"""

"""
function pclamp(x::T, l::T, u::T) where {T<:Real}
    return T(l + mod(x - l, u - l))
end
pclamp(x::Real, l::Real, u::Real) = pclamp(promote(x, l, u)...)
# pclamp(x::Integer, l::Integer, u::Integer) = pclamp(float(x), l, u)


"""
    [1] sgn(x::Real)

Returns the sign of `x` with the custom definition `sgn(0) = 1`.
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
export amu2kg, eV2J, eV2kJpmol, J2eV, kJpmol2eV, limit_acos, lng2LT, LT2lng, sgn
