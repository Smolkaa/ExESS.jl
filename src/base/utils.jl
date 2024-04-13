############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] amu2kg(amu::Real)

Converts atomic mass unit `amu` into kilo gram.
"""
amu2kg(amu::Real) = amu / typeof(amu)(AVOGADRO_CONSTANT) / 1000
amu2kg(amu::Integer) = amu / AVOGADRO_CONSTANT / 1000


"""
    [1] eV2J(eV::Real)

Converts energy in electron volt `eV` to joule.
"""
eV2J(eV::Real) = eV * typeof(eV)(ELEMENTARY_CHARGE)
eV2J(eV::Integer) = eV * ELEMENTARY_CHARGE


"""
    [1] J2eV(J::Real)

Converts energy in joule `J` to electron volt.
"""
J2eV(J::Real) = J / typeof(J)(ELEMENTARY_CHARGE)
J2eV(J::Integer) = J / ELEMENTARY_CHARGE


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
    [1] sgn(x::Real)

Returns the sign of `x` with the custom definition `sgn(0) = 1`.
"""
sgn(x::Real) = x != 0 ? sign(x) : typeof(x)(1) 


############################################################################################
#::. EXTENSIONS
############################################################################################
"""
    [1] Base.isapprox(t1::Tuple, t2::Tuple; kwargs...)

Extension of `Base.isapprox` function to work with tuples. Mostly internal use inside of
the test suite of the `ExESS` package.
"""
Base.isapprox(t1::Tuple, t2::Tuple; kwargs...) = isapprox([t1...], [t2...]; kwargs...)


"""
    [1] Base.promote_rule(::Type{Vector{BigInt}}, ::Type{Vector{T}}) where {T<:Integer}
    [2] Base.promote_rule(::Type{Vector{BigInt}}, ::Type{Vector{T}}) where {T<:Real}
    [3] Base.promote_rule(::Type{Vector{T}}, ::Type{Vector{BigInt}}) where {T<:Integer}
    [4] Base.promote_rule(::Type{Vector{T}}, ::Type{Vector{BigInt}}) where {T<:Real}

Promotion rules for `Vector{BigInt}` and `Vector{BigFloat}`. Mostly internal use inside of
the test suite of the `ExESS` package.
"""
Base.promote_rule(::Type{Vector{BigInt}}, ::Type{Vector{T}}) where {T<:Integer} = Vector{BigInt}
Base.promote_rule(::Type{Vector{BigInt}}, ::Type{Vector{T}}) where {T<:Real} = Vector{BigFloat}
Base.promote_rule(::Type{Vector{T}}, ::Type{Vector{BigInt}}) where {T<:Integer} = Vector{BigInt}
Base.promote_rule(::Type{Vector{T}}, ::Type{Vector{BigInt}}) where {T<:Real} = Vector{BigFloat}


############################################################################################
#::. EXPORTS
############################################################################################
export amu2kg, eV2J, J2eV, limit_acos, lng2LT, LT2lng, sgn