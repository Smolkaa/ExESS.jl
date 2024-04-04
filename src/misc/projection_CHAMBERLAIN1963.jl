############################################################################################
#::. FUNCTIONS
# TODO: Something does not add up here. Very close projections are exceeding 100%. Write
#       a test to check the implementation. Check publication again.
# 		UPDATE (02/04/2024): I checked the publication, the partition functions seem to be
#							 correct. The problem is that the close projections are not at
#							 at 1.0, so there is a slight jump. This aparently is carried
#							 over, as two separate projections are also not equal to one big
#							 projection, no matter the distances.
# 
# Note:	The potential energies are converted to Float64 to prevent them from becoming a 
#		`BigFloat` type, as this generally causes problems with the gamma functions.
############################################################################################
"""
    [1] projection_CHAMBERLAIN1963([S::Type{<:AbstractFloat},] r1::Real, r2::Real, M::Real, 
								   m::Real, T::Real; kwargs...)

Calculates the projection fraction for a particles of mass `m` at radial distance `r2`,
assuming a known value given at radial distance `r1` of the central body with mass `M`. The 
entire projection is based on an isothermal derivation assuming constant temperature `T`.

Assumes hydrostatic equilibrium and perfect-gas law with isotropic gas pressure, resulting
in the generalized form of the isothermal barometric law. The calculation is based on the 
equation (14) in Chamberlain (1963).

**Key-Word Arguments**

| Field    | Default Value  | Unit       | Description                          |
|:-------- | --------------:|:---------- |:------------------------------------ |
| `bal`    | `true`         | ()         | include ballistic particles          |
| `sat`    | `true`         | ()         | include satellite particles          |
| `esc`    | `true`         | ()         | include escaping particles           |

**Rerences**

- Chamberlain, 1963, "Planetary coronae and atmospheric evaporation"
- Cook et al. 2013, "New upper limits on numerous atmospheric species in the native lunar
  atmosphere"
"""
function projection_CHAMBERLAIN1963(r1::S, r2::S, M::S, m::S, T::S; bal::Bool=true, 
								    sat::Bool=true, esc::Bool=true) where {S<:AbstractFloat}
    if r1 < r2
		pot_energy_1 = Float64(GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * r1))
		pot_energy_2 = Float64(GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * r2))
		
		# pre-calculations
		L = pot_energy_2^2 / (pot_energy_1 + pot_energy_2)
		g =  gamma(1.5)
		gi = gamma_inc(1.5, pot_energy_2)[1] * g
		giL = gamma_inc(1.5, pot_energy_2 - L)[1] * g
		f = sqrt(pot_energy_1^2 - pot_energy_2^2) / pot_energy_1
	
		# partition functions (ballistic, satellite, escape)
		par = 0
		if bal; par += 2/sqrt(pi) * (gi - f*exp(-L)*giL); end
		if sat; par += 2/sqrt(pi) * (f*exp(-L)*giL); end
		if esc; par += 1/sqrt(pi) * (g - gi - f*exp(-L)*(g - giL)); end
	
		# return projection
		return S(exp(-(pot_energy_1 - pot_energy_2)) * par)

	elseif r1 > r2
		pot_energy_1 = Float64(GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * r2))
		pot_energy_2 = Float64(GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * r1))
		
		# pre-calculations
		L = pot_energy_2^2 / (pot_energy_1 + pot_energy_2)
		g =  gamma(1.5)
		gi = gamma_inc(1.5, pot_energy_2)[1] * g
		giL = gamma_inc(1.5, pot_energy_2 - L)[1] * g
		f = sqrt(pot_energy_1^2 - pot_energy_2^2) / pot_energy_1
	
		# partition functions (ballistic, satellite, escape)
		par = 0
		if bal; par += 2/sqrt(pi) * (gi - f*exp(-L)*giL); end
		if sat; par += 2/sqrt(pi) * (f*exp(-L)*giL); end
		if esc; par += 1/sqrt(pi) * (g - gi - f*exp(-L)*(g - giL)); end
	
		# return projection
		return S(exp(pot_energy_1 - pot_energy_2) / par)
	end
	return S(1)
end
function projection_CHAMBERLAIN1963(r1::Real, r2::Real, M::Real, m::Real, T::Real; kwargs...)
	return projection_CHAMBERLAIN1963(promote(r1, r2, M, m, T)...; kwargs...)
end
function projection_CHAMBERLAIN1963(r1::Integer, r2::Integer, M::Integer, m::Integer, 
									T::Integer; kwargs...)
	return projection_CHAMBERLAIN1963(float(r1),float(r2),float(M),float(m),float(T);kwargs...)
end
function projection_CHAMBERLAIN1963(S::Type{<:AbstractFloat}, args...; kwargs...)
	return S(projection_CHAMBERLAIN1963(args...; kwargs...))
end



############################################################################################
#::. EXPORTS
############################################################################################
export 
    projection_CHAMBERLAIN1963