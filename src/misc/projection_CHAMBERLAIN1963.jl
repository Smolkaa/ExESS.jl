############################################################################################
#::. FUNCTIONS
# TODO: Something does not add up here. Very close projections are exceeding 100%. Write
#       a test to check the implementation. Check publication again.
# 		UPDATE (02/04/2024): I checked the publication, the partition functions seem to be
#							 correct. The problem is that the close projections are not at
#							 at 1.0, so there is a slight jump. This aparently is carried
#							 over, as two separate projections are also not equal to one big
#							 projection, no matter the distances.
# TODO: Type stability is missing here.
############################################################################################
"""
    [1] projection_CHAMBERLAIN1963(r1::Real, r2::Real, T::Real, m::Real; 
                                   M::Real=LUNAR_MASS, bal::Bool=true, sat::Bool=true, 
                                   esc::Bool=true)

Calculates the projection fraction for a particles of mass `m` at radial distance `r2`,
assuming a known value given at radial distance `r1`. The entire projection is based on an
isothermal derivation assuming constant temperature `T`.

Assumes hydrostatic equilibrium and perfect-gas law with isotropic gas pressure, resulting
in the generalized form of the isothermal barometric law. The calculation is based on the 
equation (14) in Chamberlain (1963).


**Key-Word Arguments**

| Field    | Default Value  | Unit       | Description                          |
|:-------- | --------------:|:---------- |:------------------------------------ |
| `M`      | `LUNAR_MASS`   | (kg)       | mass of central object               |
| `bal`    | `true`         | ()         | include ballistic particles          |
| `sat`    | `true`         | ()         | include satellite particles          |
| `esc`    | `true`         | ()         | include escaping particles           |


**Rerences**

- Chamberlain, 1963, "Planetary coronae and atmospheric evaporation"
- Cook et al. 2013, "New upper limits on numerous atmospheric species in the native lunar
  atmosphere"
"""
function projection_CHAMBERLAIN1963(r1::Real, r2::Real, T::Real, m::Real; 
                                    M::Real=LUNAR_MASS, bal::Bool=true, sat::Bool=true, 
                                    esc::Bool=true)
    if r1 < r2
		pot_energy_1 = GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * r1)
		pot_energy_2 = GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * r2)
		
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
		# return par
		return exp(-(pot_energy_1 - pot_energy_2)) * par

	elseif r1 > r2
		pot_energy_1 = GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * r2)
		pot_energy_2 = GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * r1)
		
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
		return exp(pot_energy_1 - pot_energy_2) / par
	end

	return 1
end


############################################################################################
#::. EXPORTS
############################################################################################
export 
    projection_CHAMBERLAIN1963