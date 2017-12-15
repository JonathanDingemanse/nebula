#ifndef __INTERSECT_THOMAS_H_
#define __INTERSECT_THOMAS_H_

// TODO: store barrier here.

template<
	bool opt_quantum_transmission = true,
	bool opt_interface_refraction = true,
	bool opt_interface_absorption = false>
struct intersect_thomas
{
	template<typename particle_manager, typename material_manager, bool gpu_flag>
	PHYSICS void execute(material_manager& material_mgr,
		particle_manager& particle_mgr, typename particle_manager::particle_index_t particle_idx,
		nbl::util::random_generator<gpu_flag>& rng) const
	{
		using material_index_t = typename material_manager::material_index_t;

		// Get particle data from particle manager
		auto this_particle = particle_mgr[particle_idx];
		const triangle this_triangle = *(particle_mgr.get_last_triangle(particle_idx));

		// Get particle direction, normal of the intersected triangle
		auto normalised_dir = normalised(this_particle.dir);
		auto last_triangle_normal = normalised(this_triangle.get_normal());

		// Get angle between direction of motion and triangle
		const real cos_theta = dot_product(last_triangle_normal, normalised_dir);


		// determine the material index in the following way:
		//   material_idx_in represents the current material
		//   material_idx_out represents the material when passing through the interface
		material_index_t material_idx_in, material_idx_out;
		if (cos_theta > 0)
		{
			material_idx_in = this_triangle.material_in;
			material_idx_out = this_triangle.material_out;
		}
		else
		{
			material_idx_in = this_triangle.material_out;
			material_idx_out = this_triangle.material_in;
		}


		// manage special cases for electron detection, electron mirrors and terminators.
		//   DETECTOR always detects.
		//   DETECTOR_LT/GE50 detect under certain circumstances.
		//     If not detected, they pass through as if no intersection event has taken place.
		switch (material_idx_out) {
		case material_manager::DETECTOR:
			particle_mgr.detect(particle_idx);
			return;
		case material_manager::DETECTOR_LT50:
			if (this_particle.kin_energy < 50)
				particle_mgr.detect(particle_idx);
			return;
		case material_manager::DETECTOR_GE50:
			if (this_particle.kin_energy >= 50)
				particle_mgr.detect(particle_idx);
			return;
		case material_manager::TERMINATOR:
			particle_mgr.terminate(particle_idx);
			return;
		case material_manager::MIRROR:
			this_particle.dir = normalised_dir - 2*last_triangle_normal*cos_theta;
			particle_mgr[particle_idx] = this_particle;
			return;
		default:
			break;
		}


		// determine the change in energy `dU` (in eV) when passing through the interface
		// see thesis T.V. Eq. 3.136
		real dU = 0;
		if (material_mgr.is_physical(material_idx_in)) {
			dU -= material_mgr[material_idx_in].barrier;
		}
		if (material_mgr.is_physical(material_idx_out)) {
			dU += material_mgr[material_idx_out].barrier;
		}

		// determine transmission probability (only if energy suffices)
		// see thesis T.V. Eq. 3.145
		if (this_particle.kin_energy*cos_theta*cos_theta + dU > 0)
		{
			const real s = sqrtr(1 + dU / (this_particle.kin_energy*cos_theta*cos_theta));
			const real T = (opt_quantum_transmission
				? 4 * s / ((1 + s)*(1 + s))
				: 1);
			if (rng.unit() < T)
			{
				if (opt_interface_refraction)
				{
					// determine the angle of refraction
					// see thesis T.V. Eq. 3.139
					this_particle.dir = (normalised_dir - last_triangle_normal*cos_theta)
						+ s * last_triangle_normal * cos_theta;
				}

				// if there is transmission, then adjust the kinetic energy,
				// update the current material index and EXIT.
				this_particle.kin_energy += dU;

				particle_mgr.set_material_index(particle_idx, material_idx_out);
				particle_mgr[particle_idx] = this_particle;
				return;
			}
		}

		// surface absorption? (this is in accordance with Kieft & Bosch code)
		// note that the default behaviour has this feature disabled
		if (opt_interface_absorption)
		{
			if (dU < 0 && rng.unit() < expr(1 + 0.5_r*this_particle.kin_energy / dU))
			{
				particle_mgr.terminate(particle_idx);
				return;
			}
		}

		// the only remaining case is total internal reflection
		this_particle.dir = normalised_dir - 2 * last_triangle_normal*cos_theta;
		particle_mgr[particle_idx] = this_particle;
	}
};

#endif // __INTERSECT_THOMAS_H_
