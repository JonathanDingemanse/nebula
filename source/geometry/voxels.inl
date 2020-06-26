namespace nbl { namespace geometry {

template<bool gpu_flag>
CPU voxels<gpu_flag> voxels<gpu_flag>::create(std::vector<triangle> const & triangles)
{
	const real VOXEL_SIZE = 0.3; // voxel size in nanometers (0.27 nm is appr. one atom of Si)

	const int SIZE_X = 201; // horizontal size in the x direction in voxels
	const int SIZE_Y = 201; // horizontal size in the y direction in voxels
	const int SIZE_Z = 700; // vertical size in voxels

	const int SAMPLE_HEIGHT = 300; // height of the sample (length between the sample and the top of the simulation domain) in voxels

	// Hier moet ergens een grid gemaakt worden, 
		// een voor het materiaal, 
		// een voor the timestamp van de electronen, 
		// een voor de energy van het deposition electron,
		// en een voor de dz van dat electron
	
	mat_grid.reshape(SIZE_X*SIZE_Y*SIZE_Z, 0) // set the shape 

	// set the initial geometry
	for(i = 0, i < SIZE_X, i++){ 
		for(j = 0, j < SIZE_Y, j++){
			for(k = 0, k < SAMPLE_HEIGHT, k++){
				mat_grid[i + j*SIZE_X + k*SIZE_X*SIZE_Y] = -123
			}
		}
	}

	AABB_min = vec3{ 0, 0, 0 };
	AABB_max = vec3{ SIZE_X*VOXEL_SIZE, SIZE_Y*VOXEL_SIZE, SIZE_Z*VOXEL_SIZE };

	// TODO: error message
	/*if (triangles.empty())
		throw std::runtime_error("No triangles provided!");

	vec3 AABB_min = triangles[0].AABB_min();
	vec3 AABB_max = triangles[0].AABB_max();

	for (const triangle t : triangles)
	{
		const vec3 tri_min = t.AABB_min();
		const vec3 tri_max = t.AABB_max();

		AABB_min =
		{
			std::min(AABB_min.x, tri_min.x),
			std::min(AABB_min.y, tri_min.y),
			std::min(AABB_min.z, tri_min.z)
		};
		AABB_max =
		{
			std::max(AABB_max.x, tri_max.x),
			std::max(AABB_max.y, tri_max.y),
			std::max(AABB_max.z, tri_max.z)
		};
	}

	AABB_min -= vec3{ 1, 1, 1 };
	AABB_max += vec3{ 1, 1, 1 };
	*/

	return detail::voxels_factory<gpu_flag>::create(triangles, AABB_min, AABB_max);
}

template<bool gpu_flag>
CPU void voxels<gpu_flag>::destroy(voxels<gpu_flag> & geometry)
{
	detail::voxels_factory<gpu_flag>::destroy(geometry);
}

template<bool gpu_flag>
PHYSICS bool voxels<gpu_flag>::in_domain(vec3 pos)
{
	return ((pos.x > _AABB_min.x) && (pos.x < _AABB_max.x)
		&& (pos.y > _AABB_min.y) && (pos.y < _AABB_max.y)
		&& (pos.z > _AABB_min.z) && (pos.z < _AABB_max.z));
}

template<bool gpu_flag>
PHYSICS intersect_event voxels<gpu_flag>::propagate(vec3 start, vec3 direction, real distance,
	triangle const * ignore_triangle, int ignore_material) const
{
	intersect_event evt { distance, nullptr };

	k = (int) (start.x / voxel_size);
	l = (int) (start.y / voxel_size);
	m = (int) (start.z / voxel_size);
	
	start_mat = mat_grid[k + l*SIZE_X + m*SIZE_X*SIZE_Y];

	real dx = direction.x / voxel_size; // create vars for the direction elements
	real dy = direction.y / voxel_size;
	real dz = direction.z / voxel_size;

	vec3 dr = {dx, dy, dz};
	vec3 delta_S = {0, 0, 0};

	real delta_s_min = distance / voxel_size;
 
	if(dx > 0){
		delta_x = std::ceil(dx) - dx;
	}
	else{
		delta_x = dx - std::floor(dx);
	}
	delta_S.x = delta_x/dx;
	
	if(delta_S_min > delta_S.x){
		delta_S_min = delta_S.x;
	}

	if(dy > 0){
		delta_y = std::ceil(dy) - dy;
	}
	else{
		delta_y = dy - std::floor(dy);
	}
	delta_S.y = delta_y/dy;

	if(dz > 0){
		delta_z = std::ceil(dz) - dz;
	}
	else{
		delta_z = dz - std::floor(dz);
	}
	delta_S.z = delta_z/dz;
	
	while(distance > delta_s_min){

		for(int i ; i < 3 : i++){

			real delta_S_i = delta_S[i]

			if(delta_S_i < delta_s_min){
				delta_s_min = delta_S_i;

				vec3 new_pos = start + delta_s_min*dr;

				int k;
				int l;
				int m;

				switch (i)
				{
				case 0: // intersection with x-plane
					k = (int) (new_pos.x + 0.1);
					l = (int) (new_pos.y);
					m = (int) (new_pos.z);
					break;
				
				case 1: // intersection with y-plane
					k = (int) (new_pos.x);
					l = (int) (new_pos.y + 0.1);
					m = (int) (new_pos.z);
					break;

				default: // intersection with z-plane
					k = (int) (new_pos.x);
					l = (int) (new_pos.y);
					m = (int) (new_pos.z + 0.1);
					break;
				}

				if(mat_grid[ k + m*SIZE_X + l*SIZE_X*SIZE_Y] != start_mat){
					evt.isect_distance = delta_s_min * voxel_size;
					switch (i)
					{
					case 0:
						if(dx > 0){
							evt.isect_triangle = 1;
						}
						else{
							evt.isect_triangle = 2;
						}
						break;

					case 1:
						if(dy > 0){
							evt.isect_triangle = 3;
						}
						else{
							evt.isect_triangle = 4;
						}
						break;	
					
					default:
						if(dz > 0){
							evt.isect_triangle = 5;
						}
						else{
							evt.isect_triangle = 6;
						}
						break;
					}
				}
				
			}
		}
	}

	/*for (triangle_index_t i = 0; i < _N; ++i)
	{
		if (_triangles + i == ignore_triangle)
			continue;

		// retrieve triangle from global memory
		const triangle this_triangle = _triangles[i];

		int mat_idx_out;
		if (dot_product(this_triangle.get_normal(), direction) < 0)
			mat_idx_out = this_triangle.material_in;
		else
			mat_idx_out = this_triangle.material_out;

		// if the outgoing material is the same as current, nothing happens
		// if the triangle represents a detector which can't see the current
		// particle, nothing happens
		if (mat_idx_out == ignore_material)
			continue;

		const real t = this_triangle.intersect_ray(start, direction);
		if (t > 0 && t < evt.isect_distance)
		{
			evt.isect_distance = t;
			evt.isect_triangle = _triangles + i;
		}
	}*/
	return evt;
}

template<bool gpu_flag>
PHYSICS real voxels<gpu_flag>::get_max_extent() const
{
	return _max_extent;
}

template<bool gpu_flag>
inline PHYSICS vec3 voxels<gpu_flag>::AABB_min() const
{
	return _AABB_min;
}
template<bool gpu_flag>
inline PHYSICS vec3 voxels<gpu_flag>::AABB_max() const
{
	return _AABB_max;
}

template<bool gpu_flag>
CPU void voxels<gpu_flag>::set_AABB(vec3 min, vec3 max)
{
	_AABB_min = min;
	_AABB_max = max;

	const vec3 m = max - min;
	_max_extent = magnitude(m);
}


namespace detail
{
	template<>
	struct voxels_factory<false>
	{
		inline static CPU voxels<false> create(std::vector<triangle> triangles, vec3 AABB_min, vec3 AABB_max)
		{
			using voxels_t = voxels<false>;
			using triangle_index_t = voxels_t::triangle_index_t;

			voxels_t geometry;

			/*if (triangles.size() > std::numeric_limits<triangle_index_t>::max())
				throw std::runtime_error("Too many triangles in geometry");
			geometry._N = static_cast<triangle_index_t>(triangles.size());

			geometry._triangles = reinterpret_cast<triangle*>(malloc(geometry._N * sizeof(triangle)));
			for (triangle_index_t i = 0; i < triangles.size(); ++i)
			{
				geometry._triangles[i] = triangles[i];
			}*/

			geometry.set_AABB(AABB_min, AABB_max);

			return geometry;
		}

		inline static CPU void free(voxels<false> & geometry)
		{
			::free(geometry._triangles);

			geometry._triangles = nullptr;
			geometry._N = 0;
		}
	};

#if CUDA_COMPILER_AVAILABLE
	template<>
	struct voxels_factory<true>
	{
		inline static CPU voxels<true> create(std::vector<triangle> triangles, vec3 AABB_min, vec3 AABB_max)
		{
			using voxels_t = voxels<true>;
			using triangle_index_t = voxels_t::triangle_index_t;

			voxels_t geometry;

			if (triangles.size() > std::numeric_limits<triangle_index_t>::max())
				throw std::runtime_error("Too many triangles in geometry");
			geometry._N = static_cast<triangle_index_t>(triangles.size());

			// Copy triangle data to device
			cuda::cuda_new<triangle>(&geometry._triangles, geometry._N);
			cuda::cuda_mem_scope<triangle>(geometry._triangles, geometry._N,
				[&triangles](triangle* device)
			{
				for (triangle_index_t i = 0; i < triangles.size(); ++i)
					device[i] = triangles[i];
			});

			geometry.set_AABB(AABB_min, AABB_max);

			return geometry;
		}

		inline static CPU void free(voxels<true> & geometry)
		{
			cudaFree(geometry._triangles);

			geometry._triangles = nullptr;
			geometry._N = 0;
		}
	};
#endif // CUDA_COMPILER_AVAILABLE
} // namespace detail

}} // namespace nbl::geometry