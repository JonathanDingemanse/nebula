#include "voxels.h"
namespace nbl { namespace geometry {

template<bool gpu_flag>
inline voxels<gpu_flag>::voxels(real voxel_size, vec3 shape, std::vector<int> initial_geometry)
{
	_voxel_size = voxel_size;
	_AABB_min = vec3{ 0, 0, 0 };
	_AABB_max = vec3{ shape.x * voxel_size, shape.y * voxel_size, shape.z * voxel_size };

	
	_size_x = (int)shape.x;
	_size_y = (int)shape.y;;
	_size_z = (int)shape.z;
	
	const vec3 m = _AABB_max - _AABB_min;
	_max_extent = magnitude(m);

	_mat_grid.resize(int(shape.x) * int(shape.y) * int(shape.z), 0);

		if (initial_geometry.size() != _size_x * _size_y * _size_z)
		{
			//throw std::invalid_argument("initial geometry of wrong shape");
			std::clog << "Warning: initial geometry of wrong shape\n";
		}

	
	_mat_grid = initial_geometry;
}

template<bool gpu_flag>
CPU voxels<gpu_flag> voxels<gpu_flag>::create(std::vector<triangle> const & triangles)
{
	const real VOXEL_SIZE = 0.3; // voxel size in nanometers (0.27 nm is appr. one atom of Si)

	const int SIZE_X = 201; // horizontal size in the x direction in voxels
	const int SIZE_Y = 201; // horizontal size in the y direction in voxels
	const int SIZE_Z = 700; // vertical size in voxels

	const int SAMPLE_HEIGHT = 300; // height of the sample (length between the sample and the top of the simulation domain) in voxels
	
	vec3 shape = { SIZE_X, SIZE_Y, SIZE_Z };
	
	// Hier moet ergens een grid gemaakt worden, 
		// een voor het materiaal, 
		// een voor the timestamp van de electronen, 
		// een voor de energy van het deposition electron,
		// en een voor de dz van dat electron
	
	 // set the shape 

	// set the initial geometry
	std::vector<int> ini_geom;
	

	ini_geom.resize(SIZE_X * SIZE_Y * SIZE_Z, 0);
	
	for (int i = 0; i < SIZE_X; i++) {
		for (int j = 0; j < SIZE_Y; j++) {
			for (int k = 0; k < SAMPLE_HEIGHT; k++) {
				ini_geom.at( i + j * SIZE_X + k * SIZE_X * SIZE_Y) = -123;
			}
		}
	}

	for (int i = 0; i < SIZE_X; i++) {
		for (int j = 0; j < SIZE_Y; j++) {
			ini_geom.at(i + j * SIZE_X + (SAMPLE_HEIGHT + 4) * SIZE_X * SIZE_Y) = -126;
		}
	}
	
	

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
	
	voxels<false> geometry(VOXEL_SIZE, shape, ini_geom);
	return geometry;
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
	triangle const* ignore_triangle, int ignore_material) const
{
	intersect_event evt{ distance, nullptr };

	real x = start.x / _voxel_size; // create vars for the location elements
	real y = start.y / _voxel_size;
	real z = start.z / _voxel_size;
	
	real dx = direction.x; // create vars for the direction elements
	real dy = direction.y;
	real dz = direction.z;

	vec3 dr = { dx, dy, dz };
	vec3 delta_S = { 0, 0, 0 };

	real delta_s_min = distance / _voxel_size;

	int start_mat = ignore_material;

	if(start_mat == -122)
	{
		std::clog << "\n  Mirror!!!  " << start_mat << "\n";
	}

	real delta_x;
	if (dx > 0) {
		delta_x = std::ceil(x) - x;
		delta_S.x = delta_x / dx;
	}
	else if(dx < 0){
		delta_x = x - std::floor(x);
		delta_S.x = -delta_x / dx;
	}
	else // dx == 0
	{
		delta_S.x = distance / _voxel_size;
	}

	real delta_y;
	if (dy > 0) {
		delta_y = std::ceil(y) - y;
		delta_S.y = delta_y / dy;
	}
	else if(dy < 0){
		delta_y = y - std::floor(y);
		delta_S.y = -delta_y / dy;
	}
	else // dy == 0
	{
		delta_S.y = distance / _voxel_size;
	}

	real delta_z;
	if (dz > 0) {
		delta_z = std::ceil(z) - z;
		delta_S.z = delta_z / dz;
	}
	else if(dz < 0) {
		delta_z = z - std::floor(z);
		delta_S.z = -delta_z / dz;
	}
	else // dz == 0
	{
		delta_S.z = distance / _voxel_size;
	}

	while (distance / _voxel_size >= delta_s_min) {

		//std::clog << "\n" << delta_s_min ;

		// Determine minimum from delta_S
		delta_s_min = std::min(std::min(delta_S.x, delta_S.y), delta_S.z);

		if(delta_s_min < 0.000001)
		{
			std::clog << "\n0 in delta_s_min " << delta_S.x << "  " << delta_S.y << "  " << delta_S.z;
			if(delta_S.x < 0.000001)
			{
				delta_S.x += std::abs(1 / dx);
			}
			if (delta_S.y < 0.000001)
			{
				delta_S.y += std::abs(1 / dy);
			}
			if (delta_S.z < 0.000001)
			{
				delta_S.z += std::abs(1 / dz);
			}
			delta_s_min = std::min(std::min(delta_S.x, delta_S.y), delta_S.z);
		}

		int min_index = 0;
		if(delta_s_min == delta_S.y)
		{
			min_index = 1;
		}
		else if (delta_s_min == delta_S.z)
		{
			min_index = 2;
		}

		//std::clog << "   " << min_index << "   " << distance / _voxel_size;
		const int min_i = min_index;

		vec3 new_pos = start / _voxel_size + (delta_s_min + 0.001) * dr; // new position in voxels

		vec3 pos = new_pos * _voxel_size; // new position in nm, for check

		if(!((pos.x > _AABB_min.x) && (pos.x < _AABB_max.x)
			&& (pos.y > _AABB_min.y) && (pos.y < _AABB_max.y)
			&& (pos.z > _AABB_min.z) && (pos.z < _AABB_max.z))) // if out of range, return
		{
			return evt;
		}

		std::clog << "\nposition: " << new_pos.x << "  " << new_pos.y << "  " << new_pos.z;

		int k;
		int l;
		int m;

		real dx_sgn = 0;
		real dy_sgn = 0;
		real dz_sgn = 0;
		
		//std::clog << "\n" << dx_sgn << "   " << dy_sgn << "   " << dz_sgn;
		switch (min_i)
		{
		case 0: // intersection with x-plane
			dx_sgn = 0; // dx / std::abs(dx);
			k = (int)std::floor(new_pos.x + 0.1 * dx_sgn);
			l = (int)std::floor(new_pos.y);
			m = (int)std::floor(new_pos.z);
			
			if (k >= _size_x || k < 0)
			{
				return evt;
			}
			break;

		case 1: // intersection with y-plane
			dy_sgn = 0; // dy / std::abs(dy);
			k = (int)std::floor(new_pos.x);
			l = (int)std::floor(new_pos.y + 0.1 * dy_sgn);
			m = (int)std::floor(new_pos.z);
			
			if (l >= _size_y || l < 0)
			{
				return evt;
			}
			break;

		default: // intersection with z-plane
			dz_sgn = 0; //dz / std::abs(dz);
			k = (int)std::floor(new_pos.x);
			l = (int)std::floor(new_pos.y);
			m = (int)std::floor(new_pos.z + 0.1 *dz_sgn);

			if(m >= _size_z || m < 0)
			{
				return evt;
			}
			break;
		}		
		
		int new_mat = _mat_grid.at(k + l * _size_x + m * _size_x * _size_y);

		//std::clog << "   " << new_mat << "   " << start_mat;
		
		if (new_mat != start_mat) {

			std::clog << "\nintersection from " << start_mat << " to " << new_mat;
			
			evt.isect_distance = delta_s_min * _voxel_size;

			// Determine voxel side
			int	voxel_side;
			switch (min_i)
			{
			case 0:
				if (dx > 0) {
					voxel_side = 1;
				}
				else {
					voxel_side = 2;
				}
				break;

			case 1:
				if (dy > 0) {
					voxel_side = 3;
				}
				else {
					voxel_side = 4;
				}
				break;

			default:
				if (dz > 0) {
					voxel_side = 5;
				}
				else {
					voxel_side = 6;
				}
				break;
			}

			uint64_t isect_id;

			reinterpret_cast<int32_t*>(&isect_id)[0] = new_mat;
			reinterpret_cast<int32_t*>(&isect_id)[1] = voxel_side;
			evt.isect_triangle = reinterpret_cast<triangle*>(isect_id);

			return evt;
		}

		// Calculate new value of delta_S_i for next iteration
		switch (min_i)
		{
		case 0:
			delta_S.x += std::abs(1 / dx);
			break;

		case 1:
			delta_S.y += std::abs(1 / dy);
			break;

		default:
			delta_S.z += std::abs(1 / dz);
			break;
		}
	}
	return evt;
	
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

			std::vector<int> a;
			voxels_t geometry(3, {5, 5, 6}, a);

			
			/*if (triangles.size() > std::numeric_limits<triangle_index_t>::max())
				throw std::runtime_error("Too many triangles in geometry");
			geometry._N = static_cast<triangle_index_t>(triangles.size());

			geometry._triangles = reinterpret_cast<triangle*>(malloc(geometry._N * sizeof(triangle)));
			for (triangle_index_t i = 0; i < triangles.size(); ++i)
			{
				geometry._triangles[i] = triangles[i];
			}*/

			//geometry.set_AABB(AABB_min, AABB_max);

			return geometry;
		}

		inline static CPU void free(voxels<false> & geometry)
		{
			//::free(geometry._triangles);

			//geometry._triangles = nullptr;
			//geometry._N = 0;
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
