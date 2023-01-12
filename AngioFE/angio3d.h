///////////////////////////////////////////////////////////////////////
// angio3d.h
// Contains commonly used operations and utilities
///////////////////////////////////////////////////////////////////////


#pragma once

//! reflect incident vector where I is the incident vector and N is the normal 
//! vector to the surface
//see https://www.opengl.org/sdk/docs/man/html/reflect.xhtml
inline vec3d reflect(vec3d & I, vec3d & N)
{
	double temp = (2.0 * (N * I));
	return I - (N * temp);
}

// mix operation
//see https://www.opengl.org/sdk/docs/man/html/mix.xhtml
inline vec3d mix(vec3d & x, vec3d & y, double a)
{
	return x * (1.0 - a) + y * a;
}
//consider making this a template or using GLM
inline double mix(double x, double y, double a)
{
	return x * (1.0 - a) + y * a;
}

//new mix method 
//x is per, y is col_dir
// method of rotation about an axis: https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
inline vec3d mix3d(vec3d & x, vec3d & y, double a)
{
	// determine the angle between the two vectors then scale it
	double phi = a * acos(x * y);
	// determine the normal vector to the plane spanned by x and y
	vec3d normal_vec = x ^ y;
	normal_vec.unit();
	mat3d rot_mat;
	auto nx = normal_vec.x, ny = normal_vec.y, nz = normal_vec.z;
	double cp = cos(phi), sp = sin(phi);
	//! assemble the rotation matrix about the normal vector within the plane 
	//! spanned by x and y
	rot_mat[0][0] = cp + pow(nx, 2.0) * (1.0 - cp);
	rot_mat[0][1] = nx * ny * (1.0 - cp) - nz * sp;
	rot_mat[0][2] = nx * nz * (1.0 - cp) + ny * sp;
	rot_mat[1][0] = ny * nx * (1.0 - cp) + nz * sp;
	rot_mat[1][1] = cp + pow(ny, 2.0) * (1.0 - cp);
	rot_mat[1][2] = ny * nz * (1.0 - cp) - nx * sp;
	rot_mat[2][0] = nz * nx * (1.0 - cp) - ny * sp;
	rot_mat[2][1] = nz * ny * (1.0 - cp) + nx * sp;
	rot_mat[2][2] = cp + pow(nz, 2.0) * (1 - cp);

	return rot_mat * x;
}

//new mix method but the input t is an angle theta in radians.
//x is per, y is col_dir
// method of rotation about an axis: https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
inline vec3d mix3d_t(vec3d & x, vec3d & y, double t)
{
	// determine the normal vector to the plane spanned by x and y
	vec3d normal_vec = x ^ y;
	normal_vec.unit();
	mat3d rot_mat;
	auto nx = normal_vec.x, ny = normal_vec.y, nz = normal_vec.z;
	double ct = cos(t), st = sin(t);
	// assemble the rotation matrix about the normal vector within the plane spanned by x and y
	rot_mat[0][0] = ct + pow(nx, 2) * (1 - ct);	
	rot_mat[0][1] = nx * ny*(1 - ct) - nz * st;	
	rot_mat[0][2] = nx * nz*(1 - ct) + ny * st;
	rot_mat[1][0] = ny * nx*(1 - ct) + nz * st;	
	rot_mat[1][1] = ct + pow(ny, 2)*(1 - ct);	
	rot_mat[1][2] = ny * nz*(1 - ct) - nx * st;
	rot_mat[2][0] = nz * nx*(1 - ct) - ny * st;
	rot_mat[2][1] = nz * ny*(1 - ct) + nx * st;
	rot_mat[2][2] = ct + pow(nz, 2)*(1 - ct);

	return rot_mat * x;
}

// method of rotation about an axis: https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
// returns the rotation matrix instead. Can have the og function return a pair in future.
inline mat3d mix3d_t_r(vec3d & x, vec3d & y, double t)
{
	// determine the normal vector to the plane spanned by x and y
	vec3d normal_vec = x ^ y;
	normal_vec.unit();
	mat3d rot_mat;
	auto nx = normal_vec.x, ny = normal_vec.y, nz = normal_vec.z;
	double ct = cos(t), st = sin(t);
	// assemble the rotation matrix about the normal vector within the plane spanned by x and y
	rot_mat[0][0] = ct + pow(nx, 2) * (1 - ct);	
	rot_mat[0][1] = nx * ny*(1 - ct) - nz * st;	
	rot_mat[0][2] = nx * nz*(1 - ct) + ny * st;
	rot_mat[1][0] = ny * nx*(1 - ct) + nz * st;	
	rot_mat[1][1] = ct + pow(ny, 2)*(1 - ct);	
	rot_mat[1][2] = ny * nz*(1 - ct) - nx * st;
	rot_mat[2][0] = nz * nx*(1 - ct) - ny * st;	
	rot_mat[2][1] = nz * ny*(1 - ct) + nx * st;	
	rot_mat[2][2] = ct + pow(nz, 2)*(1 - ct);

	return rot_mat;
}

// Binary search used in volume/area seeders.
static size_t 
	findElement(double val, int lo, int high, double * begin, double * end)
{
	int mid = lo + (high - lo) / 2;
	if (val < begin[mid]) 
	{
		return findElement(val, lo, mid - 1, begin, end);
	} 
	else if (val > end[mid]) 
	{
		return findElement(val, mid + 1, high, begin, end);
	} 
	else 
	{
		return mid;
	}
}

inline bool sortinrev(const pair<double, int> &a, const pair<double, int> &b)
{
	return (a.first > b.first);
}