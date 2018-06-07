///////////////////////////////////////////////////////////////////////
// angio3d.h
///////////////////////////////////////////////////////////////////////


#pragma once


//see https://www.opengl.org/sdk/docs/man/html/reflect.xhtml
inline vec3d reflect(vec3d & I, vec3d & N)
{
	double temp = (2.0 * (N * I));
	return I - (N * temp);
}
//see https://www.opengl.org/sdk/docs/man/html/mix.xhtml
inline vec3d mix(vec3d & x, vec3d & y, double a)
{
	return x*(1 - a) + y*a;
}
//consider making this a template or using GLM
inline double mix(double x, double y, double a)
{
	return x*(1 - a) + y*a;
}

// Binary search used in volume/area seeders.
static size_t findElement(double val, int lo, int high, double * begin, double * end)
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