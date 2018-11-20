#include "rbf_norms.h"

/*
double gaussian_rbf_norm::rbfnorm(vec3d pos, vec3d center)
{
	double dist = (pos - center).norm();
	if(dist < cut_off_point)
	{
		return max_intpt_contribution;
	}
	double inverse_critical_radius = 1 / dist;
	inverse_critical_radius *= inverse_critical_radius;
	return exp(-inverse_critical_radius);
}
*/