#include "FiberManager.h"
#include <FEBioMech/FESPRProjection.h>
#include <FEBioMech/FEFiberMaterialPoint.h>
#include "FEAngioMaterial.h"
#include "angio3d.h"


void FiberManager::Update()
{

}

vec3d FiberManager::GetFiberAtNode(int node, double & lambda)
{
	vec3d rv;

	rv.x = fibers_at_nodes[0][node];
	rv.y = fibers_at_nodes[1][node];
	rv.z = fibers_at_nodes[2][node];
	lambda = fibers_at_nodes[3][node];
	rv.unit();
	return rv;
}
vec3d FiberManager::GetMinor1AtNode(int node, double  &lambda)
{
	vec3d rv;

	rv.x = minoraxis1_at_nodes[0][node];
	rv.y = minoraxis1_at_nodes[1][node];
	rv.z = minoraxis1_at_nodes[2][node];
	lambda = minoraxis1_at_nodes[3][node];

	rv.unit();
	return rv;
}
vec3d FiberManager::GetMinor2AtNode(int node, double  &lambda)
{
	vec3d rv;

	rv.x = minoraxis2_at_nodes[0][node];
	rv.y = minoraxis2_at_nodes[1][node];
	rv.z = minoraxis2_at_nodes[2][node];
	lambda = minoraxis2_at_nodes[3][node];

	rv.unit();
	return rv;
}

void RandomFiberInitializer::InitializeFibers(FiberManager * fman)
{
	//set the fiber vectors at the nodes of the fiber manager
	//make sure the map has been generated
	
}


void RandomFiberInitializerNonMangling::InitializeFibers(FiberManager * fman)
{
	
}

void RandomFiberInitializerPE::InitializeFibers(FiberManager * fman)
{
	
}

ExplicitDistributionsFiberInitializer::ExplicitDistributionsFiberInitializer(FEModel * model) : FiberInitializer(model)
{
	AddProperty(&alpha, "alpha");
	AddProperty(&alpha, "beta");
	AddProperty(&alpha, "gamma");
}
void ExplicitDistributionsFiberInitializer::Setup()
{
	alpha->StepToTime(0);
	beta->StepToTime(0);
	gamma->StepToTime(0);
}

void ExplicitDistributionsFiberInitializer::InitializeFibers(FiberManager * fman)
{
	
}

vec3d EllipsoidPos(double a, double b, double c, double theta, double phi)
{
	return vec3d(a * cos(theta)*cos(phi),
		b*cos(theta)*sin(phi),
		c* sin(theta));
}

void EllipsoidalFiberInitializer::InitializeFibers(FiberManager * fman)
{

}

EllipsoidalFiberInitializer::EllipsoidalFiberInitializer(FEModel * model) : FiberInitializer(model)
{

}

BEGIN_PARAMETER_LIST(EllipsoidalFiberInitializer, FiberInitializer)
ADD_PARAMETER(a, FE_PARAM_DOUBLE, "a");
ADD_PARAMETER(b, FE_PARAM_DOUBLE, "b");
ADD_PARAMETER(c, FE_PARAM_DOUBLE, "c");

ADD_PARAMETER(theta_slice , FE_PARAM_DOUBLE, "theta_slice");
ADD_PARAMETER(phi_slice, FE_PARAM_DOUBLE, "a");
END_PARAMETER_LIST();

