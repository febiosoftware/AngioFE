#include "InitialModifiers.h"
#include "FEAngioMaterialPoint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include "FEAngio.h"
#include <iostream>
#include "angio3d.h"
#include "FECore/FEDomainMap.h"
//#include <FEBioMech/FEEllipsoidalFiberDistribution.h>
//#include <FEBioMech/FEFiberDensityDistribution.h>

void InitialModifierManager::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for(int i=0; i< initial_modifiers.size();i++)
	{
		initial_modifiers[i]->ApplyModifier(angio_element, mesh, feangio);
	}
}

// Rename this to AlignedFiberRandomizer or DiscreteFiberRandomizer
// take a given angio element and give it a randomized discrete fiber direction
void FiberRandomizer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	// for each integration point in the element
	for(int i=0; i < angio_element->_elem->GaussPoints();i++)
	{
		// get the material point of the integration point
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		// get the angio material point
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		// get the elastic material point
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		// multiply the elastic material point's global orientation by a random matrix
		mat3d RandMat = feangio->unifromRandomRotationMatrix(angio_element->_rengine);
		// get local domain index of element
		angio_pt->angio_fiber_dir = RandMat*angio_pt->angio_fiber_dir;
	}
}

// take a given angio element and give it a randomized discrete fiber direction based on user input spd
void DiscreteFiberEFDRandomizer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{	
	mat3d elem_dir;
	elem_dir.setCol(0, vec3d(m_SPD.xx(), m_SPD.xy(), m_SPD.xz()));
	elem_dir.setCol(1, vec3d(m_SPD.xy(), m_SPD.yy(), m_SPD.yz()));
	elem_dir.setCol(1, vec3d(m_SPD.xz(), m_SPD.yz(), m_SPD.zz()));

	std::vector<pair<double, int>> v;
	v.push_back(pair<double, int>(elem_dir.col(0).norm(), 0));
	v.push_back(pair<double, int>(elem_dir.col(1).norm(), 1));
	v.push_back(pair<double, int>(elem_dir.col(2).norm(), 2));
	sort(v.begin(), v.end(), sortinrev);

	// store the indices
	int u1 = v[0].second;
	int u2 = v[1].second;
	int u3 = v[2].second;

	vec3d axis_0 = elem_dir.col(u1); axis_0.unit();
	vec3d axis_1 = elem_dir.col(u2); axis_1.unit();
	vec3d axis_2 = elem_dir.col(u3); axis_2.unit();
	double r0 = (elem_dir.col(u1).norm());
	double r1 = (elem_dir.col(u2).norm());
	double r2 = (elem_dir.col(u3).norm());	
	FEEllipticalDistribution E0(this->GetFEModel());
	FEEllipticalDistribution E1(this->GetFEModel());
	E0.a = r0; E0.b = r1; E0.Init();
	E1.a = r0; E1.b = r2; E1.Init();

	// for each integration point in the element
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		double theta_12 = E0.NextValue(angio_element->_rengine);
		double theta_13 = E1.NextValue(angio_element->_rengine);
		//double theta_12 = angio_element->GetEllipseAngle(r0, r1,0,2*PI,360);
		//double theta_13 = angio_element->GetEllipseAngle(r0, r2,0,2*PI,360);
		// rotate the primary direction by theta_12 about the normal between them
		vec3d axis = mix3d_t(axis_0, axis_1, theta_12); axis.unit();
		mat3d R12 = mix3d_t_r(axis_0, axis_1, theta_12);
		//angio_pt->angio_fd = mix3d_t(axis, axis_2, theta_13); angio_pt->angio_fd.unit();
		mat3d R13 = mix3d_t_r(axis, axis_2, theta_13);
		
		// get the material point of the integration point
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		// get the angio material point
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		//Rotate by the sampled direction
		angio_pt->angio_fiber_dir = R13*R12*vec3d(1, 0, 0);
	}
	// needs to be averaged between integration points for this element
}

BEGIN_FECORE_CLASS(DiscreteFiberEFDRandomizer, InitialModifier)
ADD_PARAMETER(m_SPD, "spd");
END_FECORE_CLASS();

// take a given angio element and give it a randomized discrete fiber direction based on user input spd
void DiscreteFiberEFDMatRandomizer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	std::vector<pair<double, int>> v;
	mat3d ax;
	ax.setCol(0, vec3d(m_SPD.xx(), m_SPD.xy(), m_SPD.xz()));
	ax.setCol(1, vec3d(m_SPD.xy(), m_SPD.yy(), m_SPD.yz()));
	ax.setCol(2, vec3d(m_SPD.xz(), m_SPD.yz(), m_SPD.zz()));
	v.push_back(pair<double, int>(ax.col(0).norm(), 0));
	v.push_back(pair<double, int>(ax.col(1).norm(), 1));
	v.push_back(pair<double, int>(ax.col(2).norm(), 2));
	sort(v.begin(), v.end(), sortinrev);

	// store the indices
	int i = v[0].second;
	int j = v[1].second;
	int k = v[2].second;

	vec3d axis_0 = ax.col(i); axis_0.unit();
	vec3d axis_1 = ax.col(j); axis_1.unit();
	vec3d axis_2 = ax.col(k); axis_2.unit();
	double r0 = (ax.col(i).norm());
	double r1 = (ax.col(j).norm());
	double r2 = (ax.col(k).norm());
	FEEllipticalDistribution E0(this->GetFEModel());
	FEEllipticalDistribution E1(this->GetFEModel());
	E0.a = r0; E0.b = r1; E0.Init();
	E1.a = r0; E1.b = r2; E1.Init();

	// for each integration point in the element
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		// get the material point of the integration point
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		// get the angio material point
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		double theta_12 = E0.NextValue(angio_element->_rengine);
		double theta_13 = E1.NextValue(angio_element->_rengine);
		//double theta_12 = angio_element->GetEllipseAngle(r0, r1, 0, 2*PI, 360);
		//double theta_13 = angio_element->GetEllipseAngle(r0, r2, 0, 2*PI, 360);
		// rotate the primary direction by theta_12 about the normal between them
		vec3d axis = mix3d_t(axis_0, axis_1, theta_12); axis.unit();
		mat3d R12 = mix3d_t_r(axis_0, axis_1, theta_12);
		mat3d R13 = mix3d_t_r(axis, axis_2, theta_13);
		angio_pt->angio_fiber_dir = R13*R12*vec3d(1, 0, 0);
	}
}

BEGIN_FECORE_CLASS(DiscreteFiberEFDMatRandomizer, InitialModifier)
ADD_PARAMETER(m_SPD, "spd");
END_FECORE_CLASS();

void EFDFiberInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	//// norm the initial axes
	//initial_axes_a.unit();
	//initial_axes_b.unit();
	//initial_axes_c.unit();
	//
	//// scale initial axes by values	
	//angio_element->initial_angioSPD.setCol(0, initial_axes_a*initial_spd.x);
	//angio_element->initial_angioSPD.setCol(1, initial_axes_b*initial_spd.y);
	//angio_element->initial_angioSPD.setCol(2, initial_axes_c*initial_spd.z);
	
	//double d[3]; vec3d r[3];
	//p.eig();
	/*mat3ds ax = mat3ds(p.x)
	angio_element->initial_angioSPD.setCol(0, initial_axes_a*initial_spd.x);
	angio_element->initial_angioSPD.setCol(1, initial_axes_b*initial_spd.y);
	angio_element->initial_angioSPD.setCol(2, initial_axes_c*initial_spd.z);*/


	// store the SPD
	angio_element->initial_angioSPD = m_SPDa*(3/m_SPDa.tr());
	angio_element->angioSPD = m_SPDa*(3 / m_SPDa.tr());
	angio_element->UpdateAngioFractionalAnisotropy();
}

BEGIN_FECORE_CLASS(EFDFiberInitializer, InitialModifier)
ADD_PARAMETER(m_SPDa, "spd");
//ADD_PARAMETER(initial_axes_a, "initial_axes_a");
//ADD_PARAMETER(initial_axes_b, "initial_axes_b");
//ADD_PARAMETER(initial_axes_c, "initial_axes_c");
//ADD_PARAMETER(initial_spd, "initial_spd");
END_FECORE_CLASS();

void EFDMatPointFiberInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	FESolidElement* se = angio_element->_elem;
	for (auto i = 0; i < se->GaussPoints(); i++)
	{
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		angio_pt->initial_angioSPD = m_SPD;
		angio_pt->angioSPD = m_SPD;
		//std::cout << angio_pt->angioSPD.xx() << ", " << angio_pt->angioSPD.yy() << ", " << angio_pt->angioSPD.zz() << ", " << angio_pt->angioSPD.xy() << ", " << angio_pt->angioSPD.xz() << ", " << angio_pt->angioSPD.yz() << endl;
		angio_pt->UpdateAngioFractionalAnisotropy();
		//std::cout << angio_pt->angioSPD.xx() << ", " << angio_pt->angioSPD.yy() << ", " << angio_pt->angioSPD.zz() << ", " << angio_pt->angioSPD.xy() << ", " << angio_pt->angioSPD.xz() << ", " << angio_pt->angioSPD.yz() << endl;
	}
}

BEGIN_FECORE_CLASS(EFDMatPointFiberInitializer, InitialModifier)
ADD_PARAMETER(m_SPD, "spd");
END_FECORE_CLASS();

void DensityInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint *mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint *angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);

		angio_pt->ref_ecm_density = initial_density;

	}
}

BEGIN_FECORE_CLASS(DensityInitializer, InitialModifier)
ADD_PARAMETER(initial_density, "initial_density");
END_FECORE_CLASS();

void RepulseInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint *mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint *angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);

		angio_pt->repulse_value = initial_repulse_value;
	}
}

BEGIN_FECORE_CLASS(RepulseInitializer, InitialModifier)
ADD_PARAMETER(initial_repulse_value, "initial_repulse_value");
END_FECORE_CLASS();