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
// TODO: Currently assigns both global material orientation and angio specific orientation. Will need to pick one or have separate functions.
void FiberRandomizer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	//get the FE domain
	FEDomain* Dom = dynamic_cast<FEDomain*>(angio_element->_elem->GetMeshPartition());
	//
	FEElementSet* elset = mesh->FindElementSet(Dom->GetName());
	int local_index = elset->GetLocalIndex(*angio_element->_elem);
	
	FEMaterial* Mat_a = Dom->GetMaterial();
	// assumes that materials mat_axis is already mapped which we'll need to do somewhere else.
	FEParam* matax = Mat_a->FindParameter("mat_axis");
	FEParamMat3d& p = matax->value<FEParamMat3d>();
	FEMappedValueMat3d* val = dynamic_cast<FEMappedValueMat3d*>(p.valuator());
	FEDomainMap* map = dynamic_cast<FEDomainMap*>(val->dataMap());


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
		map->setValue(local_index, i, RandMat);
	}
}

// take a given angio element and give it a randomized discrete fiber direction based on user input spa
void DiscreteFiberEFDRandomizer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	efd_axes_a.unit();
	efd_axes_b.unit();
	efd_axes_c.unit();
	mat3d elem_dir;
	elem_dir.setCol(0, efd_axes_a*efd_spa.x);
	elem_dir.setCol(1, efd_axes_b*efd_spa.y);
	elem_dir.setCol(2, efd_axes_c*efd_spa.z);

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

	// for each integration point in the element
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		// get the material point of the integration point
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		// get the angio material point
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);

		double theta_12 = angio_element->GetEllipseAngle(r0, r1,0,2*PI,360);
		double theta_13 = angio_element->GetEllipseAngle(r0, r2,0,2*PI,360);
		// rotate the primary direction by theta_12 about the normal between them
		vec3d axis = mix3d_t(axis_0, axis_1, theta_12); axis.unit();
		mat3d R12 = mix3d_t_r(axis_0, axis_1, theta_12);
		//angio_pt->angio_fd = mix3d_t(axis, axis_2, theta_13); angio_pt->angio_fd.unit();
		mat3d R13 = mix3d_t_r(axis, axis_2, theta_13);

		//get the FE domain
		FEDomain* Dom = dynamic_cast<FEDomain*>(angio_element->_elem->GetMeshPartition());
		//
		FEElementSet* elset = mesh->FindElementSet(Dom->GetName());
		int local_index = elset->GetLocalIndex(*angio_element->_elem);
		FEMaterial* Mat_a = Dom->GetMaterial();
		// assumes that materials mat_axis is already mapped which we'll need to do somewhere else.
		FEParam* matax = Mat_a->FindParameter("mat_axis");
		FEParamMat3d& p = matax->value<FEParamMat3d>();
		FEMappedValueMat3d* val = dynamic_cast<FEMappedValueMat3d*>(p.valuator());
		FEDomainMap* map = dynamic_cast<FEDomainMap*>(val->dataMap());

		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		// get local domain index of element
		mat3d temp_mat;
		temp_mat.setCol(0, axis_0); temp_mat.setCol(1, axis_1); temp_mat.setCol(2, axis_2);
		temp_mat = R13*R12*temp_mat;
		map->setValue(local_index, i, temp_mat);

		//FEElasticMaterialPoint * mat_emp = angio_pt->matPt->ExtractData<FEElasticMaterialPoint>();
		//mat_emp->m_Q = emp->m_Q;
		//// create a vessel material point and copy the elastic material point to it
		//FEElasticMaterialPoint * ves_emp = angio_pt->vessPt->ExtractData<FEElasticMaterialPoint>();
		//ves_emp->m_Q = emp->m_Q;
	}
}

BEGIN_FECORE_CLASS(DiscreteFiberEFDRandomizer, InitialModifier)
ADD_PARAMETER(efd_spa, "efd_spa");
ADD_PARAMETER(efd_axes_a, "efd_axes_a");
ADD_PARAMETER(efd_axes_b, "efd_axes_b");
ADD_PARAMETER(efd_axes_c, "efd_axes_c");
END_FECORE_CLASS();

void EFDFiberInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	// norm the initial axes
	initial_axes_a.unit();
	initial_axes_b.unit();
	initial_axes_c.unit();
	
	// scale initial axes by values	
	angio_element->initial_angioSPA.setCol(0, initial_axes_a*initial_spa.x);
	angio_element->initial_angioSPA.setCol(1, initial_axes_b*initial_spa.y);
	angio_element->initial_angioSPA.setCol(2, initial_axes_c*initial_spa.z);
	
	// store the deformed SPA
	angio_element->angioSPA = angio_element->initial_angioSPA;
	angio_element->UpdateAngioFractionalAnisotropy();
}

BEGIN_FECORE_CLASS(EFDFiberInitializer, InitialModifier)
ADD_PARAMETER(initial_axes_a, "initial_axes_a");
ADD_PARAMETER(initial_axes_b, "initial_axes_b");
ADD_PARAMETER(initial_axes_c, "initial_axes_c");
ADD_PARAMETER(initial_spa, "initial_spa");
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