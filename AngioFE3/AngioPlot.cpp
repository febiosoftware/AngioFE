#include "StdAfx.h"
#include "AngioPlot.h"
#include "FEAngio.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEModel.h>
#include "FEAngioMaterial.h"
#include "FECore/FEMesh.h"
#include <unordered_map>
#include <algorithm>
#include <FEBioMech/FEViscoElasticMaterial.h>
#include <FEBioMix/FEMultiphasic.h>
#include "FEBioMix/FEBiphasicSolute.h"
#include "FEBioMix/FETriphasic.h"
#include "FEBioMix/FEMultiphasicSolidDomain.h"
#include "FEBioMix/FEMultiphasicShellDomain.h"
#include "FECore/FEDomainMap.h"
#include <iostream>
#include <FECore/mathalg.h>
#include "angio3d.h"

extern FEAngio* pfeangio;

//-----------------------------------------------------------------------------
bool FEPlotAngioStress::Save(FEDomain& d, FEDataStream& str)
{
	// get the angio component of the material and check if it is an angio material
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	// get the solid domain from the FE domain
	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	// for each element get the stress
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		// for each gauss point get the stress
		for (int j=0; j<nint; ++j)
		{
			// get the material point of the gauss point
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			// get the angio material point of the material point to get the stress then add it to the stress sum
			FEAngioMaterialPoint * angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			mat3ds & sj = angio_mp->m_as;
			s += sj;
		}
		// average sum of gauss point stress
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();
			mat3ds & sj = matrix_elastic.m_s;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotVesselStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			// get elastic stress from vessel
			FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
			mat3ds & sj = vessel_elastic.m_s;
			
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotVesselWeight::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		double s = 0.0;
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			s += angioPt->vessel_weight;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixWeight::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		double s = 0.0;
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			s += angioPt->matrix_weight;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixTangent::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		tens4ds s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			// evaluate the tangent at the material point
			tens4ds ten = pmat->GetMatrixMaterial()->ExtractProperty<FEElasticMaterial>()->Tangent(mp);
			tens4ds sj = ten;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixViscoStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			FEViscoElasticMaterialPoint& matrix_visco_elastic = *angioPt->matPt->ExtractData<FEViscoElasticMaterialPoint>();
			
			mat3ds sj =  matrix_visco_elastic.m_se;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixElasticStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();	
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			FEElasticMaterialPoint& emp = *angioPt->matPt->Next()->ExtractData<FEElasticMaterialPoint>();

			mat3ds sj = emp.m_s;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixElastic_m_Q::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3d s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint * mp = (el.GetMaterialPoint(j));
			FEElasticMaterialPoint*  emp = mp->ExtractData<FEElasticMaterialPoint>();
			FEMesh* mesh = d.GetMesh();
			//get the FE domain
			FEElementSet* elset = mesh->FindElementSet(d.GetName());
			int local_index = elset->GetLocalIndex(el);
			//Ask Steve about this
			FEMaterial* Mat_a = d.GetMaterial()->ExtractProperty<FEElasticMaterial>();
			// assumes that materials mat_axis is already mapped which we'll need to do somewhere else.
			FEParam* matax = Mat_a->FindParameter("mat_axis");
			FEParamMat3d& p = matax->value<FEParamMat3d>();
			FEMappedValueMat3d* val = dynamic_cast<FEMappedValueMat3d*>(p.valuator());
			FEDomainMap* map = dynamic_cast<FEDomainMap*>(val->dataMap());
			mat3d sj = map->valueMat3d(emp);

			//mat3d sj = emp->m_Q;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

bool FEPlotBranches::Save(FEDomain& d, FEDataStream& str)
{
	// for each element
	for(int i=0; i < d.Elements();i++)
	{
		FEElement & elem = d.ElementRef(i);
		// get teh solid element
		FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
		if(se)
		{
			// convert solid element to angio element
			AngioElement * angio_element = pfeangio->se_to_angio_element.at(se);
			// get the branch count for the element and output to the stream
			str << static_cast<double>(angio_element->branch_count);
		}
	}
	
	return true;
};

bool FEPlotAnastamoses::Save(FEDomain& d, FEDataStream& str)
{
	for (int i = 0; i < d.Elements(); i++)
	{
		FEElement & elem = d.ElementRef(i);
		FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
		if (se)
		{
			AngioElement * angio_element = pfeangio->se_to_angio_element.at(se);
			str << static_cast<double>(angio_element->anastamoses);
		}
	}

	return true;
};

bool FEPlotSegmentLength::Save(FEDomain& d, FEDataStream& str)
{
	for (int i = 0; i < d.Elements(); i++)
	{
		FEElement & elem = d.ElementRef(i);
		FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
		if (se)
		{
			AngioElement * angio_element = pfeangio->se_to_angio_element.at(se);
			str << angio_element->global_segment_length;
		}
	}

	return true;
};

bool FEPlotRefSegmentLength::Save(FEDomain& d, FEDataStream& str)
{
	for (int i = 0; i < d.Elements(); i++)
	{
		FEElement & elem = d.ElementRef(i);
		FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
		if (se)
		{
			AngioElement * angio_element = pfeangio->se_to_angio_element.at(se);
			str << angio_element->refernce_frame_segment_length;
		}
	}

	return true;
};

//-----------------------------------------------------------------------------
bool FEPlotAngioECMDensity::Save(FEDomain& d, FEDataStream& str)
{
	//Check if the domain has an angio component i.e. make sure this is not a rigid body
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat != nullptr)
	{
		for (int i = 0; i < d.Elements(); i++)
		{
			FEElement & elem = d.ElementRef(i);
			FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
			if (se)
			{
				AngioElement* angio_element = pfeangio->se_to_angio_element.at(se);
				double den = 0.0;
				for (int j = 0; j < angio_element->_elem->GaussPoints(); j++)
					{
						FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(j);
						FEAngioMaterialPoint *angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
						FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
						den += angio_pt->ref_ecm_density* (1.0 / emp->m_J);
						}
						den /= angio_element->_elem->GaussPoints();
						str << den;
					}

			}
	}
		return true;
}

//-----------------------------------------------------------------------------
bool FEPlotAngioRepulseVal::Save(FEDomain& d, FEDataStream& str)
{
	//Check if the domain has an angio component i.e. make sure this is not a rigid body
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat != nullptr)
	{
		for (int i = 0; i < d.Elements(); i++)
		{
			FEElement & elem = d.ElementRef(i);
			FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
			if (se)
			{
				AngioElement* angio_element = pfeangio->se_to_angio_element.at(se);
				double repulse_val = 0.0;
				for (int j = 0; j < angio_element->_elem->GaussPoints(); j++)
				{
					FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(j);
					FEAngioMaterialPoint *angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
					FEElasticMaterialPoint* elastic_mp = mp->ExtractData<FEElasticMaterialPoint>();
					repulse_val += angio_mp->repulse_value * (1.0 / elastic_mp->m_J);
				}
				repulse_val /= angio_element->_elem->GaussPoints();
				str << repulse_val;
			}

		}
	}
	return true;
}

bool FEPlotangioSPD::Save(FEDomain& d, FEDataStream& str)
{
	//Check if the domain has an angio component i.e. make sure this is not a rigid body
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat != nullptr)
	{
		// for each element
		for (int i = 0; i < d.Elements(); i++)
		{
			FEElement & elem = d.ElementRef(i);
			FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
			if (se) {
				// get the pointer to the angio element
				AngioElement* angio_element = pfeangio->se_to_angio_element.at(se);

				// store the transformed spd
				// should go somewhere else...
				angio_element->UpdateSPD();
				str << angio_element->angioSPD;
			}
		}
	}
	return true;
}

bool FEPlotAngioFractionalAnisotropy::Save(FEDomain& d, FEDataStream& str)
{
	for (int i=0; i < d.Elements(); i++)
	{
		FEElement & elem = d.ElementRef(i);
		FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
		if (se) {
			// get pointer to the angio element
			AngioElement* angio_element = pfeangio->se_to_angio_element.at(se);
			// store the fractional anisotropy.
			// should go somewhere else...
			angio_element->UpdateAngioFractionalAnisotropy();
			str << angio_element->angioFA;
		}
	}
	return true;
}

bool FEPlotMatangioSPD::Save(FEDomain& d, FEDataStream& str)
{
	for (int i = 0; i < d.Elements(); i++)
	{
		FEElement & elem = d.ElementRef(i);
		FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
		if (se) {
			// get pointer to the angio element
			AngioElement* angio_element = pfeangio->se_to_angio_element.at(se);
			// vector containing the SPD for each gauss point in the element
			mat3ds SPDs_gausspts[8];
			double H[8];
			// get each gauss point's SPD
			for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
			{
				// get the angio point
				FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
				FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
				// Get the SPD   
				angio_mp->UpdateSPD();
				mat3ds temp_SPD = angio_mp->angioSPD;
				SPDs_gausspts[i] = temp_SPD*(3.0/temp_SPD.tr());
				// TODO: calculate all distances from mp to nodes then normalize. Currently assumes equidistance
				H[i] = 1.0 / angio_element->_elem->GaussPoints();
			}
			// array containing the SPD for each node in the element
			mat3ds SPDs_nodes[8];
			// Get the interpolated SPD from the shape function-weighted Average Structure Tensor
			mat3ds SPD_int = weightedAverageStructureTensor(SPDs_gausspts, H, angio_element->_elem->GaussPoints());
			SPD_int = (3.0 / SPD_int.tr())*SPD_int;
			str << SPD_int;
		}
	}
	return true;
}

bool FEPlotMatAngioFractionalAnisotropy::Save(FEDomain& d, FEDataStream& str)
{
	for (int i = 0; i < d.Elements(); i++)
	{
		FEElement & elem = d.ElementRef(i);
		FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
		if (se) {
			// get pointer to the angio element
			AngioElement* angio_element = pfeangio->se_to_angio_element.at(se);
			// vector containing the SPD for each gauss point in the element
			mat3ds SPDs_gausspts[FEElement::MAX_INTPOINTS];
			double H[FEElement::MAX_INTPOINTS ];
			// get each gauss point's SPD
			for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
			{
				// get the angio point
				FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
				FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);

				// Get the SPD
				angio_mp->UpdateSPD();
				SPDs_gausspts[i] = angio_mp->angioSPD;
				H[i] = 1.0/angio_element->_elem->GaussPoints();
			}
			// array containing the SPD for each node in the element
			mat3ds SPDs_nodes[FEElement::MAX_NODES];
			// array for the shape function values
			
			// project the spds from integration points to the nodes
			angio_element->_elem->project_to_nodes(&SPDs_gausspts[0], SPDs_nodes);
			// Get the interpolated SPD from the shape function-weighted Average Structure Tensor
			//mat3ds SPD_int = weightedAverageStructureTensor(SPDs_gausspts, H, angio_element->_elem->GaussPoints());
			mat3ds SPD_int = weightedAverageStructureTensor(SPDs_nodes, H, angio_element->_elem->GaussPoints());
			
			// get the vectors of the principal directions and sort in descending order
			std::vector<pair<double, int>> v;
			mat3d ax;
			ax.setCol(0, vec3d(SPD_int.xx(), SPD_int.xy(), SPD_int.xz()));
			ax.setCol(1, vec3d(SPD_int.xy(), SPD_int.yy(), SPD_int.yz()));
			ax.setCol(2, vec3d(SPD_int.xz(), SPD_int.yz(), SPD_int.zz()));
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

			// calculate the fractional anisotropy
			double angioFA_int = sqrt(0.5)*(sqrt(pow(r0 - r1, 2) + pow(r1 - r2, 2) + pow(r2 - r0, 2)) / (sqrt(pow(r0, 2) + pow(r1, 2) + pow(r2, 2))));
			str << angioFA_int;
		}
	}
	return true;
}

bool FEPlotEFDFiberDirection::Save(FEDomain& d, FEDataStream& str)
{
	FEMesh * mesh = pfeangio->GetMesh();
	for (int i = 0; i < d.Elements(); i++)
	{
		FEElement & elem = d.ElementRef(i);
		FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
		if (se)
		{
			AngioElement * angio_element = pfeangio->se_to_angio_element.at(se);
			vec3d primary_dir;
			std::vector<quatd> gauss_data;
			std::vector<quatd> nodal_data;

			for (int i = 0; i< angio_element->_elem->GaussPoints(); i++)
			{
				FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
				FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
				vec3d axis(1, 0, 0);

				//Ask Steve about this
				//get the FE domain
				FEDomain* Dom = dynamic_cast<FEDomain*>(angio_element->_elem->GetMeshPartition());
				//
				FEElementSet* elset = mesh->FindElementSet(Dom->GetName());
				int local_index = elset->GetLocalIndex(*angio_element->_elem);

				FEAngioMaterial* Mat_ang = Dom->GetMaterial()->ExtractProperty<FEAngioMaterial>();
				FEMaterial * Mat_a = Mat_ang->GetMatrixMaterial();
				mat3d m_Q = Mat_a->GetLocalCS(*mp);

				axis = emp->m_F * m_Q * axis;
				primary_dir += axis;
				//gauss_data.push_back({ axis });
			}
			primary_dir = primary_dir / 8; primary_dir.unit();
			str << primary_dir;
		}
	}
	return true;
}

bool FEPlotPrimaryVesselDirection::Save(FEDomain& d, FEDataStream& str)
{
	FEMesh * mesh = pfeangio->GetMesh();
	for (int i = 0; i < d.Elements(); i++)
	{
		FEElement & elem = d.ElementRef(i);
		FESolidElement *se = dynamic_cast<FESolidElement*>(&elem);
		if (se)
		{
			AngioElement * angio_element = pfeangio->se_to_angio_element.at(se);
			vec3d primary_dir;
			for (int j = 0; j < angio_element->grown_segments.size(); j++)
			{
				primary_dir += angio_element->grown_segments[j]->Direction(mesh)* angio_element->grown_segments[j]->Length(mesh);
			}
			str << primary_dir.norm();
		}
	}
	return true;
}

//stopgap plot implementation

//-----------------------------------------------------------------------------
// find the local solute ID, given a global ID. If the material is not a 
// biphasic-solute, triphasic, or multiphasic material, this returns -1.
int GetLocalSoluteID(FEMaterial* pm, int nsol)
{
	// figure out the solute ID to export. This depends on the material type.
	int nsid = -1;
	FEBiphasicSolute* psm = dynamic_cast<FEBiphasicSolute*> (pm);
	if (psm)
	{
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (psm->GetSolute()->GetSoluteID() == nsol);
		if (!present) return false;
		nsid = 0;
	}

	FETriphasic* ptm = dynamic_cast<FETriphasic*> (pm);
	if (ptm)
	{
		// Check if this solute is present in this specific triphasic mixture
		if (ptm->m_pSolute[0]->GetSoluteID() == nsol) nsid = 0;
		else if (ptm->m_pSolute[1]->GetSoluteID() == nsol) nsid = 1;
	}

	FEMultiphasic* pmm = dynamic_cast<FEMultiphasic*> (pm);
	if (pmm)
	{
		// Check if this solute is present in this specific multiphasic mixture
		for (int i = 0; i<pmm->Solutes(); ++i)
			if (pmm->GetSolute(i)->GetSoluteID() == nsol) { nsid = i; break; }
	}
	return nsid;
}