#include "StdAfx.h"
#include "ECMInitializer.h"
#include <FECore/FEMesh.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngioMaterialBase.h"
#include "FEAngio.h"

void ECMInitializer::updateECMdensity(FEAngioMaterialBase* mat)
{
	std::vector<int> matls;
	matls.emplace_back(mat->GetID_ang());
	FEMesh * mesh = mat->m_pangio->GetMesh();
	//break even on core in field model
	mat->m_pangio->ForEachElement([&](FESolidElement & elem, FESolidDomain & d)
	{
		//these will hold the natural coordinates once the project to nodes is complete 
		double nr[FEElement::MAX_NODES];
		double ns[FEElement::MAX_NODES];
		double nt[FEElement::MAX_NODES];
		//these hold the natural coordinates of the integration points (r,s,t)
		double gr[FEElement::MAX_NODES];
		double gs[FEElement::MAX_NODES];
		double gt[FEElement::MAX_NODES];
		//TODO: if needed get FEBIO to expose the vectors that contain these to avoid this copy
		for (int i = 0; i < elem.Nodes(); i++)
		{
			gr[i] = elem.gr(i);
			gs[i] = elem.gs(i);
			gt[i] = elem.gt(i);
		}

		elem.project_to_nodes(gr, nr);
		elem.project_to_nodes(gs, ns);
		elem.project_to_nodes(gt, nt);

		// For each node in the element...
		for (int j = 0; j<elem.Nodes(); ++j)
		{
			// get the node
			int nnum = elem.m_node[j];
			nnum = mesh->Node(nnum).GetID();
			// get the ecm density and collagen fiber
			//TODO

			/*
			//clamp n* to [1,-1]
			nr[j] = min(max(nr[j], -1), 1);
			ns[j] = min(max(ns[j], -1), 1);
			nt[j] = min(max(nt[j], -1), 1);
			*/

			//round to nearest integer
			nr[j] = round(nr[j]);
			ns[j] = round(ns[j]);
			nt[j] = round(nt[j]);

			// Calculate the deformation gradient tensor and jacobian at the node
			mat3d F;
			double Jacob;
			try {
				Jacob = d.defgrad(elem, F, nr[j], ns[j], nt[j]);
			}
			catch(NegativeJacobian nj)
			{
				Jacob = 1;
				F = mat3d::identity();
			}
			

			//make sure the function is differentiable and preserves orientation
			assert(Jacob > 0.0);

		}
	}, matls);
}

void ECMInitializerConstant::seedECMDensity(FEAngioMaterialBase * mat)
{
	std::vector<int> matls;
	matls.emplace_back(mat->GetID_ang());
	mat->m_pangio->ForEachNode([&](FENode & node)
	{

	}, matls);
}

void ECMInitializerSpecified::seedECMDensity(FEAngioMaterialBase * mat)
{
	std::vector<int> matls;
	matls.emplace_back(mat->GetID_ang());
	FEMesh * mesh = mat->m_pangio->GetMesh();
	mat->m_pangio->ForEachElement([this, mat, mesh](FESolidElement & se, FESolidDomain & sd)
	{
		double den[FEElement::MAX_INTPOINTS];
		double anis[FEElement::MAX_INTPOINTS];
		double pden[FEElement::MAX_INTPOINTS];
		double panis[FEElement::MAX_INTPOINTS];
		for (int n = 0; n<se.GaussPoints(); ++n)
		{
			// generate a coordinate transformation at this integration point
			FEMaterialPoint* mpoint = se.GetMaterialPoint(n);
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mpoint);
			den[n] = angioPt->m_D;
			anis[n] = angioPt->m_DA;
		}
		se.project_to_nodes(den, pden);
		se.project_to_nodes(anis, panis);
		for (int k = 0; k < se.Nodes(); k++)
		{
		}
	}, matls);
}

void ECMInitializerNoOverwrite::seedECMDensity(FEAngioMaterialBase * mat)
{
	std::vector<int> matls;
	matls.emplace_back(mat->GetID_ang());
	mat->m_pangio->ForEachNode([&](FENode & node)
	{

	}, matls);
}