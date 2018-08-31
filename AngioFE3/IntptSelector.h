#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FEElement.h>
#include "FEAngio.h"

class IntptSelector : public FEMaterial
{
public:
	explicit IntptSelector(FEModel* pfem) : FEMaterial(pfem) {}
	virtual void SelectIntPtr(std::vector<FEMaterialPoint*> & int_pts, FESolidElement * se, vec3d natc, FEAngio* feangio, FEMesh* mesh) = 0;
protected:
	double safety_factor = 1.1;
};


class NarrowIntptSelector : public IntptSelector
{
public:
	explicit NarrowIntptSelector(FEModel* pfem) : IntptSelector(pfem) {}
	void SelectIntPtr(std::vector<FEMaterialPoint*> & int_pts, FESolidElement * se, vec3d natc, FEAngio* feangio, FEMesh* mesh)  override;

};

class AdjacentElementIntptSelector : public IntptSelector
{
public:
	explicit AdjacentElementIntptSelector(FEModel* pfem) : IntptSelector(pfem) {}
	void SelectIntPtr(std::vector<FEMaterialPoint*> & int_pts, FESolidElement * se, vec3d natc, FEAngio* feangio, FEMesh* mesh)  override;

};