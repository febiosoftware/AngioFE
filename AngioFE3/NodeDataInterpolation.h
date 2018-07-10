#pragma once
#include <FECore/FEMaterial.h>
#include "FEAngioMaterialPoint.h"
#include "VariableInterpolation.h"

class FEAngio;
class FEAngioMaterial;

class NodeDataInterpolation :public FEMaterial
{
public:
	NodeDataInterpolation(FEModel * pfem) :FEMaterial(pfem) { }
	virtual double & ValueReference(FEMaterialPoint * mp) = 0;
	virtual const char * GetDataName() const = 0;
	void PrepValues(FEAngio* angio, FEMesh* mesh, FEAngioMaterial* angio_mat);
protected:
	int node_set_id = -1;
	std::unordered_map<int, double> node_id_to_values;
	int interpolation_mode = 0;//0 use shape functions to interpolate to gauss points, 1 average all nodes in element and set that value at all gauss points
	DECLARE_PARAMETER_LIST();
};

class NodeDataInterpolationManager :public FEMaterial
{
public:
	NodeDataInterpolationManager(FEModel * pfem) :FEMaterial(pfem) { AddProperty(&node_data_interpolation_vals, "node_interpolation_value", false); }
	bool Init() override { return FEMaterial::Init(); }
	void DoInterpolations(FEAngio * angio, FEMesh* mesh, FEAngioMaterial* angio_mat);
private:
	FEVecPropertyT<NodeDataInterpolation> node_data_interpolation_vals;
};

class DensityValuesNodeDataInterpolation : public NodeDataInterpolation
{
public:
	DensityValuesNodeDataInterpolation(FEModel * pfem) :NodeDataInterpolation(pfem) {}
	bool Init() override { return FEMaterial::Init(); }
	double & ValueReference(FEMaterialPoint * mp) override;
	const char * GetDataName() const override;
};


