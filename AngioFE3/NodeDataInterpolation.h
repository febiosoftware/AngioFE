#pragma once
#include <FECore/FEMaterial.h>
#include "FEAngioMaterialPoint.h"
#include "VariableInterpolation.h"

class FEAngio;
class FEAngioMaterial;

//! the ability to set certain values on a per node basis
class NodeDataInterpolation :public FEMaterial
{
public:
	//! constructor
	NodeDataInterpolation(FEModel * pfem) :FEMaterial(pfem) { }
	//! returns the value from a material point
	virtual double & ValueReference(FEMaterialPoint * mp) = 0;
	//! the name to which the data is bound
	virtual const char * GetDataName() const = 0;
	//! do the initialization of values before the first step
	void PrepValues(FEAngio* angio, FEMesh* mesh, FEAngioMaterial* angio_mat);
protected:
	//! id of the nodeset that is effected
	int node_set_id = -1;
	//! map of node ids to the values that are associated with them
	std::unordered_map<int, double> node_id_to_values;
	//! 0 use shape functions to interpolate to gauss points, 1 average all nodes in element and set that value at all gauss points
	int interpolation_mode = 0;
	//! parameter list
	DECLARE_FECORE_CLASS();
};

//! applies all node data interpolation values
class NodeDataInterpolationManager :public FEMaterial
{
public:
	//! constructor
	NodeDataInterpolationManager(FEModel * pfem) :FEMaterial(pfem) { AddClassProperty(this, &node_data_interpolation_vals, "node_interpolation_value", FEProperty::Optional); }
	//! performs initialization
	bool Init() override { return FEMaterial::Init(); }
	//! do the interpolation for the given element
	void DoInterpolations(FEAngio * angio, FEMesh* mesh, FEAngioMaterial* angio_mat);
private:
	std::vector<NodeDataInterpolation*>	node_data_interpolation_vals;	//!< pointers to elastic materials
	//NodeDataInterpolation* node_data_interpolation_vals;
};

//! set ecm density values on a per node basis
class DensityValuesNodeDataInterpolation : public NodeDataInterpolation
{
public:
	//! constructor
	DensityValuesNodeDataInterpolation(FEModel * pfem) :NodeDataInterpolation(pfem) {}
	//! performs initialization
	bool Init() override { return FEMaterial::Init(); }
	//! returns the density values from a material point
	double & ValueReference(FEMaterialPoint * mp) override;
	//! the name to which the data is bound
	const char * GetDataName() const override;
};

//! set repulse values on a per node basis(replaces bouncy boundary condition)
class RepulseValuesNodeDataInterpolation : public NodeDataInterpolation
{
public:
	//! constructor
	RepulseValuesNodeDataInterpolation(FEModel * pfem) :NodeDataInterpolation(pfem) {}
	//! performs initialization
	bool Init() override { return FEMaterial::Init(); }
	//! returns the repulse value from a material point
	double & ValueReference(FEMaterialPoint * mp) override;
	//! the name to which the data is bound
	const char * GetDataName() const override;
};

