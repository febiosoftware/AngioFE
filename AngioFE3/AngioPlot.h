#pragma once
#include <FECore/FEPlotData.h>

//base class for computationally intensive plot variables that cannot be cached
//-----------------------------------------------------------------------------
class CIFEDomainData : public FEDomainData
{
public:
	//! constructor
	explicit CIFEDomainData(Var_Type t, Storage_Fmt s) : FEDomainData(t, s) {}
	//! update any cached data as this fucntion
	virtual void Update(FEModel* model) = 0;

};

//-----------------------------------------------------------------------------
class FEPlotAngioStress : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotAngioStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	//! plot angio stress
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotMatrixStress : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	//! plot matrix stress
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotMatrixViscoStress : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixViscoStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	//! plot matrix visco elastic stress
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotMatrixElasticStress : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixElasticStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	//! plot matrix elastic stress
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotMatrixElastic_m_Q : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixElastic_m_Q(FEModel* pfem) : FEDomainData(PLT_MAT3F, FMT_ITEM){}
	//! plot matrix material orientation
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotVesselStress : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotVesselStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	//! plot vessel stress
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotVesselWeight : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotVesselWeight(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	//! plot vessel weight
	bool Save(FEDomain& d, FEDataStream& str) override;
};
//-----------------------------------------------------------------------------
class FEPlotMatrixWeight : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixWeight(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	//! plot matrix weight
	bool Save(FEDomain& d, FEDataStream& str) override;
};
//-----------------------------------------------------------------------------
class FEPlotMatrixTangent : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixTangent(FEModel* pfem) : FEDomainData(PLT_TENS4FS, FMT_ITEM){}
	//! plot matrix tangent
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotBranches : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotBranches(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	//! plot the number of branches per element
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotAnastamoses : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotAnastamoses(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	//! plot the number of anastamoses per element
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotSegmentLength : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotSegmentLength(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	//! plot the segment length per element
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotRefSegmentLength : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotRefSegmentLength(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	//! plot the segment length per element in the reference frame
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioECMDensity : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotAngioECMDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	//! plot the ecm density per element
	bool Save(FEDomain& d, FEDataStream& str)  override;
};

//-----------------------------------------------------------------------------
class FEPlotPrimaryVesselDirection : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotPrimaryVesselDirection(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM) {}
	//! plot the primary direction that vessels grow per element
	bool Save(FEDomain& d, FEDataStream& str)  override;
private:
	//maps element id to segment direction in the reference configuration
	std::unordered_map<int, vec3d> accumulated_seg_direction;
};
