#pragma once
#include <FECore/FEPlotData.h>

//! plots the stress from the stress policy
class FEPlotAngioStress : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotAngioStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	//! plot angio stress
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//! Plot the componet of stress from the matrix material
class FEPlotMatrixStress : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	//! plot matrix stress
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//! plot the visco component of the matrix material
class FEPlotMatrixViscoStress : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixViscoStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	//! plot matrix visco elastic stress
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//! Plot the elastic component of the matrix material
class FEPlotMatrixElasticStress : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixElasticStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	//! plot matrix elastic stress
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//! plot the material orientation of the matrix elastic material
class FEPlotMatrixElastic_m_Q : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixElastic_m_Q(FEModel* pfem) : FEDomainData(PLT_MAT3F, FMT_ITEM){}
	//! plot matrix material orientation
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//! plot the stress of the vessel material
class FEPlotVesselStress : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotVesselStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	//! plot vessel stress
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//! plot the weight of the vessel material
class FEPlotVesselWeight : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotVesselWeight(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	//! plot vessel weight
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//! plot the weight of the matrix material
class FEPlotMatrixWeight : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixWeight(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	//! plot matrix weight
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//! Plot the matrix material's tangent
class FEPlotMatrixTangent : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotMatrixTangent(FEModel* pfem) : FEDomainData(PLT_TENS4FS, FMT_ITEM){}
	//! plot matrix tangent
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//! plot the branches per element
class FEPlotBranches : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotBranches(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	//! plot the number of branches per element
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//! plot the anastamoses per element
class FEPlotAnastamoses : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotAnastamoses(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	//! plot the number of anastamoses per element
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//! plot the segment length per element
class FEPlotSegmentLength : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotSegmentLength(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	//! plot the segment length per element
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//! plot the segment length per element in the reference frame
class FEPlotRefSegmentLength : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotRefSegmentLength(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	//! plot the segment length per element in the reference frame
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//! plot the ecm density per element
class FEPlotAngioECMDensity : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotAngioECMDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	//! plot the ecm density per element
	bool Save(FEDomain& d, FEDataStream& str)  override;
};

//! plot the repulse value per element
class FEPlotAngioRepulseVal : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotAngioRepulseVal(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	//! plot the repulse value per element
	bool Save(FEDomain& d, FEDataStream& str)  override;
};

class FEPlotAngioSPA : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotAngioSPA(FEModel* pfem) : FEDomainData(PLT_MAT3F, FMT_ITEM) {}
	//! plot the semi principal axes per element
	bool Save(FEDomain&d, FEDataStream& str) override;
};

class FEPlotAngioFractionalAnisotropy : public FEDomainData
{
public:
	//! constructor
	explicit FEPlotAngioFractionalAnisotropy(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	//! plot the semi principal axes per element
	bool Save(FEDomain&d, FEDataStream& str) override;
};

//! plot the primary direction that vessels grow per element
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
