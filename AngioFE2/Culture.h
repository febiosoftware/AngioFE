#pragma once

#include <list>
#include "BC.h"
#include "Segment.h"
#include <vector>
#include <fstream>
#include <random>
#include <set>
#include <algorithm>

//-----------------------------------------------------------------------------
class FEAngio;
class SimulationTime;
class GridPoint;
class FEAngioMaterial;
class BouncyBC;

//-----------------------------------------------------------------------------
typedef list<Segment> SegmentList;
typedef list<Segment>::iterator SegIter;
typedef list<Segment>::const_iterator ConstSegIter;
typedef list<Segment::TIP*>	SegmentTipList;
typedef list<Segment::TIP*>::iterator TipIter;
typedef list<Segment::TIP*>::const_iterator ConstTipIter;

const int MAXPARAMSIZE = 256;


//contains the material parameters and the accessor fucntions when the variables are made dependent on time
class CultureParameters
{
public:
	CultureParameters();                                                         

	//not currently exposed
	//These parameters govern the amount of growth at each time step for the microvessels
	//          g(t) = y0 + a/(1+e^-(t-x0)/b)
	double	m_y0 = -19.1278;   // Bottom of sigmoid curve
	double	m_culture_a = 1900.0;    // Distance between top and bottom of the curve (a + y0 = top of curve)
	double	m_x0 = 4.9474;   // Time point at which t is halfway between top & bottom
	double	m_culture_b = 1.4549;    // Steepness of the curve

	//the initial material parameters a,b,N
	double sprout_s_mag = 3.72e-12;
	double sprout_s_range = 0.004;
	double sprout_s_width = 2.0;

	//dont expose
	double	m_initial_vessel_length = 0.0;	// initial vessel length
	double	m_initial_branch_probability = 0.0;     // Probability that initial segments will branch
	bool m_bzfibflat = false;//flatten the fibers in the z direction

	//parameters adjustable by the user for vessel network morphology
	// Length adjuster, scales the length of new segments
	// allows fine tuning of total vascular length produced by the simulation.
	double m_length_adjustment = 1.0;
	// If the shortest distance vector between a segment and any other segment is less than this value
	//TODO: I think making this a percentage of the growth length makes more sense.
	double m_anastomosis_distance = 25;
	bool	m_branch= false; //whether the vessels are allowed to branch
	bool	m_anastomosis = false; //whether the vessels are allowed to fuse together
	// vessel_width - Diameter of microvessels (Default: 7 um)
	double m_vessel_width = 7.0;
	char m_boundary_condition_type[MAXPARAMSIZE];//currently s for stop, or b for bouncy 

	//parameters for ECM density/alignment
	int m_matrix_condition = 0; // flag indicating how the collagen fibers are oriented initially ( 0 = random classic mode, 1 multimaterial mode,  3 = along local element direction)
	int ecm_control = 0; //flag indicating how the ecm density and anisotropy are initialized: 0 constant mode, 1 specified mode, 2 no overwrite 
	double m_matrix_density = 3.0;//mg/ml 3.0 default

	//symmety
	vec3d m_symmetry_plane = vec3d(0.0, 0.0, 0.0); //symmetry plane

	//uncategoriezed
	int m_composite_material = 1;//whether or not the material is a composite material
	
	//TODO: doesn't appear to be the same as sprout_s_mag
	double m_sprout_force = 1.0;

	int m_number_fragments = 0;//number of fragments to seed in the material

	//TODO: not currently exposed
	int tries_per_segment = 10;// number of times the algorithm will retry finding an initial segment
	//when all weights of the elemetns are equal 1 will work but this needs increased in porportion to the maximun difference in element weights

	vec3d	vessel_orient_weights = vec3d(0.9090909090, 0.9090909090, 0.0);			// W.x = Weight for collagen orientation
	//W.y = Weight for previous vessel direction
	//W.z is unused
	// Microvessel growth rate is modeled using the Boltzmann sigmoid curve 

	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);

	double min_segment_length = 0.1;

	int fragment_seeder = 1;//0 classic fragment seeder, 1 multidomian fragment seeder, 2 by volume fragment seeder, 3 from file

	int branching_scheme = 1;

	char vessel_file[MAXPARAMSIZE];//only used if seeding from file

	int angio_boundary_type = 1;// 0 same as non angio materials, 1 pass through
	int angio_boundary_groups = 1;//each bit in this parameter fragments can only travel between the groups they are in


	double density_gradient_threshold = 0.01;//set this to be higher to turn off the direction change on encountering different densities

	double average_length_to_branch_point = 240.0;//sets the average length before a branch is encountered
	double std_deviation = 5.0;
	//consider allowing the user to specify the distribution and any paramters associated with the distribution that is used

	const bool io = true;
	friend class FEParamContainer;
	friend class FEAngioMaterial;

	double GetBranchProbility(double dt)const{ return m_branch_chance * dt; }
	double GetWeightInterpolation(double dt) const { return m_weight_interpolation * dt; }
private:
	double	m_branch_chance = 0.1;    // Probability of forming a new branch
	double m_weight_interpolation = 0.4; //used to control the amount of affect of the collagen direction vs the previous direction for an approximation over 1t

};
//used for generating initial segments
struct SegGenItem
{
	bool operator<(SegGenItem const & item)const{ return weight < item.weight; }
	double weight;
	FEDomain * domain;
	int ielement;
};
class Culture;
class FragmentSeeder
{
public:
	virtual bool SeedFragments(SimulationTime& time, Culture * culture)=0;
	FragmentSeeder(CultureParameters * cp, FEAngio & angio) : m_angio(angio), culture_params(cp) {}
	virtual ~FragmentSeeder(){}
protected:
	FEAngio & m_angio;
	CultureParameters * culture_params;
	std::vector<FEDomain *> domains;
	virtual bool createInitFrag(Segment& Seg, SegGenItem & item, Culture * culture);
};

class ClassicFragmentSeeder : public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	ClassicFragmentSeeder(CultureParameters * cp, FEAngio & angio) : FragmentSeeder(cp, angio) {}
private:
	// Seed an initial fragment within the grid
	bool createInitFrag(Segment& Seg);
};
class MultiDomainFragmentSeeder : public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	MultiDomainFragmentSeeder(CultureParameters * cp, FEAngio & angio) : FragmentSeeder(cp, angio) {}
private:
};

class MDByVolumeFragmentSeeder : public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	MDByVolumeFragmentSeeder(CultureParameters * cp, FEAngio & angio) : FragmentSeeder(cp, angio) {}
private:
};

class MDAngVessFileFragmentSeeder :public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	MDAngVessFileFragmentSeeder(CultureParameters * cp, FEAngio & angio) : FragmentSeeder(cp, angio) {}
private:
	std::ifstream infile;
};
//this class encapsulates the branching behavior for a given culture(instance)
//the static methods/members synchronize the instances where needed, these should only be called once from the model
class FragmentBranching
{
public:
	class BranchPoint
	{
	public:
		BranchPoint(double emt, double ept, Segment * p, double pctp, int prior, FragmentBranching * fb) :emerge_time(emt),
			epoch_time(ept), parent(p), percent_of_parent(pctp), priority(prior), brancher(fb)
		{
			assert(emerge_time >= -1.0);
			assert(epoch_time >= -1.0);
			assert(parent != nullptr);
			assert(brancher != nullptr);
		}
		~BranchPoint(){}

		double emerge_time;//time at which this branch begins to grow
		double epoch_time;//time when this branch point was created
		Segment * parent;
		double percent_of_parent;//0-1.0 determines where on the parent the branch will sprout from 0 is 100% contribution from tip(0)
		int priority;//it there is a tie in time this will break it consistently needs to be athe same for the brnachpoints between runs with equivalent paramters and unique among branch points
		FragmentBranching * brancher;//used to get the rng needed for this segment

		bool branch; //used in the timeline to tell if the generated point is generating a branch

		//include utility if the other relational operators are needed
		//should allos the set to be iterated over from low to high times
		bool operator< (const BranchPoint& rhs) const
		{
			if (emerge_time < rhs.emerge_time)
				return true;
			else if (emerge_time > rhs.emerge_time)
				return false;
			else
				return priority < rhs.priority;
		}
		bool operator== (const BranchPoint & rhs) const
		{
			return emerge_time == rhs.emerge_time && parent == rhs.parent && percent_of_parent == rhs.percent_of_parent && priority == rhs.priority;
		}
		
	};
	//allows Branch points to be sorted w/i a set by epoch time
	struct BranchPointEpochCompare
	{
		bool operator()(const BranchPoint & lhs, const BranchPoint & rhs) const
		{
			if (lhs.epoch_time < rhs.epoch_time)
				return true;
			else if (lhs.epoch_time > rhs.epoch_time)
				return false;
			else
				return lhs.priority < rhs.priority;
		}
	};
	//allows Branch points to be sorted w/i set sorted by the lowest of the times
	struct BranchPointTimeFloorCompare
	{
		bool operator()(const BranchPoint & lhs, const BranchPoint & rhs) const
		{
			return (std::min(lhs.epoch_time, lhs.emerge_time) < std::min(rhs.epoch_time, rhs.emerge_time));
		}
	};


	FragmentBranching(Culture *cp)
	{
		culture = cp;
		fragment_branchers.push_back(this);
	}
	virtual ~FragmentBranching()
	{
		assert(std::find(fragment_branchers.begin(), fragment_branchers.end(), this) != fragment_branchers.end());
		fragment_branchers.erase(std::find(fragment_branchers.begin(), fragment_branchers.end(),this));
	}
	//Grows a segment from a branchpoint
	virtual void GrowSegment(std::set<BranchPoint>::iterator bp) = 0;


	void GrowSegments();
	//Grows a segment from an active tip
	//this function should be time independent
	//thsi grows a segment from an active tip and sets the BranchPoints needed to generate any needed rng
	virtual void GrowSegment(Segment::TIP * tip, double starttime, double grow_time) = 0;

	virtual void UpdateSegmentBranchDistance(std::set<BranchPoint>::iterator bp) = 0;

	//modify the length to brnach if needed and set the time of birth of the segment
	virtual void PostProcess(Segment & seg)=0;

	virtual double GetLengthToBranch() = 0;

	//Grow is a synchronized grow operation for all FragmentBranchers
	static void Grow();

	//call once to have this class record a timeline of the branch points it uses
	static void OutputTimeline(){ create_timeline = true; }

	friend class Fileout;
protected:
	Culture * culture;
	double current_time;
	static std::vector<FragmentBranching *> fragment_branchers;//consider making this private
	static std::set<BranchPoint> branch_points;//used for creating the branches
	static std::set<BranchPoint, BranchPointEpochCompare> parentgen;//used to generate the next branchpoint for the parent segment
	static bool create_timeline;
	static std::multiset<BranchPoint, FragmentBranching::BranchPointTimeFloorCompare> timeline;
};
//fragments branch as they grow
class ForwardFragmentBranching :public FragmentBranching
{
public:
	ForwardFragmentBranching(Culture * cp) : FragmentBranching(cp)
	{
		//in the future have away to set these as paramters/ have classes dedicated to statistical distributions
		length_to_branch_point = normal_distribution<double>(40.0, 2.0);
	}
	void GrowSegment(std::set<BranchPoint>::iterator bp) override;
	void GrowSegment(Segment::TIP * tip, double starttime, double grow_time) override;
	void UpdateSegmentBranchDistance(std::set<BranchPoint>::iterator bp) override;
	void PostProcess(Segment & seg) override;
	double GetLengthToBranch() override;
private:
	normal_distribution<double> length_to_branch_point;
};
//fragments determine the branch points as they grow but there is a time delay to when they start growing
class PsuedoDeferedFragmentBranching :public FragmentBranching
{
	PsuedoDeferedFragmentBranching(Culture * cp) : FragmentBranching(cp){}
	void GrowSegment(std::set<BranchPoint>::iterator bp) override;
	void GrowSegment(Segment::TIP * tip, double starttime, double grow_time) override;
	void UpdateSegmentBranchDistance(std::set<BranchPoint>::iterator bp) override;
};
//for testing and replacing existing no branch functionality
class NoFragmentBranching : public FragmentBranching
{
public:
	NoFragmentBranching(Culture * cp) : FragmentBranching(cp){}
	void GrowSegment(std::set<BranchPoint>::iterator bp) override{}
	void GrowSegment(Segment::TIP * tip, double starttime, double grow_time) override;
	void UpdateSegmentBranchDistance(std::set<BranchPoint>::iterator bp) override{}
	void PostProcess(Segment & seg) override {}
	double GetLengthToBranch() override{ return 1000000.0; }
};
//the interface for all operations that modify the growth direction
//some examples of this include deflection due to gradient and change in direction for anastamosis
//in the future this can be used to change the direction based on vegf concentrations
class GrowDirectionModifier
{
public:
	GrowDirectionModifier(Culture * c, int p): culture(c), priority(p){}
	virtual ~GrowDirectionModifier(){}
	virtual vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) = 0;
	//used to sort these by priority
	int Priority() const { return priority; }
protected:
	Culture * culture;
	int priority;

};
//will ignore the previous direction and generate the direction a segmetn should grow based on collagen direction
class DefaultGrowDirectionModifier : public GrowDirectionModifier
{
public:
	DefaultGrowDirectionModifier(Culture * c, int p) : GrowDirectionModifier(c, p){}
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;
};

//this class changes the grow direction if the segment is a branch
class BranchGrowDirectionModifier : public GrowDirectionModifier
{
public:
	BranchGrowDirectionModifier(Culture * c, int p) : GrowDirectionModifier(c, p){}
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;
};

//the class modifies the grow dierction if the gradeint is above a given threshold
class GradientGrowDirectionModifier : public GrowDirectionModifier
{
public:
	GradientGrowDirectionModifier(Culture * c, int p) : GrowDirectionModifier(c, p){}
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;
};
//modifies the direction a segment grows based on its proximity to other segments
class AnastamosisGrowDirectionModifier : public GrowDirectionModifier
{
public:
	AnastamosisGrowDirectionModifier(Culture * c, int p) : GrowDirectionModifier(c, p){}
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;
};

//-----------------------------------------------------------------------------
// The CULTURE class contains all the functions that describe how the 
// SEGMENT class is to grow and orient itself. These functions are the 
// rule-set that arrange the line segments into vascular networks, 
// mimicking angiogenesis. 
class Culture  
{
public:
	Culture(FEAngio& angio, FEAngioMaterial * matl, CultureParameters * cp);
	virtual ~Culture();

	// initialize
	bool Init();

	// Reposition the vessels based on the FE solution
	void Update();

	// get the total vessel length
	double TotalVesselLength() const { return m_total_length; }
	
	// Find the density-based length scale factor at a point of the grid
	double FindDensityScale(const GridPoint& pt) const;
	
	// Add a segment to the culture (without checking BCs).
	// This assumes that the segment is valid.
	void AddSegment(Segment& seg);

	// Add a new segment to the culture.
	// This will apply BCs to the new segment and may result in 
	// the addition of several new segments. 
	//forcing a segment will ignore checks for minimum segment length
	void AddNewSegment(Segment& seg);

	//returns the segments that were added since the last call to AddNewSegment
	const std::vector<Segment *> & RecentSegments() const { return recents; }

	// return the number of segments
	int Segments() const { return m_nsegs; }

	// get the segment list
	const SegmentList& GetSegmentList() const { return m_frag; }

	// return the active segment list
	const SegmentTipList& GetActiveTipList() const
	{
		return m_active_tips;
	}
	//returns the active tips sorted by x position, the sorting key is arbitrary but should give consistent results between runs
	//consider sorting by distance from the origin instead
	const SegmentTipList & GetActiveSortedTipList()
	{
		m_active_tips.sort([](Segment::TIP * t0, Segment::TIP * t1){return t0->pos().x < t1->pos().x; });
		return m_active_tips;
	}

	// get the total number of active tips
	int ActiveTips() const { return static_cast<int>(m_active_tips.size()); }


	void ChangeBC(FEAngio & angio, int bcset);

	void CreateBranchingForce(Segment& seg);

	// Find the active tips
	void FindActiveTips();

	// Grow a segment
	Segment GrowSegment(Segment::TIP& it,double starttime, double grow_time, bool branch = false, bool bnew_vessel = false);

	//returns the segment length given the time parameters this will be further adjusted by ECM density
	double SegmentLength(double starttime, double grow_time) const;

	// create a branch
	void BranchSegment(Segment::TIP& it, double starttime, double grow_time);


private:	
	// fuse segments (i.e. anastomosis)
	void FuseVessels();

	// Create a new segment connecting two existing segments that are fusing through anastomosis
	static Segment ConnectSegment(Segment& it, Segment& it2, int k, int kk);

	


	int					m_nsegs;				// Counter that stores in current number of Segments within the simulation domain
	double				m_total_length;			// Total vascular length within the domain (sum of the length of all Segments) (in um)
	list<Segment>		m_frag;					// vessel fragments
	list<Segment::TIP*> m_active_tips;			// list of active tips
	CultureParameters * m_cultParams;
	BC *		bc;
	FEAngio&	m_angio;
	
	FragmentSeeder * fseeder = nullptr;
	
	std::vector<Segment *> recents;//used to hold the segments added by the most recent call to AddNewSegment these segments will be ordered first to last

	std::vector<GrowDirectionModifier *> grow_direction_modifiers;//these will be applied in order to determine the direction a tip should grow

public:
	FragmentBranching * fbrancher = nullptr;
	double	m_vess_length;	// new segment length
	int		m_num_vessel;   // Counter that indicates the next vessel ID number of new Segments
	int		m_num_branches;		// Counter indicating the number of branches formed during the simulation
	FEAngioMaterial * m_pmat = nullptr;
	
	int m_num_anastom;	// Counter indicating the number of anastomoses formed during the simulation
	int m_num_zdead;

	

	friend class BouncyBC;
	friend class BC;
	friend class StopBC;
	friend class MBC;
	friend class PassThroughMBC;
};
