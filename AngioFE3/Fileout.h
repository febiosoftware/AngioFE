///////////////////////////////////////////////////////////////////////
// Fileout.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The FILEOUT class writes the output files that contain the results
// of the simulation.
///////////////////////////////////////////////////////////////////////

#pragma once

#include "StdAfx.h"
//#include "zlib.h"

class FEAngio;

//! class that contains the streams for doing output
class Fileout
{
public:
	//! constructor
	Fileout(FEAngio& angio);
	virtual ~Fileout();
	//! print the status to the console
	void printStatus(FEAngio& angio, double time);
	//! save the vessel state
	void save_vessel_state(FEAngio& angio);
	//! saves all segments not just recent ones, should only be used in initialization
	void bulk_save_vessel_state(FEAngio& angio);
	//! save the active tips
	void save_active_tips(FEAngio& angio) const;
	//! save the timeline of branch points
	void save_timeline(FEAngio& angio);

	//! timing statistics
	void save_feangio_stats(FEAngio& angio);
	//! save the final vascular network as a csv
	static void save_final_vessel_csv(FEAngio & angio);

private:
	int getBranchCount(FEAngio& angio);
	// get per angio material
	std::vector <int> getBranchCount_pm(FEAngio& angio);
	int getSegmentCount(FEAngio& angio);
	// get per angio material
	std::vector <int> getSegmentCount_pm(FEAngio& angio);
	double getSegmentLength(FEAngio& angio, double time);
	// get per angio material
	std::vector <double> getSegmentLength_pm(FEAngio& angio, double time);
	int getTipCount(FEAngio& angio);
	// get per angio material
	std::vector <int> getTipCount_pm(FEAngio& angio);
	std::ofstream logstream;

	FILE*	m_stream4 = 0;	// active tips
	FILE*  vessel_state_stream=0;
	FILE*  feangio_state_stream=0;
#ifndef NDEBUG
	std::ofstream vessel_csv_stream;
#endif
};
