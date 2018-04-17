///////////////////////////////////////////////////////////////////////
// Fileout.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Fileout.h"
#include "Segment.h"
#include "FEAngio.h"
#include <iostream>
#include <FECore/FEModel.h>
#include "FEAngioMaterial.h"
#include <unordered_set>


using namespace std;

//-----------------------------------------------------------------------------
Fileout::Fileout(FEAngio& angio)
{
    logstream.open("out_log.csv");
	//write the headers
	logstream << "Time,Material,Segments,Total Length,Vessels,Branch Points,Anastamoses,Active Tips,Sprouts" << endl;

	vessel_state_stream = fopen("out_vess_state.ang2" , "wb");//check the parameters consider setting the compression level
	unsigned int magic = 0xfdb97531;
	unsigned int version = 1;
	unsigned int num_bitmasks = angio.m_fem->Materials()/32 + 1;
	fwrite(&magic, sizeof(unsigned int), 1, vessel_state_stream);
	fwrite(&version, sizeof(unsigned int), 1, vessel_state_stream);
	fwrite(&num_bitmasks, sizeof(unsigned int), 1, vessel_state_stream);
	for(int j=0; j < num_bitmasks;j++)
	{
		unsigned int c_bitmask = 0;
		unsigned int place_holder = 1;
		for (int i = 0; i < 32; i++)
		{
			int index = i + j * 32;
			if (index == angio.m_fem->Materials())
				break;
			FEMaterial* mat = angio.m_fem->GetMaterial(index);
			FEAngioMaterial * angio_mat = angio.GetAngioComponent(mat);
			if(angio_mat)
			{
				c_bitmask |= place_holder;
			}
			place_holder = place_holder << 1;
		}
		fwrite(&c_bitmask, sizeof(unsigned int), 1, vessel_state_stream);
	}

	m_stream4 = fopen("out_active_tips.csv", "wt");		// active tips
	fprintf(m_stream4, "%-5s,%-12s,%-12s,%-12s\n", "State", "X", "Y", "Z");

	feangio_state_stream = fopen("angio_stats.csv", "wt");
	fprintf(feangio_state_stream, "%-64s,%-64s,%-64s,%-64s,%-64s,%-64s,%-64s,%-64s\n",
		"Timestep", "Grow Segments","Update Sprout Stress Scaling",
		"Update Sprout Stress", "Adjust Mesh Stiffness", "GDM Update",
		"ECM Update", "Material Update");
}

//-----------------------------------------------------------------------------
Fileout::~Fileout()
{
    logstream.close();
	fclose(vessel_state_stream);
	fclose(feangio_state_stream);
}

//-----------------------------------------------------------------------------
void Fileout::printStatus(FEAngio& angio)
{
	//now updated to a form which can be easily consumed by excel
	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{

	}
}

//-----------------------------------------------------------------------------
// Save microvessel position at the current time point
void Fileout::save_vessel_state(FEAngio& angio)
{
	

}

//-----------------------------------------------------------------------------
// Save active points
void Fileout::save_active_tips(FEAngio& angio) const
{	
	
}
void Fileout::save_timeline(FEAngio& angio)
{
#ifndef NDEBUG
	auto iter = FragmentBranching::timeline.begin();

	//consider how these shoudl be named to avoid collisions
	FILE * timeline_file = fopen("timeline.csv", "wt");
	assert(timeline_file);
	fprintf(timeline_file, "emerge_time,epoch_time,percent_of_parent,priority,callsite,branch\n");

	while (iter != FragmentBranching::timeline.end())
	{
		fprintf(timeline_file, "%-12.7f,%-12.7f,%-12.7f,%d\n", iter->emerge_time, iter->epoch_time, iter->percent_of_parent, iter->priority);
		++iter;
	}
	fclose(timeline_file);
#endif
}

void Fileout::save_final_vessel_csv(FEAngio & angio)
{


}

void Fileout::save_feangio_stats(FEAngio& angio)
{
	const int STR_SIZE = 64;
	char grow_time[STR_SIZE];
	char update_stress_scaling_time[STR_SIZE];
	char update_stress_time[STR_SIZE];
	char mesh_stiffnesss_time[STR_SIZE];
	char gdm_update_time[STR_SIZE];
	char ecm_update_time[STR_SIZE];
	char material_update_time[STR_SIZE];
	angio.grow_timer.time_str(grow_time);
	angio.update_sprout_stress_scaling_timer.time_str(update_stress_scaling_time);
	angio.mesh_stiffness_timer.time_str(mesh_stiffnesss_time);
	angio.update_gdms_timer.time_str(gdm_update_time);
	angio.update_ecm_timer.time_str(ecm_update_time);
	angio.material_update_timer.time_str(material_update_time);
	fprintf(feangio_state_stream, "%-12.7f,%-64s,%-64s,%-64s,%-64s,%-64s,%-64s,%-64s\n",
		angio.m_time.t, grow_time, update_stress_scaling_time, update_stress_time,
		mesh_stiffnesss_time, gdm_update_time, ecm_update_time, material_update_time);
	fflush(feangio_state_stream);
}

void Fileout::save_winfiber(FEAngio& angio)
{

}


