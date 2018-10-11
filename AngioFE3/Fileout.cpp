///////////////////////////////////////////////////////////////////////
// Fileout.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Fileout.h"
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
	logstream << "Time,Segments,Total Length,Vessels,Branches,Anastamoses,Active Tips" << endl;

#ifndef NDEBUG
	vessel_csv_stream.open("vessels.csv");
	vessel_csv_stream << "x0,y0,z0,x1,y1,z1,time " << endl;
#endif

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
	fprintf(feangio_state_stream, "%-64s,%-64s,%-64s,%-64s\n",
		"Timestep", "Growth Process","Branch Policy Update", "Update Stress");
}

//-----------------------------------------------------------------------------
Fileout::~Fileout()
{
    logstream.close();
	fclose(vessel_state_stream);
	fclose(feangio_state_stream);
}

//-----------------------------------------------------------------------------
void Fileout::printStatus(FEAngio& angio, double time)
{
	//just store the status in a csv
	//logstream << "Time,Segments,Total Length,Vessels,Branchs,Anastamoses,Active Tips" << endl;

	logstream << angio.GetFEModel()->GetTime().currentTime << ","<< getSegmentCount(angio) <<"," << getSegmentLength(angio,time)<< ",," << getBranchCount(angio) << ",," << getTipCount(angio) << std::endl;
}

void PrintSegment(vec3d r0, vec3d r1)
{
#ifndef NDEBUG
	std::cout << "r0: " << r0.x << "," << r0.y << "," << r0.z << std::endl;
	std::cout << "r1: " << r1.x << "," << r1.y << "," << r1.z << std::endl;
#endif
}

//-----------------------------------------------------------------------------
// Save microvessel position at the current time point
void Fileout::save_vessel_state(FEAngio& angio)
{
	unsigned int new_seg_count = 0;
	for(int i=0; i <angio.angio_elements.size();i++)
	{
		new_seg_count += angio.angio_elements[i]->recent_segments.size();
	}

	//write segcount and time
	fwrite(&new_seg_count, sizeof(unsigned int), 1, vessel_state_stream);
	auto ti = angio.GetFEModel()->GetTime();
	float rtime = static_cast<float>(ti.currentTime);
	fwrite(&rtime, sizeof(float), 1, vessel_state_stream) == sizeof(float);

	unsigned int crc_segcount = 0;
	for (int i = 0; i <angio.angio_elements.size(); i++)
	{
		for(int j=0; j <  angio.angio_elements[i]->recent_segments.size();j++)
		{
			crc_segcount++;
			Tip * back_tip = angio.angio_elements[i]->recent_segments[j]->back;
			Tip * front_tip = angio.angio_elements[i]->recent_segments[j]->front;
			vec3d r0 = angio.ReferenceCoordinates(back_tip);
			vec3d r1 = angio.ReferenceCoordinates(front_tip);

			PrintSegment(r0, r1);
			//consider checking fwrites return values
			float cur = static_cast<float>(r0.x);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r0.y);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r0.z);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);

			cur = static_cast<float>(r1.x);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r1.y);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r1.z);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);

			r1 = angio.ReferenceCoordinates(front_tip);
#ifndef NDEBUG
			vessel_csv_stream << r0.x << "," << r0.y << "," << r0.z << "," << r1.x << "," << r1.y << "," << r1.z << "," << rtime << std::endl;
#endif
		}
	}
	assert(crc_segcount == new_seg_count);
}

void Fileout::bulk_save_vessel_state(FEAngio& angio)
{
	unsigned int new_seg_count = 0;
	for (int i = 0; i <angio.angio_elements.size(); i++)
	{
		new_seg_count += angio.angio_elements[i]->grown_segments.size();
	}

	//write segcount and time
	fwrite(&new_seg_count, sizeof(unsigned int), 1, vessel_state_stream);
	auto ti = angio.GetFEModel()->GetTime();
	float rtime = static_cast<float>(ti.currentTime);
	fwrite(&rtime, sizeof(float), 1, vessel_state_stream) == sizeof(float);

	unsigned int crc_segcount = 0;
	for (int i = 0; i <angio.angio_elements.size(); i++)
	{
		for (int j = 0; j < angio.angio_elements[i]->grown_segments.size(); j++)
		{
			crc_segcount++;
			Tip * back_tip = angio.angio_elements[i]->grown_segments[j]->back;
			Tip * front_tip = angio.angio_elements[i]->grown_segments[j]->front;
			vec3d r0 = angio.ReferenceCoordinates(back_tip);
			vec3d r1 = angio.ReferenceCoordinates(front_tip);

			PrintSegment(r0, r1);
			//consider checking fwrites return values
			float cur = static_cast<float>(r0.x);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r0.y);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r0.z);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);

			cur = static_cast<float>(r1.x);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r1.y);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r1.z);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);

			r1 = angio.ReferenceCoordinates(front_tip);
#ifndef NDEBUG
			vessel_csv_stream << r0.x << "," << r0.y << "," << r0.z << "," << r1.x << "," << r1.y << "," << r1.z << "," << rtime << std::endl;
#endif
		}
	}
	assert(crc_segcount == new_seg_count);
}

//-----------------------------------------------------------------------------
// Save active points
void Fileout::save_active_tips(FEAngio& angio) const
{	
	
}
void Fileout::save_timeline(FEAngio& angio)
{
#ifndef NDEBUG

#endif
}

void Fileout::save_final_vessel_csv(FEAngio & angio)
{


}

void Fileout::save_feangio_stats(FEAngio& angio)
{
	FETimeInfo & ti = angio.m_fem->GetTime();
	const int STR_SIZE = 64;
	char grow_time[STR_SIZE];
	char update_branch_policy[STR_SIZE];
	char update_stress[STR_SIZE];
	angio.grow_timer.time_str(grow_time);
	angio.update_branch_policy_timestep_timer.time_str(update_branch_policy);
	angio.update_angio_stress_timer.time_str(update_stress);
	fprintf(feangio_state_stream, "%-12.7f,%-64s,%-64s,%-64s\n",
		ti.currentTime, grow_time, update_branch_policy, update_stress);
	fflush(feangio_state_stream);
}

void Fileout::save_winfiber(FEAngio& angio)
{

}

int Fileout::getBranchCount(FEAngio& angio)
{
	int count = 0;
	for(int i=0; i < angio.angio_elements.size();i++)
	{
		count += angio.angio_elements[i]->branch_count;
	}
	return count;
}

int Fileout::getSegmentCount(FEAngio& angio)
{
	int count = 0;
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		count += angio.angio_elements[i]->grown_segments.size();
	}
	return count;
}

double Fileout::getSegmentLength(FEAngio& angio, double time)
{
	FEMesh * mesh = angio.GetMesh();
	double len= 0;
#pragma omp parallel for
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		double temp = angio.angio_elements[i]->GetLengthAtTime(mesh, time);
#pragma omp critical
		len += temp;
	}
	return len;
}

int Fileout::getTipCount(FEAngio& angio)
{
	int count = 0;
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		for(auto iter = angio.angio_elements[i]->next_tips.begin(); iter != angio.angio_elements[i]->next_tips.end();++iter)
		{
			count += iter->second.size();
		}
	}
	return count;
}


