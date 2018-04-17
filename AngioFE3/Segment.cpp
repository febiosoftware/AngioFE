#include "StdAfx.h"
#include "Segment.h"
#include <FECore/vec3d.h>


//-----------------------------------------------------------------------------
Segment::TIP::TIP()
{
	bactive = false;
	nseed = -1;
	nvessel = -1;
	parent = nullptr;
}

//-----------------------------------------------------------------------------
Segment::Segment()
{
	// initialize data
	m_length = 0;                                        
	m_TofBirth = 0;
	
	// initialize IDs
	m_nid = 0;
	m_nseed   = -1;
    m_nvessel = -1;
    
	// initialize status flags
	m_nflag = 0;

	tip(0).parent = this;
	tip(1).parent = this;
}
Segment::Segment(const Segment &obj)
{
	//can i just do what ever would have been done before for the rest of the parameters
	memcpy(this, &obj, sizeof(Segment)); //-V598
	//update the parent pointers
	tip(0).parent = this;
	tip(1).parent = this;
}
void Segment::operator = (Segment& seg)
{
	//can i just do what ever would have been done before for the rest of the parameters
	memcpy(this, &seg, sizeof(Segment)); //-V598
										 //update the parent pointers
	tip(0).parent = this;
	tip(1).parent = this;
}

//-----------------------------------------------------------------------------
Segment::~Segment()
{

}
