// former Flow.cxx of HepMC
#include "Flow.H"
#include "Parton.H"
//#include "Blob.H"


namespace APHYTOOLS {
  long Flow::qcd_counter = 600;
}

using namespace APHYTOOLS;

Flow::Flow( Parton * _owner ) : m_owner(_owner) { m_code[1]=0; m_code[2]=0; }

Flow::Flow( const Flow & inflow ) : 
  m_owner(inflow.m_owner)
{
  *this = inflow;
}

Flow::~Flow() { 
  m_code.clear(); 
}
	
/////////////
// Friends //
/////////////

std::ostream& APHYTOOLS::operator<<( std::ostream& ostr, const Flow& f ) {
  ostr << f.m_code.size();
  for ( std::map<int,int>::const_iterator i = f.m_code.begin();
	i != f.m_code.end(); ++i ) {
    ostr << " " << (*i).first << " " << (*i).second;
  }
  return ostr;
}


















