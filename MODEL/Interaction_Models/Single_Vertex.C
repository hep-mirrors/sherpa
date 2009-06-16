#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

// Constructors and Destructors
Single_Vertex::Single_Vertex()
{ t = 0; nleg=3; cpl.resize(4);}

Single_Vertex::Single_Vertex(const Single_Vertex& v): 
  t(0)
{ 
  *this=v; 
}

Single_Vertex::~Single_Vertex()
{
  for (size_t i(0);i<Lorentz.size();++i) delete Lorentz[i];
}

Complex Single_Vertex::Coupling(size_t i) const
{
  return cpl[i].Value();
}
 
// Operators
Single_Vertex& Single_Vertex::operator=(const Single_Vertex& v) 
{
  for (size_t i(0);i<Lorentz.size();++i) delete Lorentz[i];
  Lorentz=std::vector<Lorentz_Function*>();
    
      if (this!=&v) {
	for (short int i=0;i<4;i++) in[i]  = v.in[i];
	cpl.clear();
	for (size_t j=0;j<v.cpl.size();j++) cpl.push_back(v.cpl[j]);
	
	nleg = v.nleg;
	Str  = v.Str;
	on   = v.on;
	t=v.t;
	Color.resize(v.Color.size());
	for (size_t i(0);i<Color.size();++i)
	  Color[i] = v.Color[i];
	Lorentz.resize(v.Lorentz.size());
	for (size_t i(0);i<Lorentz.size();++i)
	  Lorentz[i] = v.Lorentz[i]->GetCopy();
      }
      return *this;
    }


 
const bool Single_Vertex::operator==(const Single_Vertex& probe) 
{
  switch (nleg) // different checks for 3-leg and 4-leg vertices
  {
    case 4: return (probe.nleg == 4) &&
       	           (in[0] == probe.in[0]) &&
                   (in[1] == probe.in[1]) &&
                   (in[2] == probe.in[2]) &&
                   (in[3] == probe.in[3]);
    case 3:  return (probe.nleg == 3) &&
	            (in[0] == probe.in[0]) &&
                    (in[1] == probe.in[1]) &&
                    (in[2] == probe.in[2]);
    default: return 0; 
  }
}


ostream& MODEL::operator<<(ostream& s, const Single_Vertex& sv)
{
  return s<<'('<<sv.in[0]<<','<<sv.in[1]<<','<<sv.in[2]<<','<<sv.in[3]
          <<") with cpl["<<sv.Coupling(0)<<','<<sv.Coupling(1)<<','<<sv.Coupling(2)<<','<<sv.Coupling(3)<<']'
	  <<" is "<<((sv.on) ? "on" : "off");
}


// MPI stuff
void MODEL::Single_Vertex2MPI(const Single_Vertex * v , MPI_Single_Vertex & mpi_v) {
  
    
  if (!v) {
    for (int i=0;i<4;++i)
      mpi_v.m_fl[i] = 0;
    
    return; 
  }
  /*
  Lorentz_Function2MPI(v->Lorentz,mpi_v.m_lf);
  Color_Function2MPI(v->Color,mpi_v.m_cf);
  */

  for (int i=0;i<4;++i)
    mpi_v.m_fl[i] = int(v->in[i]);

  /*  
  for (int i=0;i<7;i+=2) {
    mpi_v.m_cpl[i]   = real(v->cpl[i/2]);
    mpi_v.m_cpl[i+1] = imag(v->cpl[i/2]);
  }
  */
}


Single_Vertex * MODEL::MPI2Single_Vertex(const MPI_Single_Vertex & mpi_v ) {

  Single_Vertex * v ;
  
  v = new Single_Vertex();
  
  for (int i=0;i<4;++i) {
    v->in[i] = Flavour((kf_code)(abs(mpi_v.m_fl[i])));
    if (mpi_v.m_fl[i]<0) v->in[i]=v->in[i].Bar();
  }

  /*
  v->Lorentz = MODEL::MPI2Lorentz_Function(mpi_v.m_lf);
  v->Color   = MODEL::MPI2Color_Function(mpi_v.m_cf);
  
  for (int i=0;i<7;i+=2) {
    (v->cpl[i/2]) = Complex(mpi_v.m_cpl[i],mpi_v.m_cpl[i+1]);
  }
  */
  
  return v;

}

std::ostream & MODEL::operator<<(std::ostream & s, const MPI_Single_Vertex & sv) {
  s<<sv.m_fl[0]<<","<<sv.m_fl[1]<<","<<sv.m_fl[2]<<","<<sv.m_fl[3];
  return s;
}
