#include "MALARIC/Tools/Amplitude.H"

#include "MALARIC/Shower/Kernel.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace MALARIC;
using namespace ATOOLS;

Amplitude::Amplitude(Cluster_Amplitude *const a):
  m_mec(0), m_t(0.0), m_t0(0.0), p_ampl(a)
{
}

Amplitude::~Amplitude()
{
  for (const_iterator it(begin());
       it!=end();++it) delete *it;
}

void Amplitude::Add(Parton *const p)
{
  push_back(p);
}

void Amplitude::Remove(Parton *const p)
{
  if (back()!=p) Abort(); 
  pop_back();
  delete p;
}

ATOOLS::Cluster_Amplitude *Amplitude::GetAmplitude() const
{
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  ampl->CopyFrom(p_ampl,1);
  for (const_iterator it(begin());it!=end();++it)
    ampl->CreateLeg((*it)->Mom(),(*it)->Flav(),
		    ColorID((*it)->Col().m_i,(*it)->Col().m_j));
  return ampl;
}

namespace MALARIC {

  std::ostream &operator<<(std::ostream &s,const Amplitude &a)
  {
    Vec4D p;
    int c[4]={0,0,0};
    s<<"("<<&a<<"): t = "<<a.T()<<", t0 = "<<a.T0()
     <<", nlo = "<<ID(a.ClusterAmplitude()->NLO())
     <<", flag = "<<ID(a.ClusterAmplitude()->Flag())
     <<", mec = "<<a.MEC()
     <<" {\n  "<<a.Split()<<"\n";
    for (Amplitude::const_iterator
	   it(a.begin());it!=a.end();++it) {
      msg_Indent();
      p+=(*it)->Mom();
      ++c[(*it)->Col().m_i];
      --c[(*it)->Col().m_j];
      s<<**it<<"\n";
    }
    return s<<"  \\sum p = "<<p
	    <<", \\sum c = ("<<c[1]
	    <<","<<c[2]<<","<<c[3]<<")\n}";
  }

}
