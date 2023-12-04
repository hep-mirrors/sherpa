#include "ALARIC++/Tools/Amplitude.H"

#include "ALARIC++/Shower/Kernel.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace ALARIC;
using namespace ATOOLS;

Amplitude::Amplitude(Cluster_Amplitude *const a,
		     std::vector<Amplitude*> *const all):
  m_t(0.0), m_t0(0.0), p_ampl(a), p_all(all)
{
}

Amplitude::~Amplitude()
{
  for (const_iterator it(begin());
       it!=end();++it) delete *it;
}

void Amplitude::Reduce()
{
  for (iterator it(begin());it!=end();++it) {
    (*it)->SetOut(0,NULL);  
    (*it)->SetOut(1,NULL);
  }
  bool kill(false);
  for (Amplitude_Vector::iterator
	 ait(p_all->begin());ait!=p_all->end();) {
    if (!kill) ++ait;
    else { delete *ait; ait=p_all->erase(ait); }
    if (*(ait-1)==this) kill=true;
  }
}

void Amplitude::Add(Parton *const p)
{
  push_back(p);
  if (front()->Out(0)==NULL) return;
  Parton *cp(new Parton(front()->Out(0)->Ampl(),
			p->Flav(),p->Mom(),
			p->Col(),p->Hel()));
  cp->SetId(cp->Counter());
  p->SetOut(0,cp);
  cp->SetIn(p);
  front()->Out(0)->Ampl()->Add(cp);
}

void Amplitude::Remove(Parton *const p)
{
  if (p->Out(0)) p->Out(0)->Ampl()->Remove(p->Out(0));
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

namespace ALARIC {

  std::ostream &operator<<(std::ostream &s,const Amplitude &a)
  {
    Vec4D p;
    int c[4]={0,0,0};
    s<<"("<<&a<<"): t = "<<a.T()<<", t0 = "<<a.T0()
     <<", nlo = "<<ID(a.ClusterAmplitude()->NLO())
     <<", flag = "<<ID(a.ClusterAmplitude()->Flag())
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
