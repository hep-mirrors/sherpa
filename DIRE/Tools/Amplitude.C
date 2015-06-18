#include "DIRE/Tools/Amplitude.H"

#include "DIRE/Shower/Kernel.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace DIRE;
using namespace ATOOLS;

Amplitude::Amplitude(Cluster_Amplitude *const a):
  m_t(0.0), m_t0(0.0), p_ampl(a)
{
}

Amplitude::~Amplitude()
{
  for (const_iterator it(begin());
       it!=end();++it) delete *it;
}

int Amplitude::Construct()
{
  if (m_s.p_c==NULL) return 1;
  for (iterator it(begin());it!=end();++it)
    (*it)->Out(0)->SetMom((*it)->Mom());
  Splitting cs(m_s);
  cs.p_n=m_s.p_c->Out(1);
  cs.p_c=m_s.p_c->Out(0);
  cs.p_s=m_s.p_s->Out(0);
  cs.p_sk->LF()->SetLimits(cs);
  cs.p_sk->LF()->Compute(cs);
  return cs.p_sk->LF()->Construct(cs,1);
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
  if (back()!=p) abort(); 
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

namespace DIRE {

  std::ostream &operator<<(std::ostream &s,const Amplitude &a)
  {
    Vec4D p;
    int c[4]={0,0,0};
    s<<"("<<&a<<"): t = "<<a.T()<<", t0 = "<<a.T0()
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
