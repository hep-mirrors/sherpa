#include "DIM/Main/MCatNLO.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Exception.H"

using namespace DIM;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

MCatNLO::MCatNLO(const NLOMC_Key &key):
  NLOMC_Base("Dire"), p_ms(NULL)
{
  m_subtype=1;
  p_mcatnlo = new Shower();
  p_gamma = new Gamma(this,p_mcatnlo);
  p_mcatnlo->SetGamma(p_gamma);
  p_mcatnlo->Init(key.p_model,key.p_isr,key.p_reader);
  m_psmode=key.p_reader->Get<int>
    ("NLO_CSS_PSMODE",0,"PS mode",METHOD);
  m_wcheck=key.p_reader->Get<int>
    ("NLO_CSS_WEIGHT_CHECK",0,"Weight check",METHOD);
  for (int i(0);i<2;++i) m_kt2min[i]=p_mcatnlo->TMin(i);
}

MCatNLO::~MCatNLO()
{
  delete p_mcatnlo;
  delete p_gamma;
}

int MCatNLO::GeneratePoint(Cluster_Amplitude *const ampl)
{
  DEBUG_FUNC(this);
  p_mcatnlo->SetVariations(p_variationweights);
  m_weight=1.0;
  CleanUp();
  PrepareShower(ampl);
  if (p_rampl->NLO()&4) return 1;
  unsigned int nem=0;
  int stat(p_mcatnlo->Evolve(*m_ampls.back(),nem));
  m_weight*=p_mcatnlo->GetWeight();
  if (m_wcheck && dabs(m_weight)>m_maxweight) {
    m_maxweight=dabs(m_weight);
    std::string rname="direnlo.random.dat";
    if (ATOOLS::msg->LogFile()!="")
      rname=ATOOLS::msg->LogFile()+"."+rname;
    ATOOLS::ran->WriteOutSavedStatus(rname.c_str());
    std::ofstream outstream(rname.c_str(),std::fstream::app);
    outstream<<std::endl;
    outstream<<"# Wrote status for weight="<<m_weight<<" in event "
	     <<rpa->gen.NumberOfGeneratedEvents()+1<<std::endl;
    outstream.close();
  }
  if (stat!=1) return stat;
  if (nem) {
    Cluster_Amplitude *rampl(GetRealEmissionAmplitude());
    rampl->SetNext(ampl);
    size_t idnew(rampl->IdNew());
    rampl->SetIdNew(0);
    const Splitting &s(p_mcatnlo->LastSplitting());
    while (rampl->Next()) {
      rampl=rampl->Next();
      for (size_t i(0);i<rampl->Legs().size();++i) {
	size_t cid(rampl->Leg(i)->Id());
	if (cid&(1<<(s.p_c->Id()-1))) {
	  for (size_t j(0);j<rampl->Legs().size();++j)
	    if (rampl->Leg(j)->K()==cid)
	      rampl->Leg(j)->SetK(cid|idnew);
	  rampl->Leg(i)->SetId(cid|idnew);
	  if (rampl->Prev()->Prev()==NULL) {
	    rampl->Leg(i)->SetK(1<<(s.p_s->Id()-1));
	    ampl->Prev()->SetIdNew(idnew);
	  }
	  break;
	}
      }
    }
  }
  return stat;
}

ATOOLS::Cluster_Amplitude *MCatNLO::
GetRealEmissionAmplitude(const int mode)
{
  int nmax(p_rampl->Legs().back()->NMax());
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  ampl->CopyFrom(p_rampl,1);
  ampl->SetProcs(p_rampl->Procs<void>());
  ampl->SetIdNew(1<<(m_ampls.back()->size()-1));
  for (Amplitude::const_iterator pit(m_ampls.back()->begin());
       pit!=m_ampls.back()->end();++pit) {
    ampl->CreateLeg((*pit)->Mom(),(*pit)->Flav(),
		    ColorID((*pit)->Col().m_i,(*pit)->Col().m_j),
		    1<<((*pit)->Id()-1));
    ampl->Legs().back()->SetNMax(nmax);
  }
  ampl->SetKT2(p_mcatnlo->LastSplitting().m_t);
  ampl->SetMuQ2(p_rampl->KT2());
  Process_Base::SortFlavours(ampl);
  return ampl;
}

void MCatNLO::CleanUp()
{
  for (Amplitude_Vector::const_iterator it(m_ampls.begin());
       it!=m_ampls.end();++it) delete *it;
  m_ampls.clear();
}

Amplitude *MCatNLO::Convert
(Cluster_Amplitude *const campl,
 std::map<Cluster_Leg*,Parton*> &lmap)
{
  Amplitude *ampl(new Amplitude(campl));
  ampl->SetT(campl->KT2());
  if (campl->Prev()) ampl->SetT0(campl->Prev()->KT2());
  for (size_t i(0);i<campl->Legs().size();++i) {
    Cluster_Leg *cl(campl->Leg(i));
    Parton *p(new Parton(ampl,cl->Flav(),cl->Mom(),
			 Color(cl->Col().m_i,cl->Col().m_j)));
    ampl->push_back(p);
    p->SetId(p->Counter());
    if (i<campl->NIn()) p->SetBeam(1+(cl->Mom()[3]>0.0));
    lmap[cl]=p;
  }
  msg_Debugging()<<*ampl<<"\n";
  return ampl;
}

bool MCatNLO::PrepareShower
(Cluster_Amplitude *const ampl,const bool &soft)
{
  DEBUG_FUNC(soft);
  p_rampl=ampl;
  p_ms=ampl->MS();
  p_mcatnlo->SetMS(p_ms);
  std::map<Cluster_Leg*,Parton*> lmap;
  m_ampls.push_back(Convert(ampl,lmap));
  std::string pname(Process_Base::GenerateName(p_rampl));
  const IDip_Set &iinfo((*p_rampl->IInfo<StringIDipSet_Map>())[pname]);
  for (size_t i(0);i<m_ampls.back()->size();++i) {
    Parton *c((*m_ampls.back())[i]);
    msg_Debugging()<<"spectators for "
		   <<c->Flav()<<"("<<c->Id()<<") ... ";
    for (size_t j(0);j<m_ampls.back()->size();++j) {
      Parton *s((*m_ampls.back())[j]);
      if (iinfo.find(IDip_ID(c->Id()-1,s->Id()-1))!=iinfo.end()) {
	msg_Debugging()<<s->Flav()<<"("<<s->Id()<<") ";
	c->S().push_back(s);
      }
    }
    msg_Debugging()<<"-> "<<c->S().size()<<" dipole(s)\n";
  }
  if (ampl->NIn()+ampl->Leg(2)->NMax()==
      ampl->Legs().size()+1) m_ampls.back()->SetJF(NULL);
  Cluster_Amplitude *campl(ampl);
  while (campl->Next()) campl=campl->Next();
  return true;
}

double MCatNLO::KT2(const ATOOLS::NLO_subevt &sub,
		    const double &x,const double &y,const double &Q2)
{
  double mi2(sqr(sub.p_real->p_fl[sub.m_i].Mass()));
  double mj2(sqr(sub.p_real->p_fl[sub.m_j].Mass()));
  double mk2(sqr(sub.p_real->p_fl[sub.m_k].Mass()));
  double mij2(sqr(sub.p_fl[sub.m_ijt].Mass()));
  if (sub.m_ijt>=2) {
    double t;
    if (sub.m_kt>=2) t=(Q2-mi2-mj2-mk2)*y*(1.0-y);
    else t=(-Q2+mi2+mj2+mk2)*(1.0-y*(Q2-mi2-mj2-mk2)/(Q2-mij2-mk2));
    if (sub.p_real->p_fl[sub.m_i].IsGluon()) {
      if (!sub.p_real->p_fl[sub.m_j].IsGluon()) return t*x;
      // approximate, need to split g->gg kernel
      return t*Min(1.0-x,x);
    }
    else {
      if (sub.p_real->p_fl[sub.m_j].IsGluon()) return t*(1.0-x);
      // approximate, need to split g->qq kernel
      return t*Min(1.0-x,x);
    }
  }
  if (sub.m_ijt<2 && sub.m_kt>=2) {
    return (-Q2+mi2+mj2+mk2)*y*(1.0-x);
  }
  if (sub.m_ijt<2 && sub.m_kt<2) {
    return (Q2-mi2-mj2-mk2)*y*(1.0-x-y);
  }
  THROW(fatal_error,"Implement me");
}

DECLARE_GETTER(MCatNLO,"MC@NLO_Dire",NLOMC_Base,NLOMC_Key);

NLOMC_Base *Getter<NLOMC_Base,NLOMC_Key,MCatNLO>::
operator()(const NLOMC_Key &key) const
{
  return new MCatNLO(key);
}

void Getter<NLOMC_Base,NLOMC_Key,MCatNLO>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Dire MC@NLO generator"; 
}
