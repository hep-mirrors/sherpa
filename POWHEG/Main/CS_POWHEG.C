#include "POWHEG/Main/CS_POWHEG.H"

#include "POWHEG/Main/CS_Gamma.H"
#include "POWHEG/Showers/Splitting_Function_Base.H"
#include "PHASIC++/Process/POWHEG_Process.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace POWHEG;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;
using namespace std;

CS_POWHEG::CS_POWHEG(PDF::ISR_Handler *const _isr,MODEL::Model_Base *const model,
		     Data_Reader *const _dataread) : 
  NLOMC_Base("POWHEG_CSS"), p_isr(_isr), 
  p_powheg(NULL), p_cluster(NULL), p_gamma(NULL)
{
  m_maxem=_dataread->GetValue<int>("PH_CSS_MAXEM",1);
  double zhth(_dataread->GetValue<double>("PH_CSS_ZH_THRESHOLD",100.0));
  double ktres(_dataread->GetValue<double>("PH_CSS_KT_RESOLUTION",4.0));
  p_powheg = new Shower(_isr,0,_dataread);
  p_next = new All_Singlets();
  p_cluster = new CS_Cluster_Definitions(p_powheg,1);
  p_gamma = new CS_Gamma(this,p_powheg,p_cluster);
  p_powheg->SetGamma(p_gamma);
  p_gamma->SetZHParams(zhth,ktres);
  m_kt2min=p_powheg->GetSudakov()->PT2Min();
}

CS_POWHEG::~CS_POWHEG() 
{
  CleanUp();
  if (p_powheg) delete p_powheg;
  if (p_cluster) delete p_cluster;
  if (p_gamma) delete p_gamma;
  delete p_next;
}

void CS_POWHEG::CleanUp()
{
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
    if (*sit) delete *sit;
  }
  m_allsinglets.clear();
}

int CS_POWHEG::GeneratePoint(Cluster_Amplitude *const ampl) 
{
  CleanUp();
  PreparePOWHEG(ampl);
  m_nem=0;
  int stat(PerformPOWHEG(m_maxem,m_nem));
  if (m_nem) {
    Cluster_Amplitude *rampl(GetRealEmissionAmplitude());
    rampl->SetNext(ampl);
    rampl->SetRBMap(NULL);
    size_t idnew(rampl->IdNew());
    rampl->SetIdNew(0);
    Parton * const* last(p_powheg->GetLast());
    while (rampl->Next()) {
      rampl=rampl->Next();
      rampl->SetRBMap(NULL);
      for (size_t i(0);i<rampl->Legs().size();++i) {
	rampl->Leg(i)->SetNMax(rampl->Leg(i)->NMax()+1);
	size_t cid(rampl->Leg(i)->Id());
	if (cid&last[0]->Id()) {
	  for (size_t j(0);j<rampl->Legs().size();++j)
	    if (rampl->Leg(j)->K()==cid)
	      rampl->Leg(j)->SetK(cid|idnew);
	  rampl->Leg(i)->SetId(cid|idnew);
	  if (rampl->Prev()->Prev()==NULL) {
	    rampl->Leg(i)->SetK(last[2]->Id());
	    rampl->SetIdNew(rampl->Leg(i)->Id());
	  }
	  break;
	}
      }
    }
  }
  return stat;
}

int CS_POWHEG::PerformPOWHEG(const size_t &maxem,size_t &nem)
{
  std::set<Parton*> nxs;
  Singlet *last(*(m_allsinglets.end()-1));
  std::string pname(Process_Base::GenerateName(p_rampl));
  const IDip_Set &iinfo((*p_rampl->IInfo<StringIDipSet_Map>())[pname]);
  for (Singlet::iterator cit(last->begin());cit!=last->end();++cit) {
    msg_Debugging()<<"filling partner list for "<<(*cit)->GetFlavour()
		   <<ID((*cit)->Id())<<" ... ";
    for (Singlet::iterator pit(last->begin());pit!=last->end();++pit) {
      if (iinfo.find(IDip_ID((*cit)->Idx(),(*pit)->Idx()))!=iinfo.end()) {
	int scc(((*cit)->Idx()<p_rampl->NIn()?(*cit)->GetFlavour().Bar():
		 (*cit)->GetFlavour()).StrongCharge());
	int scp(((*pit)->Idx()<p_rampl->NIn()?(*pit)->GetFlavour().Bar():
		 (*pit)->GetFlavour()).StrongCharge());
	if (abs(scc)==3 && abs(scp)==3 && scc!=-scp) continue;
	msg_Debugging()<<(*pit)->GetFlavour()<<ID((*pit)->Id())<<" ";
	(*cit)->Specs().push_back(*pit);
      }
    }
    msg_Debugging()<<"-> "<<(*cit)->Specs().size()<<" dipole(s)\n";
  }
  p_gamma->SetOn(1);
  p_powheg->SetRB((*m_allsinglets.begin())->RBMap());
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
    msg_Debugging()<<"before powheg step\n";
    msg_Debugging()<<**sit;
    size_t pem(nem);
    if (!p_powheg->EvolveShower(*sit,maxem,nem)) return 0;
    msg_Debugging()<<"after powheg step with "<<nem-pem<<" emission(s)\n";
    msg_Debugging()<<**sit;
    msg_Debugging()<<"\n";
  }
  return 1;
}

bool CS_POWHEG::PreparePOWHEG(Cluster_Amplitude *const ampl)
{
  CleanUp();
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  p_rampl=ampl;
  p_ms=ampl->MS();
  p_next->clear();
  m_allsinglets.clear();
  Cluster_Amplitude *campl(ampl);
  msg_Debugging()<<*campl<<"\n";
  std::map<Parton*,Cluster_Leg*> lmap;
  std::map<Cluster_Leg*,Parton*> pmap;
  Singlet *sing(TranslateAmplitude(campl,pmap,lmap));
  m_allsinglets.push_back(sing);
  p_next->push_back(sing);
  msg_Debugging()<<"\nSinglet lists:\n\n";
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
    (*sit)->SetJF(ampl->JF<PHASIC::Jet_Finder>());
    (*sit)->SetShower(p_shower);
    (*sit)->SetAll(p_next);
    msg_Debugging()<<**sit;
    msg_Debugging()<<"\n";
  }
  msg_Debugging()<<"}\n";
  p_powheg->SetMS(p_ms);
  return true;
}

Singlet *CS_POWHEG::TranslateAmplitude
(Cluster_Amplitude *const ampl,
 std::map<Cluster_Leg*,Parton*> &pmap,std::map<Parton*,Cluster_Leg*> &lmap)
{
  PHASIC::Jet_Finder *jf(ampl->JF<PHASIC::Jet_Finder>());
  Singlet *singlet(new Singlet());
  singlet->SetMS(p_ms);
  singlet->SetRBMap(ampl->RBMap());
  singlet->SetProcs(ampl->Procs<void>());
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *cl(ampl->Leg(i));
    if (cl->Flav().IsHadron() && cl->Id()&((1<<ampl->NIn())-1)) continue;
    bool is(cl->Id()&((1<<ampl->NIn())-1));
    Particle p(1,is?cl->Flav().Bar():cl->Flav(),is?-cl->Mom():cl->Mom());
    if (cl->Col().m_i>0 || cl->Col().m_j>0) {
      if (is) {
	p.SetFlow(2,cl->Col().m_i);
	p.SetFlow(1,cl->Col().m_j);
      }
      else {
	p.SetFlow(1,cl->Col().m_i);
	p.SetFlow(2,cl->Col().m_j);
      }
    }
    Parton *parton(new Parton(&p,is?pst::IS:pst::FS));
    pmap[cl]=parton;
    lmap[parton]=cl;
    parton->SetIdx(i);
    parton->SetId(cl->Id());
    parton->SetRFlow();
    parton->SetKin(p_powheg->KinScheme());
    if (is) {
      if (Vec3D(p.Momentum())*Vec3D(rpa->gen.PBeam(0))>0.) {
	parton->SetXbj(p.Momentum()[0]/rpa->gen.PBeam(0)[0]);
	parton->SetBeam(0);
      }
      else { 
	parton->SetXbj(p.Momentum()[0]/rpa->gen.PBeam(1)[0]);
	parton->SetBeam(1);
      }
    }
    parton->SetStart(sqr(rpa->gen.Ecms()));
    double ktveto2(jf?jf->Ycut()*sqr(rpa->gen.Ecms()):parton->KtStart());
    double ktmax2(ampl->Legs().size()-ampl->NIn()==
		  ampl->Leg(2)->NMax()?parton->KtStart():0.0);
    parton->SetKtMax(ktmax2);
    parton->SetVeto(ktveto2);
    singlet->push_back(parton);
    parton->SetSing(singlet);
  }
  return singlet;
}

ATOOLS::Cluster_Amplitude *CS_POWHEG::GetRealEmissionAmplitude()
{
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  Singlet *sing(*(m_allsinglets.end()-1));
  ampl->CopyFrom(p_rampl,1);
  ampl->SetRBMap(sing->RBMap());
  ampl->SetProcs(sing->Procs<void>());
  ampl->SetIdNew(1<<(sing->size()-1));
  for (Singlet::const_iterator
	 it(sing->begin());it!=sing->end();++it) {
    if ((*it)->GetType()==pst::IS)
      ampl->CreateLeg
	(-(*it)->Momentum(),(*it)->GetFlavour().Bar(),
	 ColorID((*it)->GetFlow(1),(*it)->GetFlow(2)),
	 (*it)->Id()?(*it)->Id():ampl->IdNew());
    ampl->Legs().back()->SetStat(1);
    ampl->Legs().back()->SetNMax(p_rampl->Leg(2)->NMax()+1);
  }
  for (Singlet::const_iterator
	 it(sing->begin());it!=sing->end();++it) {
    if ((*it)->GetType()==pst::FS)
      ampl->CreateLeg
	((*it)->Momentum(),(*it)->GetFlavour(),
	 ColorID((*it)->GetFlow(1),(*it)->GetFlow(2)),
	 (*it)->Id()?(*it)->Id():ampl->IdNew());
    ampl->Legs().back()->SetStat(1);
    ampl->Legs().back()->SetNMax(p_rampl->Leg(2)->NMax()+1);
  }
  ampl->SetKT2(p_powheg->GetLast()[0]->KtStart());
  Process_Base::SortFlavours(ampl);
  return ampl;
}

void CS_POWHEG::AddRBPoint(ATOOLS::Cluster_Amplitude *const ampl)
{
  p_gamma->AddRBPoint(ampl);
}

SH_Pair CS_POWHEG::SHSplit(const double &B,const double &Qij2,
			   const RB_Data *rbd) const
{
  if (rbd->m_ktres<1.0) return SH_Pair(1.0,0.0);
  double S(Min(1.0,B/rbd->m_bmax));
  double H(Min(1.0,Qij2/sqr(rpa->gen.Ecms()/rbd->m_ktres)));
  return SH_Pair(S,H);
}

namespace PDF {

  DECLARE_GETTER(CSS_POWHEG_Getter,"POWHEG_CSS",NLOMC_Base,NLOMC_Key);

  NLOMC_Base *CSS_POWHEG_Getter::operator()(const NLOMC_Key &key) const
  {
    return new CS_POWHEG(key.p_isr,key.p_model,key.p_read);
  }

  void CSS_POWHEG_Getter::PrintInfo(std::ostream &str,const size_t width) const
  { 
    str<<"The CSS POWHEG generator"; 
  }

}
