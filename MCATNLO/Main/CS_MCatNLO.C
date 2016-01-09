#include "MCATNLO/Main/CS_MCatNLO.H"

#include "MCATNLO/Main/CS_Gamma.H"
#include "MCATNLO/Showers/Splitting_Function_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace MCATNLO;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;
using namespace std;

CS_MCatNLO::CS_MCatNLO(PDF::ISR_Handler *const _isr,MODEL::Model_Base *const model,
		     Data_Reader *const _dataread) : 
  NLOMC_Base("MC@NLO_CSS"), p_isr(_isr), 
  p_mcatnlo(NULL), p_cluster(NULL), p_gamma(NULL)
{
  m_psmode=_dataread->GetValue<int>("NLO_CSS_PSMODE",0);
  if (m_psmode) msg_Info()<<METHOD<<"(): Set PS mode "<<m_psmode<<".\n";
  m_maxweight=_dataread->GetValue<double>("NLO_CSS_MAXWEIGHT",1.0e3);
  if (m_maxweight!=1.0e3) msg_Info()<<METHOD<<"(): Set max weight "<<m_maxweight<<".\n";
  m_maxem=_dataread->GetValue<int>("NLO_CSS_MAXEM",1);
  SF_Lorentz::SetKappa(_dataread->GetValue<double>("DIPOLE_KAPPA",2.0/3.0));

  p_mcatnlo = new Shower(_isr,0,_dataread);
  p_next = new All_Singlets();
  p_cluster = new CS_Cluster_Definitions(p_mcatnlo,1);
  p_gamma = new CS_Gamma(this,p_mcatnlo,p_cluster);
  p_gamma->SetOEF(_dataread->GetValue<double>("CSS_OEF",3.0));
  p_mcatnlo->SetGamma(p_gamma);
  m_kt2min=p_mcatnlo->GetSudakov()->ISPT2Min();
}

CS_MCatNLO::~CS_MCatNLO() 
{
  CleanUp();
  if (p_mcatnlo) delete p_mcatnlo;
  if (p_cluster) delete p_cluster;
  if (p_gamma) delete p_gamma;
  delete p_next;
}

void CS_MCatNLO::CleanUp()
{
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
    if (*sit) delete *sit;
  }
  m_allsinglets.clear();
}

int CS_MCatNLO::GeneratePoint(Cluster_Amplitude *const ampl) 
{
  DEBUG_FUNC("");
  for (double qfac(1.0);;qfac*=10.0) {
  m_nem=0;
  m_weight=1.0;
  CleanUp();
  PrepareMCatNLO(ampl);
  int stat(PerformMCatNLO(m_maxem,m_nem,qfac));
  if (dabs(m_weight)>m_maxweight) {
    msg_Info()<<METHOD<<"(): Weight is "<<m_weight
	      <<". Retry with charge factor "<<qfac*10.0<<".\n";
    continue;
  }
  if (m_nem) {
    Cluster_Amplitude *rampl(GetRealEmissionAmplitude());
    rampl->SetNext(ampl);
    size_t idnew(rampl->IdNew());
    rampl->SetIdNew(0);
    Parton * const* last(p_mcatnlo->GetLast());
    while (rampl->Next()) {
      rampl=rampl->Next();
      for (size_t i(0);i<rampl->Legs().size();++i) {
	rampl->Leg(i)->SetNMax(rampl->Leg(i)->NMax());
	size_t cid(rampl->Leg(i)->Id());
	if (cid&last[0]->Id()) {
	  for (size_t j(0);j<rampl->Legs().size();++j)
	    if (rampl->Leg(j)->K()==cid)
	      rampl->Leg(j)->SetK(cid|idnew);
	  rampl->Leg(i)->SetId(cid|idnew);
	  if (rampl->Prev()->Prev()==NULL) {
	    rampl->Leg(i)->SetK(last[2]->Id());
	    ampl->Prev()->SetIdNew(idnew);
	  }
	  break;
	}
      }
    }
  }
  return stat;
  }
  THROW(fatal_error,"Internal error");
  return false;
}

int CS_MCatNLO::PerformMCatNLO(const size_t &maxem,size_t &nem,const double &qfac)
{
  SF_Coupling::SetQFac(qfac);
  std::set<Parton*> nxs;
  Singlet *last(*(m_allsinglets.end()-1));
  std::string pname(Process_Base::GenerateName(p_rampl));
  const IDip_Set &iinfo((*p_rampl->IInfo<StringIDipSet_Map>())[pname]);
  for (Singlet::iterator cit(last->begin());cit!=last->end();++cit) {
    msg_Debugging()<<"filling partner list for "<<(*cit)->GetFlavour()
		   <<ID((*cit)->Id())<<" ... ";
    for (Singlet::iterator pit(last->begin());pit!=last->end();++pit) {
      if (iinfo.find(IDip_ID((*cit)->Idx(),(*pit)->Idx()))!=iinfo.end()) {
	if (m_psmode &&
	    !(((*cit)->GetFlow(1) &&
	       (*cit)->GetFlow(1)==(*pit)->GetFlow(2)) ||
	      ((*cit)->GetFlow(2) &&
	       (*cit)->GetFlow(2)==(*pit)->GetFlow(1)))) continue;
	msg_Debugging()<<(*pit)->GetFlavour()<<ID((*pit)->Id())<<" ";
	(*cit)->Specs().push_back(*pit);
      }
    }
    if ((*cit)->GetFlavour().StrongCharge()==8 &&
	(*cit)->Specs().size()<2) SF_Coupling::SetQFac(2.0*qfac);
    msg_Debugging()<<"-> "<<(*cit)->Specs().size()<<" dipole(s)\n";
  }
  p_gamma->SetOn(1);
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
    msg_Debugging()<<"before mc@nlo step\n";
    msg_Debugging()<<**sit;
    size_t pem(nem);
    if (!p_mcatnlo->EvolveShower(*sit,maxem,nem)) return 0;
    m_weight*=p_mcatnlo->Weight();
    msg_Debugging()<<"after mc@nlo step with "<<nem-pem
		   <<" emission(s), w = "<<m_weight<<"\n";
    msg_Debugging()<<**sit;
    msg_Debugging()<<"\n";
  }
  return 1;
}

bool CS_MCatNLO::PrepareMCatNLO(Cluster_Amplitude *const ampl)
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
  p_mcatnlo->SetMS(p_ms);
  return true;
}

Singlet *CS_MCatNLO::TranslateAmplitude
(Cluster_Amplitude *const ampl,
 std::map<Cluster_Leg*,Parton*> &pmap,std::map<Parton*,Cluster_Leg*> &lmap)
{
  double muQ2(ampl->MuQ2());
  for (Cluster_Amplitude *campl(ampl);
       campl->Next();campl=campl->Next())
    if (campl->Next()->OrderQCD()<campl->OrderQCD()) {
      muQ2=campl->KT2();
      break;
    }
  PHASIC::Jet_Finder *jf(ampl->JF<PHASIC::Jet_Finder>());
  Singlet *singlet(new Singlet());
  singlet->SetMS(p_ms);
  singlet->SetProcs(ampl->Procs<void>());
  singlet->SetMuR2(ampl->MuR2());
  CI_Map col(ampl->ColorMap());
  col[0]=0;
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
    CI_Map::const_iterator ci(col.find(parton->GetFlow(1)));
    CI_Map::const_iterator cj(col.find(parton->GetFlow(2)));
    if (ci!=col.end()) parton->SetMEFlow(1,ci->second);
    else parton->SetMEFlow(1,0);
    if (cj!=col.end()) parton->SetMEFlow(2,cj->second);
    else parton->SetMEFlow(2,0);
    parton->SetKin(p_mcatnlo->KinScheme());
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
    parton->SetStart(muQ2);
    double ktveto2(jf?jf->Ycut()*sqr(rpa->gen.Ecms()):parton->KtStart());
    double ktmax2(ampl->Legs().size()-ampl->NIn()+1==
		  ampl->Leg(0)->NMax()?parton->KtStart():0.0);
    parton->SetKtMax(ktmax2);
    parton->SetVeto(ktveto2);
    singlet->push_back(parton);
    parton->SetSing(singlet);
  }
  return singlet;
}

ATOOLS::Cluster_Amplitude *CS_MCatNLO::
GetRealEmissionAmplitude(const int mode)
{
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  Singlet *sing(*(m_allsinglets.end()-1));
  ampl->CopyFrom(p_rampl,1);
  ampl->SetProcs(sing->Procs<void>());
  ampl->SetIdNew(1<<(sing->size()-1));
  for (Singlet::const_iterator
	 it(sing->begin());it!=sing->end();++it) {
    if ((*it)->GetType()==pst::IS) {
      ampl->CreateLeg
	(-(*it)->Momentum(),(*it)->GetFlavour().Bar(),
	 mode==0?ColorID((*it)->GetFlow(1),(*it)->GetFlow(2)):
	 ColorID((*it)->GetMEFlow(1),(*it)->GetMEFlow(2)),
	 (*it)->Id()?(*it)->Id():ampl->IdNew());
      ampl->Legs().back()->SetNMax
	(p_rampl->IdLeg((*it)->Id()?(*it)->Id():1)->NMax());
    }
  }
  for (Singlet::const_iterator
	 it(sing->begin());it!=sing->end();++it) {
    if ((*it)->GetType()==pst::FS) {
      ampl->CreateLeg
	((*it)->Momentum(),(*it)->GetFlavour(),
	 mode==0?ColorID((*it)->GetFlow(1),(*it)->GetFlow(2)):
	 ColorID((*it)->GetMEFlow(1),(*it)->GetMEFlow(2)),
	 (*it)->Id()?(*it)->Id():ampl->IdNew());
      ampl->Legs().back()->SetNMax
	(p_rampl->IdLeg((*it)->Id()?(*it)->Id():1)->NMax());
    }
  }
  ampl->SetKT2(p_mcatnlo->GetLast()[3]->KtTest());
  ampl->SetNewCol(p_mcatnlo->GetLast()[3]->Color().m_new);
  Process_Base::SortFlavours(ampl);
  return ampl;
}

void CS_MCatNLO::AddRBPoint(ATOOLS::Cluster_Amplitude *const ampl)
{
  p_gamma->AddRBPoint(ampl);
}

DECLARE_GETTER(CS_MCatNLO,"MC@NLO_CSS",NLOMC_Base,NLOMC_Key);

NLOMC_Base *ATOOLS::Getter<NLOMC_Base,NLOMC_Key,CS_MCatNLO>::
operator()(const NLOMC_Key &key) const
{
  return new CS_MCatNLO(key.p_isr,key.p_model,key.p_read);
}

void ATOOLS::Getter<NLOMC_Base,NLOMC_Key,CS_MCatNLO>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The CSS MC@NLO generator"; 
}
