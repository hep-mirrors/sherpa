#include "MCATNLO/Main/CS_MCatNLO.H"

#include "MCATNLO/Main/CS_Gamma.H"
#include "MCATNLO/Showers/Splitting_Function_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "MODEL/Main/Model_Base.H"
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

CS_MCatNLO::CS_MCatNLO(PDF::ISR_Handler *const _isr,
		       MODEL::Model_Base *const model,
		       Data_Reader *const _dataread) : 
  NLOMC_Base("MC@NLO_CSS"), p_isr(_isr), 
  p_powheg(NULL), p_as(NULL), p_cluster(NULL), p_gamma(NULL)
{
  m_psmode=_dataread->GetValue<int>("NLO_CSS_PSMODE",0);
  m_maxem=_dataread->GetValue<int>("NLO_CSS_MAXEM",1);
  m_scale2fac  = _dataread->GetValue<double>("CSS_SHOWER_SCALE2_FACTOR",-1.);
  if (m_scale2fac>0. && m_scale2fac!=1.) {
    p_as = (MODEL::Running_AlphaS*)model->GetScalarFunction("alpha_S");
  }
  SF_Lorentz::SetKappa(_dataread->GetValue<double>("DIPOLE_KAPPA",2.0/3.0));

  p_powheg = new Shower(_isr,0,_dataread);
  p_next = new All_Singlets();
  p_cluster = new CS_Cluster_Definitions(p_powheg,1);
  p_gamma = new CS_Gamma(this,p_powheg,p_cluster);
  p_gamma->SetOEF(_dataread->GetValue<double>("CSS_OEF",9.0));
  p_powheg->SetGamma(p_gamma);
  m_kt2min=p_powheg->GetSudakov()->ISPT2Min();
}

CS_MCatNLO::~CS_MCatNLO() 
{
  CleanUp();
  if (p_powheg) delete p_powheg;
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
  m_nem=0;
  m_weight=1.0;
  CleanUp();
  PrepareMCatNLO(ampl);
  int stat(PerformMCatNLO(m_maxem,m_nem));
  if (m_nem) {
    Cluster_Amplitude *rampl(GetRealEmissionAmplitude());
    rampl->SetNext(ampl);
    size_t idnew(rampl->IdNew());
    rampl->SetIdNew(0);
    Parton * const* last(p_powheg->GetLast());
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

int CS_MCatNLO::PerformMCatNLO(const size_t &maxem,size_t &nem)
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
	if (m_psmode &&
	    !(((*cit)->GetFlow(1) &&
	       (*cit)->GetFlow(1)==(*pit)->GetFlow(2)) ||
	      ((*cit)->GetFlow(2) &&
	       (*cit)->GetFlow(2)==(*pit)->GetFlow(1)))) continue;
	msg_Debugging()<<(*pit)->GetFlavour()<<ID((*pit)->Id())<<" ";
	(*cit)->Specs().push_back(*pit);
      }
    }
    msg_Debugging()<<"-> "<<(*cit)->Specs().size()<<" dipole(s)\n";
  }
  p_gamma->SetOn(1);
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
    msg_Debugging()<<"before powheg step\n";
    msg_Debugging()<<**sit;
    size_t pem(nem);
    if (!p_powheg->EvolveShower(*sit,maxem,nem)) return 0;
    m_weight*=p_powheg->Weight();
    msg_Debugging()<<"after powheg step with "<<nem-pem
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
  p_powheg->SetMS(p_ms);
  // p_as is alphaS in case we do a shower variation it is != NULL.
  if (p_as && m_scale2fac!=-1. && m_scale2fac!=1.) {
    double mu2=(ampl->Q2());
    double newasmu2((*p_as)(m_scale2fac*mu2));
    // fixes an operator (*p_as)(qt^2,true) to be with new Lambda2
    p_as->FixShowerLambda2(mu2,newasmu2,p_as->Nf(mu2),p_as->Order());
    p_powheg->SetCouplingMax();
    // msg_Out()<<METHOD<<" tests alphaS(order = "<<p_as->Order()<<", "
    // 	     <<"mu = "<<sqrt(mu2)<<"): \n"
    // 	     <<" mu alphaS(mu) alphaS(mu,new)\n"
    // 	     <<"  5. "<<(*p_as)(25.)<<"  "<<p_as->AlphaS(25.,true)<<"\n"
    // 	     <<" 10. "<<(*p_as)(100.)<<"  "<<p_as->AlphaS(100.,true)<<"\n"
    // 	     <<" 25. "<<(*p_as)(625.)<<"  "<<p_as->AlphaS(625.,true)<<"\n"
    // 	     <<" 50. "<<(*p_as)(2500.)<<"  "<<p_as->AlphaS(2500.,true)<<"\n"
    // 	     <<"100. "<<(*p_as)(10000.)<<"  "<<p_as->AlphaS(10000.,true)<<"\n";
  }

  return true;
}

Singlet *CS_MCatNLO::TranslateAmplitude
(Cluster_Amplitude *const ampl,
 std::map<Cluster_Leg*,Parton*> &pmap,std::map<Parton*,Cluster_Leg*> &lmap)
{
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
    parton->SetStart(ampl->Q2());
    double ktveto2(jf?jf->Ycut()*sqr(rpa->gen.Ecms()):parton->KtStart());
    double ktmax2(ampl->Legs().size()-ampl->NIn()+1==
		  ampl->Leg(2)->NMax()?parton->KtStart():0.0);
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
    if ((*it)->GetType()==pst::IS)
      ampl->CreateLeg
	(-(*it)->Momentum(),(*it)->GetFlavour().Bar(),
	 mode==0?ColorID((*it)->GetFlow(1),(*it)->GetFlow(2)):
	 ColorID((*it)->GetMEFlow(1),(*it)->GetMEFlow(2)),
	 (*it)->Id()?(*it)->Id():ampl->IdNew());
    ampl->Legs().back()->SetStat(1);
    ampl->Legs().back()->SetNMax(p_rampl->Leg(2)->NMax());
  }
  for (Singlet::const_iterator
	 it(sing->begin());it!=sing->end();++it) {
    if ((*it)->GetType()==pst::FS)
      ampl->CreateLeg
	((*it)->Momentum(),(*it)->GetFlavour(),
	 mode==0?ColorID((*it)->GetFlow(1),(*it)->GetFlow(2)):
	 ColorID((*it)->GetMEFlow(1),(*it)->GetMEFlow(2)),
	 (*it)->Id()?(*it)->Id():ampl->IdNew());
    ampl->Legs().back()->SetStat(1);
    ampl->Legs().back()->SetNMax(p_rampl->Leg(2)->NMax());
  }
  ampl->SetKT2(p_rampl->KT2());
  ampl->SetNewCol(p_powheg->GetLast()[3]->Color().m_new);
  Process_Base::SortFlavours(ampl);
  return ampl;
}

void CS_MCatNLO::AddRBPoint(ATOOLS::Cluster_Amplitude *const ampl)
{
  p_gamma->AddRBPoint(ampl);
}

namespace PDF {

  DECLARE_GETTER(CSS_MCatNLO_Getter,"MC@NLO_CSS",NLOMC_Base,NLOMC_Key);

  NLOMC_Base *CSS_MCatNLO_Getter::operator()(const NLOMC_Key &key) const
  {
    return new CS_MCatNLO(key.p_isr,key.p_model,key.p_read);
  }

  void CSS_MCatNLO_Getter::PrintInfo(std::ostream &str,const size_t width) const
  { 
    str<<"The CSS MC@NLO generator"; 
  }

}
