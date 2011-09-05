#ifndef APACIC_Main_Apacic_H
#define APACIC_Main_Apacic_H

#include "PDF/Main/Shower_Base.H"
#include "AddOns/Apacic++/Showers/Initial_State_Shower.H"
#include "AddOns/Apacic++/Showers/Final_State_Shower.H"
#include "AddOns/Apacic++/Main/Tree.H"

#include "PDF/Main/ISR_Handler.H"

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "AddOns/Apacic++/Showers/Jet_Veto.H"
#include "ATOOLS/Phys/Cluster_Definitions_Base.H"

namespace ATOOLS { 
  class Cluster_Leg;
  class Data_Reader; 
}

namespace APACIC {

  class Apacic : public PDF::Shower_Base {
  private:

    Initial_State_Shower *p_inishower;
    Final_State_Shower   *p_finshower;

    Jet_Veto *p_jetveto;

    Tree **p_initrees;
    Tree  *p_fintree;

    ATOOLS::Poincare m_cfcms;

    bool m_isron, m_fsron, m_showers, m_last_ljv, m_did_isshower;
    int  m_wmode;

    double m_scale, m_fcs;

    ATOOLS::Cluster_Definitions_Base *p_cluster;
    ATOOLS::Cluster_Amplitude *p_ampl, *p_rampl;

    std::map<size_t,double> m_tmap;

    double HardScale(const ATOOLS::Cluster_Amplitude *const ampl);

    double CouplingWeight(const ATOOLS::Cluster_Leg *const ij,
			  const size_t &oqcd,
			  const double &kt2,const double &kt2r) const;

    double PseudoShowerWeight(ATOOLS::Cluster_Amplitude *const ampl);

    bool FillTrees(ATOOLS::Cluster_Amplitude *ampl);
    void SetFSProps(Knot *const mo,Knot *const left,Knot *const right);

    int TrialEmission(const int mode);

  public :
    
    // constructor
    Apacic(PDF::ISR_Handler *const isr,MODEL::Model_Base *const model,
	   ATOOLS::Data_Reader *const dataread);

    // destructor
    ~Apacic(); 

    // member functions
    int  PerformShowers();
    int  PerformDecayShowers();
    int  TrialEmission();

    double CouplingWeight(ATOOLS::Cluster_Amplitude *const ampl);
    double TrialWeight(ATOOLS::Cluster_Amplitude *const ampl);

    int  PerformShowers(const size_t &nem);
    bool ExtractPartons(ATOOLS::Blob_List *const blist);

    void CleanUp();
    void OutputTrees();

    // inline functions
    inline Tree  *FinTree() const  { return p_fintree;  }
    inline Tree **IniTrees() const { return p_initrees; }

    inline Final_State_Shower   *FinShower() const { return p_finshower; }
    inline Initial_State_Shower *IniShower() const { return p_inishower; }

    ATOOLS::Cluster_Definitions_Base * GetClusterDefinitions();
    ATOOLS::Cluster_Amplitude *GetRealEmissionAmplitude();
    bool PrepareShower(ATOOLS::Cluster_Amplitude *const ampl);
    double CalculateWeight(ATOOLS::Cluster_Amplitude *const ampl);
    
    std::string GetKT2(const std::string &jm2) const;

  };// end of class Apacic

}// end of namespace APACIC

#endif

#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "AddOns/Apacic++/Main/Apacic_Cluster_Definitions.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

#ifdef PROFILE__all
#define PROFILE__Apacic
#endif
#ifdef PROFILE__Apacic
#include "prof.hh"
#else 
#define PROFILE_HERE
#define PROFILE_LOCAL(LOCALNAME)
#endif

using namespace APACIC;
using namespace PDF;
using namespace ATOOLS;

Apacic::Apacic(ISR_Handler *const isr,MODEL::Model_Base *const model,
	       Data_Reader *const dataread):
  Shower_Base("Apacic"), p_inishower(NULL), p_finshower(NULL), 
  p_jetveto(NULL), p_initrees(NULL), p_fintree(NULL), 
  m_showers(false), p_cluster(NULL), p_ampl(NULL)
{
  rpa->gen.AddCitation
    (1,"Apacic is published under \\cite{Kuhn:2000dk,Krauss:2005re}.");
  Splitting_Function::SetKFactorScheme
    (ToType<int>(rpa->gen.Variable("S_KFACTOR_SCHEME","0"))&1);        
  m_fcs=dataread->GetValue<double>("APACIC_FCMODE",0.0);
  m_wmode=dataread->GetValue<int>("APACIC_SWMODE",0);
  m_fsron=bool(dataread->GetValue<int>("FSR_SHOWER",1));
  m_isron = 0;
  if (isr && isr->On()>0) m_isron=1;
  m_isron=bool(dataread->GetValue<int>("ISR_SHOWER",m_isron));
  if ((rpa->gen.Beam1().IsHadron() || rpa->gen.Beam2().IsHadron())
      && (m_fsron^m_isron)) 
    THROW(fatal_error,"Shower must be enabled for hadronic initial state.");
  if (m_fsron) {
    p_jetveto = new Jet_Veto();
    p_fintree   = new Tree();
    p_finshower = new Final_State_Shower(model,dataread);
    p_finshower->SetJetVeto(p_jetveto);
    p_finshower->SetWMode(m_wmode);
    m_showers=true;
  }
  if (m_isron) {
    if (isr->On()==0) THROW(fatal_error,"ISR must be enabled.");
    p_initrees  = new Tree*[2];
    for (int i=0;i<2;i++) p_initrees[i] = new Tree();
    p_inishower = new Initial_State_Shower(isr,p_finshower,model,dataread);
    m_showers=true;
  }
  m_scale=dataread->GetValue<double>("SHOWER_SCALE",-1.0);
  if (m_scale>0.0) {
    msg_Error()<<om::bold<<METHOD<<"(): SEVERE WARNING {\n"<<om::reset<<om::red
	       <<"  SHOWER_SCALE sets a fixed starting scale for Apacic++.\n"
	       <<"  It is used for internal testing purposes only.\n"
	       <<"  Users must never run Sherpa in this mode.\n"<<om::reset
	       <<om::bold<<"}"<<om::reset<<std::endl;
  }
}
  
Apacic::~Apacic() 
{
  if (p_fintree) delete p_fintree;
  if (p_initrees) {
    for (short unsigned int i(0);i<2;++i) delete p_initrees[i];
    delete [] p_initrees;
  }
  if (p_inishower) delete p_inishower;
  if (p_finshower) delete p_finshower;
  if (p_jetveto)   delete p_jetveto;
  if (p_cluster)   delete p_cluster;
}

int Apacic::PerformShowers() 
{
  m_did_isshower=true;
  return PerformShowers(std::numeric_limits<size_t>::max());
}

int Apacic::PerformDecayShowers() 
{
  m_did_isshower=false;
  return PerformShowers(std::numeric_limits<size_t>::max());
}

int Apacic::PerformShowers(const size_t &nem) 
{
  if (!m_showers) return 1;
  p_finshower->Sudakov()->SetNMax(nem);
  p_finshower->Sudakov()->ResetN();
  if (m_isron) {
    p_initrees[0]->Store();
    p_initrees[1]->Store();
  }
  if (m_fsron) p_fintree->Store();
    if (msg_LevelIsDebugging()) {
      msg_Out()<<"Apacic::PerformShowers : Before showering."<<std::endl;
      OutputTrees();
    }
    if (m_isron && m_did_isshower) {
      p_inishower->InitShowerPT(p_initrees[0]->GetRoot()->maxpt2);
      if (p_inishower->PerformShower(p_initrees,p_fintree)==0) {
	if (m_fsron) p_fintree->ClearStore();
	p_initrees[0]->ClearStore();
	p_initrees[1]->ClearStore();
	return 0;
      }
    }
    if (m_fsron) {
      int res(p_finshower->PerformShower(p_fintree));
      if (res!=1) {
	if (m_isron) {
	  p_initrees[0]->ClearStore();
	  p_initrees[1]->ClearStore();
	}
	p_fintree->ClearStore();
	return res;
      }
      p_finshower->SetAllColours(p_fintree->GetRoot());
    }
  if (m_isron) {
    p_initrees[0]->ClearStore();
    p_initrees[1]->ClearStore();
  }
  if (m_fsron) p_fintree->ClearStore();
  p_fintree->CheckMomentumConservation();
  if (m_isron) {
    p_initrees[0]->CheckMomentumConservation();
    p_initrees[1]->CheckMomentumConservation();
  }
  msg_Debugging()<<"kinematics check passed"<<std::endl;
  int number(0);
  Vec4D sum_fs(p_finshower->GetMomentum(p_fintree->GetRoot(),number));
  if (number<0) {
    msg_Error()<<METHOD<<"(..): Four Momentum not conserved. Abort."
	       <<std::endl;
    return 0;
  }
  m_last_ljv=false;
  if (msg_LevelIsDebugging()) {
    msg_Out()<<"Apacic::PerformShowers : After showering."<<std::endl;
    OutputTrees();
  }
  return 1;
}

void Apacic::CleanUp() 
{
  if (m_fsron) p_fintree->Reset();
  if (m_isron) for (int i=0;i<2;i++) p_initrees[i]->Reset();
}

bool Apacic::ExtractPartons(Blob_List *const bl) 
{
  Blob *sb(bl->FindLast(btp::Shower));
  if (sb==NULL) THROW(fatal_error,"No Shower blob");
  sb->SetTypeSpec("APACIC++2.0");
  for (int i=0;i<sb->NInP();++i) 
    sb->InParticle(i)->SetStatus(part_status::decayed);
  for (int i=0;i<sb->NOutP();++i) 
    sb->OutParticle(i)->SetStatus(part_status::decayed);
  sb->SetStatus(blob_status::needs_beams |
		blob_status::needs_harddecays |
		blob_status::needs_hadronization);
  if (m_fsron) {
    if (p_fintree->CheckStructure(true)) 
      p_finshower->ExtractFinalState(sb,p_fintree->GetRoot());
    else return false;
  }
  if (m_isron && m_did_isshower) {
    for (int i=0;i<2;i++) {
      if (p_initrees[i]->CheckStructure(true)) 
	p_inishower->SingleExtract(sb,p_initrees[i]->GetInitiator());
      else return false;
    }
  }
  return true;
}

void Apacic::OutputTrees() 
{
  if (m_fsron) {
    p_fintree->CheckMomentumConservation();
    p_finshower->OutputTree(p_fintree);
  }
  if (m_isron) {
    p_initrees[0]->CheckMomentumConservation();
    p_initrees[1]->CheckMomentumConservation();
    p_inishower->OutputTree(p_initrees[0]);
    p_inishower->OutputTree(p_initrees[1]);
  }
}

ATOOLS::Cluster_Definitions_Base * Apacic::GetClusterDefinitions() 
{
  if (p_cluster==NULL) p_cluster = new Apacic_Cluster_Definitions();
  return p_cluster;
}

double Apacic::HardScale(const Cluster_Amplitude *const ampl)
{
  double mu2(0.0);
  for (size_t i(0);i<ampl->Legs().size();++i) {
    if (ampl->Leg(i)->Col().m_i==0 && ampl->Leg(i)->Col().m_j==0) continue;
    for (size_t j(i+1);j<ampl->Legs().size();++j) {
      if (ampl->Leg(j)->Col().m_i==0 && ampl->Leg(j)->Col().m_j==0) continue;
      mu2=Max(mu2,(ampl->Leg(i)->Mom()+ampl->Leg(j)->Mom()).Abs2());
      msg_Debugging()<<"hs: i = "<<i<<", j = "<<j<<", \\mu = "<<sqrt(mu2)<<"\n";
    }
  }
  return mu2;
}

bool Apacic::PrepareShower(Cluster_Amplitude *const ampl)
{
  p_jetveto->SetJetFinder(ampl->JF<PHASIC::Jet_Finder>());
  p_rampl=ampl;
  if (m_fsron) p_finshower->SetMS(ampl->MS());
  if (m_isron) p_inishower->SetMS(ampl->MS());
  // reset momenta to simple sum for tree filling
  Cluster_Amplitude *ref(ampl);
  std::map<size_t,Vec4D> pmap;
  for (size_t i(0);i<ref->Legs().size();++i) 
    pmap[ref->Leg(i)->Id()]=ref->Leg(i)->Mom();
  double x1(ref->X1()), x2(ref->X2());
  double shat((ref->Leg(0)->Mom()+ref->Leg(1)->Mom()).Abs2());
  msg_Debugging()<<"x1 = "<<x1<<", x2 = "<<x2<<", Q' = "<<sqrt(shat)<<"\n";
  do {
    double cshat((ref->Leg(0)->Mom()+ref->Leg(1)->Mom()).Abs2());
    for (size_t i(0);i<ref->Legs().size();++i) {
      size_t id(ref->Leg(i)->Id());
      if (pmap.find(id)==pmap.end()) {
	Cluster_Amplitude *prev(ref->Prev());
	for (size_t j(0);j<prev->Legs().size();++j)
	  if (prev->Leg(j)->Id()&id) pmap[id]+=prev->Leg(j)->Mom();
	if (id&3) {
	  if (i==0) x1*=cshat/shat;
	  else if (i==1) x2*=cshat/shat;
	  else THROW(fatal_error,"Invalid clustering");
	  shat=cshat;
	  msg_Debugging()<<ID(id)<<": x1 = "<<x1<<", x2 = "
			 <<x2<<", Q' = "<<sqrt(shat)<<"\n";
	}
      }
      ref->Leg(i)->SetMom(pmap[id]);
    }
    ref->SetX1(x1);
    ref->SetX2(x2);
    if (msg_LevelIsDebugging()) ref->Print();
    ref=ref->Next();
  } while (ref);
  return FillTrees(ampl);
}

bool Apacic::FillTrees(Cluster_Amplitude *ampl)
{
  int dir(1);
  if (dabs(ampl->Leg(0)->Mom()[3])>dabs(ampl->Leg(1)->Mom()[3]))
    dir=ampl->Leg(0)->Mom()[3]>0?1:-1;
  else dir=ampl->Leg(1)->Mom()[3]>0?-1:1;
  msg_Debugging()<<"dir = "<<dir<<"\n";
  std::map<size_t,Knot*> kmap;
  double q2(ampl->MuF2()), s(0.0);
  while (ampl->Next()) {
    ampl=ampl->Next();
    for (size_t i(0);i<ampl->Legs().size();++i) {
      const Cluster_Leg *cl(ampl->Leg(i));  
      q2=Max(q2,dabs(cl->Mom().Abs2()));
    }
  }
  q2=Max(q2,HardScale(ampl));
  msg_Debugging()<<"largest ps scale "<<sqrt(q2)<<"\n";
  if (m_scale>0.0) q2=m_scale;
  Knot *mo(p_fintree->NewKnot());
  // temporarily, will not work for decay chains
  mo->qjv=sqrt(4.0*ampl->MuR2());
  mo->shower=0;
  double tmin(q2);
  if (p_rampl->Leg(2)->NMax()+p_rampl->NIn()
      -p_rampl->Legs().size()>0) {
    tmin=sqr(mo->qjv);
  }
  else {
    if (ampl->Prev()) mo->qjv=sqrt(ampl->Prev()->KT2QCD());
  }
  while (ampl) {
    size_t isid((1<<ampl->NIn())-1);
    Knot *left(NULL), *right(NULL);
    for (size_t i(0);i<ampl->Legs().size();++i) {
      const Cluster_Leg *cl(ampl->Leg(i));
      if (kmap.find(cl->Id())==kmap.end()) {
	Flavour fl(cl->Id()&isid?cl->Flav().Bar():cl->Flav());
	Vec4D mom(cl->Id()&isid?-cl->Mom():cl->Mom());
	Particle p(1,fl,mom);
	if (i<ampl->NIn()) {
	  p.SetFlow(1,cl->Col().m_j);
	  p.SetFlow(2,cl->Col().m_i);
	}
	else {
	  p.SetFlow(1,cl->Col().m_i);
	  p.SetFlow(2,cl->Col().m_j);
	}
	Knot *k(NULL);
	if (p_initrees) {
	  if (cl->Id()&3) k=p_initrees[i]->NewKnot(&p);
	  else k=p_fintree->NewKnot(&p);
	}
	else k=p_fintree->NewKnot(&p);
	if (m_wmode==0) k->asme=(*MODEL::as)(ampl->MuR2());
	k->kn_id=cl->Id();
	k->part->SetInfo('H');
	k->part->SetStatus(part_status::active);
	if (IdCount(cl->Id())==1) k->tout=sqr(ampl->MS()->Mass(fl));
	else {
	  k->tout=cl->Mom().Abs2();
	  tmin=Min(tmin,k->tout);
	}
	k->E2=sqr(mom[0]);
	if (ampl->Legs().size()>ampl->NIn()+2) k->t=-1.0;
	else k->t=cl->Id()&isid?-q2:q2;
	if (cl->Q2Shower()>=0.0) k->t=cl->Q2Shower();
	if (m_tmap.find(cl->Id())!=m_tmap.end())
	  k->t=cl->Id()&isid?-m_tmap[cl->Id()]:m_tmap[cl->Id()];
	k->stat=3; k->costh=-1.; k->didkin=true;
	k->shower=cl->Stat()==0?2:1;
	if (isid>1 && i<ampl->NIn()) {
	  k->dir=i==1?dir:-dir;
	  if (i==0) k->x=ampl->X1();
	  else k->x=ampl->X2();
	  if (ampl->Next()) {
	    if (i==0) k->z=k->zs=ampl->Next()->X1()/ampl->X1();
	    else k->z=k->zs=ampl->Next()->X2()/ampl->X2();
	  }
 	  k->pt2lcm=ampl->MuF2();
	  if (ampl->Legs().size()==4) {
	    k->qjv=sqrt(4.0*ampl->MuR2());
	    k->part->SetInfo('G');
	  }
	}
	if ((cl->Id()&isid)==0) {
	  if (left) right=k;
	  else left=k;
	}
	if (ampl->Legs().size()>ampl->NIn()+2) {
	  bool found(false);
	  if (cl->Id()&3 && IdCount(cl->Id())>1) k->stat=0;
	  for (std::map<size_t,Knot*>::const_iterator 
		 kit(kmap.begin());kit!=kmap.end();++kit)
	    if (kit->first&cl->Id()) {
	      kit->second->t=kit->second->part->Momentum().Abs2();
	      if (kit->first&isid) {
		if (kit->second->prev==NULL || 
		    kit->second->prev->left==NULL) {
		  if (kit->second->prev) {
		    kit->second->prev->left=k;
		    k->prev=kit->second->prev;
		  }
		  else {
		    kit->second->prev=k;
		    k->right=kit->second;
		  }
		  found=true;
		  break;
		}
	      }
	      else {
		if (kit->second->left==NULL || 
		    kit->second->right==NULL) {
		  k->prev=kit->second;
		  if (!kit->second->left) {
		    kit->second->left=k;
		  }
		  else {
		    kit->second->right=k;
		    SetFSProps(kit->second,kit->second->left,k);
		  }
		  found=true;
		  break;
		}
	      }
	    }
	  if (!found) THROW(fatal_error,"No mother knot");
	}
	kmap[cl->Id()]=k;
      }
    }
    if (ampl->Legs().size()==ampl->NIn()+2) {
      if (left==NULL || right==NULL) THROW(fatal_error,"No FS knots");
      mo->part->SetMomentum(left->part->Momentum()+right->part->Momentum());
      m_cfcms=Poincare(mo->part->Momentum());
      mo->part->SetInfo('f');
      SetFSProps(mo,left,right);
      s=(ampl->Leg(0)->Mom()+ampl->Leg(1)->Mom()).Abs2();
    }
    ampl=ampl->Prev();
  }
  for (std::map<size_t,Knot*>::const_iterator kit(kmap.begin());
       kit!=kmap.end();++kit) kit->second->tmax=tmin;
  p_fintree->CheckMomentumConservation();
  if (p_initrees) {
    p_initrees[0]->CheckMomentumConservation();
    p_initrees[1]->CheckMomentumConservation();
  }
  mo->stat=0;
  return true;
}

void Apacic::SetFSProps
(Knot *const mo,Knot *const left,Knot *const right)
{
  Vec4D pm(m_cfcms*mo->part->Momentum());
  Vec4D pl(m_cfcms*left->part->Momentum());
  Vec4D pr(m_cfcms*right->part->Momentum());
  mo->part->SetStatus(part_status::decayed);
  mo->t=mo->maxpt2=mo->tout=mo->part->Momentum().Abs2();
  mo->costh=-1; mo->didkin=true; 
  mo->thcrit=pl.Theta(pr);
  mo->zs=mo->z=pl[0]/pm[0];
  mo->E2=sqr(pm[0]);
  left->E2=sqr(pl[0]);
  right->E2=sqr(pr[0]);
  left->prev=right->prev=mo;
  mo->left=left; 
  mo->right=right;
  double slt(left->t), srt(right->t); 
  p_finshower->EstablishRelations(mo,left,right);
  if (mo->prev) {
    mo->t=mo->prev->t;
    mo->stat=3;
    mo->shower=2;
  }
  if (slt>=0.0) left->t=slt;
  if (srt>=0.0) right->t=srt;
}

double Apacic::CalculateWeight(Cluster_Amplitude *const ampl)
{
  if (m_fcs>0.0) m_fcs=ampl->MuR2()>0.0?ampl->MuR2():1.0;
  if (m_wmode==1) return PseudoShowerWeight(ampl);
  return 1.0;
}

double Apacic::CouplingWeight(const ATOOLS::Cluster_Leg *const ij,
			      const size_t &oqcd,const double &kt2,
			      const double &kt2r) const
{
  double ckt2(ij->Id()&3?kt2:kt2/4.0);
  double asc((*MODEL::as)(m_fcs>0.0?m_fcs:ckt2));
  double asr((*MODEL::as)(kt2r));
  msg_Debugging()<<((ij->Id()&3)?"is":"fs")<<" as weight "
		 <<" (\\alpha_s("<<sqrt(m_fcs>0.0?m_fcs:ckt2)
		 <<")/\\alpha_s("<<sqrt(kt2r)<<"))^O_{as} = ( "
		 <<asc<<" / "<<asr<<" ) ^ "<<oqcd<<" = "
		 <<pow(asc/asr,oqcd)<<"\n";
  return pow(asc/asr,oqcd);
}

double Apacic::PseudoShowerWeight(Cluster_Amplitude *const ampl)
{
  msg_Debugging()<<METHOD<<"(): {\n";
  // emulate sudakov weight
  msg_Indent();
  if (p_ampl) {
    while (p_ampl->Prev()) p_ampl=p_ampl->Prev();
    p_ampl->Delete();
    p_ampl=NULL;
  }
  // calculate splitting weight
  double wgt(1.0);
  {
    msg_Indent();
    std::map<size_t,Cluster_Leg*> legs;
    Cluster_Amplitude *ref(ampl);
    while (ref->Next()) ref=ref->Next();
    for (size_t i(0);i<ref->Legs().size();++i)
      legs[ref->Leg(i)->Id()]=ref->Leg(i);
    while (ref->Prev()) {
      ref=ref->Prev();
      bool split(false);
      for (size_t i(0);i<ref->Legs().size();++i) {
	size_t idi(ref->Leg(i)->Id());
	if (legs.find(idi)==legs.end()) {
	  for (size_t j(i+1);j<ref->Legs().size();++j) {
	    size_t idj(ref->Leg(j)->Id());
	    if (legs.find(idi+idj)!=legs.end()) {
	      msg_Debugging()<<"split "<<ID(idi+idj)
			     <<" -> "<<ID(idi)<<","<<ID(idj)<<" {\n";
	      msg_Indent();
	      double ckt2(ref->KT2QCD());
	      size_t coqcd(ref->OrderQCD()-ref->Next()->OrderQCD());
	      if (coqcd) wgt*=CouplingWeight
		(legs[idi+idj],coqcd,ckt2,ref->MuR2());
	      legs[idi]=ref->Leg(i);
	      legs[idj]=ref->Leg(j);
	      split=true;
	      break;
	    }
	  }
	  msg_Debugging()<<"}\n";
	  break;
	}
      }
      if (!split) THROW(fatal_error,"Internal error");
    } 
  }
  msg_Debugging()<<"} -> w = "<<wgt<<"\n";
  return wgt;
}

ATOOLS::Cluster_Amplitude *Apacic::GetRealEmissionAmplitude()
{
  return p_ampl;
}

int Apacic::TrialEmission()
{
  return TrialEmission(0);
}

int Apacic::TrialEmission(const int mode)
{
  /*
    stopping the shower works only for e+e- with apacic,
    because the fucking structure is so stupid that it
    does not evolve simultaneously initial and final state.
    thus we abort if a user tries to do it.
  */
  if (m_isron) THROW(fatal_error,"Invalid shower generator");
  int fsos(p_finshower->Sudakov()->OrderingScheme());
  if (mode==0) {
    p_finshower->Sudakov()->SetRBMax(p_rampl->RBMax());
    p_finshower->Sudakov()->SetOrderingScheme(0);
  }
  p_finshower->Sudakov()->SetFCS(m_fcs);
  int res(PerformShowers(1));
  p_finshower->Sudakov()->SetFCS(0.0);
  p_finshower->Sudakov()->SetOrderingScheme(fsos);
  p_finshower->Sudakov()->SetRBMax(1.0);
  if (res!=1) return -1;
  if (p_ampl) {
    while (p_ampl->Prev()) p_ampl=p_ampl->Prev();
    p_ampl->Delete();
  }
  p_ampl = p_rampl->Copy();
  p_ampl=p_ampl->InitNext();
  p_ampl->SetNIn(p_rampl->NIn());
  for (size_t i(0);i<p_ampl->NIn();++i) {
    p_ampl->CreateLeg(p_rampl->Leg(i)->Mom(),p_rampl->Leg(i)->Flav(),
		      p_rampl->Leg(i)->Col(),p_rampl->Leg(i)->Id());
    p_ampl->Legs().back()->SetStat(1);
  }
  p_finshower->ExtractFinalState(p_ampl,p_fintree->GetRoot());
  if (msg_LevelIsDebugging()) {
    msg_Out()<<METHOD<<"(): (R/B)_{max} = "
	     <<p_rampl->RBMax()<<" {\n";
    Cluster_Amplitude *ampl(p_ampl);
    while (ampl) {
      ampl->Print();
      ampl=ampl->Prev();
    }
  }
  msg_Debugging()<<"}\n";
  return p_rampl->Legs().size()<p_ampl->Legs().size();
}

double Apacic::TrialWeight(ATOOLS::Cluster_Amplitude *const ampl)
{
  double wgt(0.0);
  // double as((*MODEL::as)(rpa->gen.CplScale()));
  const Knot *mo(p_finshower->Sudakov()->LastEmission());
  const Knot *d1(mo->left), *d2(mo->right);
  if (mo->part->Flav()==d1->part->Flav()) std::swap<const Knot*>(d1,d2);
  msg_Debugging()<<METHOD<<"(): "<<d1->part->Flav()
		 <<", (R/B)_{max} = "<<p_rampl->RBMax()<<" {\n";
  for (size_t j(0);j<p_ampl->Legs().size();++j) {
    if (p_ampl->Leg(j)->Mom()==d1->part->Momentum()) continue;
    /*
    SK_Map::const_iterator fita(m_sfs.find(p_ampl->Leg(j)->Flav()));
    if (fita==m_sfs.end()) continue;
    std::map<Flavour,Splitting_Kernel*>::const_iterator 
      fitb(fita->second.find(d1->part->Flav()));
    if (fitb==fita->second.end()) continue;
    Splitting_Kernel *k(fitb->second);
    Vec4D pa(p_ampl->Leg(j)->Mom()+d1->part->Momentum());
    double t(pa.Abs2()), z(d1->part->Momentum()[0]/pa[0]);
    double qijk(sqrt(mo->prev->t)), x(2.0*p_ampl->Leg(j)->Mom()[0]/qijk);
    double w(8.0*M_PI*as/t*(*k)(z)*(1.0-z)/x);
    msg_Debugging()<<"  "<<ID(p_ampl->Leg(j)->Id())<<": "<<k->GetA()
		   <<"->("<<k->GetB()<<","<<k->GetC()<<"), \\sqrt{t} = "
		   <<sqrt(t)<<", z = "<<z<<", qijk = "<<qijk
		   <<" -> w = "<<w<<"\n";
    wgt+=w;
    */
  }
  msg_Debugging()<<"} w = "<<wgt<<"\n";
  return wgt*p_rampl->RBMax();
}

double Apacic::CouplingWeight(ATOOLS::Cluster_Amplitude *const ampl)
{
  return 0.0;
}

std::string Apacic::GetKT2(const std::string &jm2) const
{
  return "0.25*"+jm2;
}

namespace PDF {

  DECLARE_GETTER(Apacic_Getter,"Apacic",Shower_Base,Shower_Key);

  Shower_Base *Apacic_Getter::operator()(const Shower_Key &key) const
  {
    return new Apacic(key.p_isr,key.p_model,key.p_read);
  }

  void Apacic_Getter::PrintInfo(std::ostream &str,const size_t width) const
  { 
    str<<"The Apacic shower"; 
  }

}
