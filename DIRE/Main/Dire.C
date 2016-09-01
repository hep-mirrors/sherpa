#include "PDF/Main/Shower_Base.H"

#include "DIRE/Shower/Shower.H"
#include "DIRE/Shower/Cluster.H"
#include "DIRE/Main/Color_Setter.H"
#include "DIRE/Tools/Amplitude.H"
#include "ATOOLS/Phys/Blob_List.H"

namespace DIRE {

  class Dire: public PDF::Shower_Base {
  private:

    Shower *p_shower;

    Cluster      *p_clus;
    Color_Setter *p_cs;

    Amplitude_Vector m_ampls;

    ATOOLS::Mass_Selector *p_ms;

    int    m_reco, m_wcheck;
    double m_maxweight;

    void RecoCheck(Amplitude *const a,int swap) const;

    void SetColors(ATOOLS::Cluster_Amplitude *ampl) const;

    Amplitude *Convert(ATOOLS::Cluster_Amplitude *const campl,
		       std::map<ATOOLS::Cluster_Leg*,Parton*> &lmap);

    void ExtractParton(ATOOLS::Blob *const bl,Parton *const p);

  public:

    Dire(const PDF::Shower_Key &key);

    ~Dire();

    int  PerformShowers();
    int  PerformDecayShowers();

    bool ExtractPartons(ATOOLS::Blob_List *const bl);
    void CleanUp();

    PDF::Cluster_Definitions_Base *GetClusterDefinitions();

    bool PrepareShower(ATOOLS::Cluster_Amplitude *const ampl,
		       const bool & soft=false);

    double CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
		  const ATOOLS::Flavour &flk,const int type,
		  const int cpl,const double &mu2) const;

  };// end of class Dire

}// end of namespace DIRE

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"

#include <algorithm>

using namespace DIRE;
using namespace PDF;
using namespace ATOOLS;

Dire::Dire(const Shower_Key &key):
  Shower_Base("Dire"), p_cs(NULL), p_ms(NULL),
  m_maxweight(1.0)
{
  p_shower = new Shower();
  p_clus = new Cluster(p_shower);
  p_shower->Init(key.p_model,key.p_isr,key.p_read);
  int csmode=key.p_read->GetValue<int>("CSS_CSMODE",0);
  if (csmode) p_cs = new Color_Setter(csmode);
  m_reco=key.p_read->GetValue<int>("CSS_RECO_CHECK",0);
  if (m_reco) msg_Info()<<METHOD<<"(): Reco check "<<m_reco<<"\n";
  m_wcheck=key.p_read->GetValue<int>("CSS_WEIGHT_CHECK",0);
  if (m_wcheck) msg_Info()<<METHOD<<"(): Weight check "<<m_wcheck<<"\n";
}

Dire::~Dire()
{
  if (p_cs) delete p_cs;
  delete p_clus;
  delete p_shower;
}

int Dire::PerformShowers()
{
  DEBUG_FUNC(this);
  m_weight=1.0;
  unsigned int nem=0;
  for (Amplitude_Vector::const_iterator
	 it(m_ampls.begin());it!=m_ampls.end();++it) {
    int stat(p_shower->Evolve(**it,m_weight,nem));
    m_weight*=p_shower->GetWeight();
    if (stat!=1) return stat;
  }
  if (m_wcheck && dabs(m_weight)>m_maxweight) {
    m_maxweight=dabs(m_weight);
    std::string rname="dire.random."+rpa->gen.Variable("RNG_SEED")+".dat";
    if (ATOOLS::msg->LogFile()!="")
      rname=ATOOLS::msg->LogFile()+"."+rname;
    ATOOLS::ran->WriteOutSavedStatus(rname.c_str());
    std::ofstream outstream(rname.c_str(),std::fstream::app);
    outstream<<std::endl;
    outstream<<"# Wrote status for weight="<<m_weight<<" in event "
	     <<rpa->gen.NumberOfGeneratedEvents()+1<<std::endl;
    outstream.close();
  }
  return 1;
}

int Dire::PerformDecayShowers()
{
  DEBUG_FUNC(this);
  return PerformShowers();
}

void Dire::ExtractParton(Blob *const bl,Parton *const p)
{
  Particle *sp=p->Beam()?
    new Particle(-1,p->Flav().Bar(),-p->Mom(),'I'):
    new Particle(-1,p->Flav(),p->Mom(),'F');
  sp->SetNumber(0);
  sp->SetFinalMass(p_ms->Mass(p->Flav()));
  if (p->Beam()==0) {
    sp->SetFlow(1,p->Col().m_i);
    sp->SetFlow(2,p->Col().m_j);
    bl->AddToOutParticles(sp);
  }
  else {
    sp->SetFlow(1,p->Col().m_j);
    sp->SetFlow(2,p->Col().m_i);
    sp->SetBeam(p->Beam()-1);
    bl->AddToInParticles(sp);
  } 
}

bool Dire::ExtractPartons(Blob_List *const bl)
{
  Blob *b(bl->FindLast(btp::Shower));
  if (b==NULL) THROW(fatal_error,"No Shower blob");
  b->SetTypeSpec("DIRE");
  for (int i=0;i<b->NInP();++i) 
    b->InParticle(i)->SetStatus(part_status::decayed);
  for (int i=0;i<b->NOutP();++i) 
    b->OutParticle(i)->SetStatus(part_status::decayed);
  b->SetStatus(blob_status::needs_beams |
	       blob_status::needs_hadronization);
  bool nois(b->NOutP()==0);
  for (Amplitude_Vector::const_iterator
	 it(m_ampls.begin());it!=m_ampls.end();++it)
    for (Amplitude::const_iterator
	   pit((*it)->begin());pit!=(*it)->end();++pit) {
      if ((*pit)->Beam()&&nois) continue;
      if ((*pit)->Out(0)==NULL) ExtractParton(b,*pit);
    }
  return true;
}

void Dire::CleanUp()
{
  for (Amplitude_Vector::const_iterator it(m_ampls.begin());
       it!=m_ampls.end();++it) delete *it;
  m_ampls.clear();
}

Cluster_Definitions_Base *Dire::GetClusterDefinitions()
{
  return p_clus;
}

Amplitude *Dire::Convert
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
    if (i<campl->NIn()) p->SetBeam(cl->Id()&3);
    lmap[cl]=p;
  }
  msg_Debugging()<<*ampl<<"\n";
  return ampl;
}

bool Dire::PrepareShower
(Cluster_Amplitude *const ampl,const bool &soft)
{
  DEBUG_FUNC(soft);
  p_ms=ampl->MS();
  p_shower->SetMS(p_ms);
  Cluster_Amplitude *campl(ampl);
  while (campl->Next()) campl=campl->Next();
  SetColors(campl);
  double Q2(campl->MuQ2());
  std::map<Cluster_Leg*,Parton*> lmap;
  for (;campl;campl=campl->Prev()) {
    Amplitude *ampl(Convert(campl,lmap));
    m_ampls.push_back(ampl);
    if (campl->NIn()+campl->Leg(2)->NMax()==
	campl->Legs().size()) ampl->SetJF(NULL);
    Cluster_Amplitude *lampl(campl->Next());
    if (lampl) {
      int ic=-1,jc=-1,kc=-1;
      Cluster_Leg *lij(NULL);
      Cluster_Leg *nl(campl->IdLeg(campl->IdNew()));
      for (size_t i(0);i<lampl->Legs().size()&&lij==NULL;++i)
	if (lampl->Leg(i)->K()) lij=lampl->Leg(i);
      if (lij==NULL) THROW(fatal_error,"Invalid PS tree");
      for (size_t i(0);i<lampl->Legs().size();++i) {
	Cluster_Leg *cl(lampl->Leg(i));
	Parton *cp(lmap[cl]);
	for (size_t j(0);j<campl->Legs().size();++j) {
	  Cluster_Leg *dl(campl->Leg(j));
	  if (cl->Id()&dl->Id()) {
	    if (dl->Id()==lij->K()) kc=j;
	    Parton *dp(lmap[dl]);
	    if (cp->Out(0)) {
	      if (cp->Out(1)) 
		THROW(fatal_error,"Invalid PS tree");
	      if (cl==lij) (dl==nl?jc:ic)=j;
	      cp->SetOut(1,dp);
	      dp->SetIn(cp);
	    }
	    else {
	      if (cl==lij) (dl==nl?jc:ic)=j;
	      cp->SetOut(0,dp);
	      dp->SetIn(cp);
	    }
	  }
	}
      }
      if (ic<0 || jc<0 || kc<0)
	THROW(fatal_error,"Invalid PS tree");
      double ws, mu2;
      int flip(jc<ic), swap(jc<campl->NIn() && flip);
      if (swap) std::swap<int>(ic,jc);
      int type((ic<campl->NIn()?1:0)|(kc<campl->NIn()?2:0));
      Splitting s=p_clus->KT2
	     (campl->Leg(ic),campl->Leg(jc),campl->Leg(kc),
	      lij->Flav(),campl->Kin(),type,1|(swap?2:0),ws,mu2);
      s.p_s=lmap[lampl->IdLeg(lij->K())];
      s.p_c=lmap[lij];
      (*----m_ampls.end())->SetSplit(s);
      if (!flip || swap) RecoCheck(*----m_ampls.end(),swap);
    }
  }
  m_ampls.front()->SetT(Q2);
  return true;
}

void Dire::SetColors(ATOOLS::Cluster_Amplitude *ampl) const
{
  if (p_cs==NULL || !p_cs->SetColors(ampl)) {
    std::vector<int> tids, atids;
    for (size_t i(0);i<ampl->Legs().size();++i)
      if (ampl->Leg(i)->Flav().StrongCharge()>0) {
	tids.push_back(i);
	if (ampl->Leg(i)->Flav().StrongCharge()==8)
	  atids.push_back(i);
      }
      else if (ampl->Leg(i)->Flav().StrongCharge()<0) {
	atids.push_back(i);
      }
    while (true) {
      std::random_shuffle(atids.begin(),atids.end(),*ran);
      size_t i(0);
      for (;i<atids.size();++i) if (atids[i]==tids[i]) break;
      if (i==atids.size()) break;
    }
    for (size_t i(0);i<tids.size();++i) {
      int cl(Flow::Counter());
      ampl->Leg(tids[i])->SetCol(ColorID(cl,ampl->Leg(tids[i])->Col().m_j));
      ampl->Leg(atids[i])->SetCol(ColorID(ampl->Leg(atids[i])->Col().m_i,cl));
    }
  }
  for (Cluster_Amplitude *campl(ampl->Prev());campl;campl=campl->Prev()) {
    Cluster_Amplitude *next(campl->Next()); 
    Cluster_Leg *lij=NULL;
    for (size_t i(0);i<next->Legs().size();++i)
      if (next->Leg(i)->K()) {
	lij=next->Leg(i);
	break;
      }
    if (lij==NULL) THROW(fatal_error,"Invalid amplitude");
    Cluster_Leg *li=NULL, *lj=NULL;
    for (size_t i(0);i<campl->Legs().size();++i) {
      if (campl->Leg(i)->Id()&lij->Id()) {
	if (li==NULL) li=campl->Leg(i);
	else if (lj==NULL) lj=campl->Leg(i);
	else THROW(fatal_error,"Invalid splitting");
      }
      else {
	campl->Leg(i)->SetCol(next->IdLeg(campl->Leg(i)->Id())->Col());
      }
    }
    Cluster_Amplitude::SetColours(lij,li,lj);
  }
}

double Dire::CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
		    const ATOOLS::Flavour &flk,const int type,
		    const int cpl,const double &mu2) const
{
  return 1.0;
}

void Dire::RecoCheck(Amplitude *const a,int swap) const
{
  if (!(m_reco&1) || a->Split().p_c==NULL) return;
  DEBUG_FUNC(a);
  Amplitude *next(a->Split().p_c->Out(0)->Ampl());
  int ic=-1, jc=-1, kc=-1;
  Vec4D pi, pj, pk;
  for (size_t i(0);i<next->size();++i) {
    if ((*next)[i]==a->Split().p_c->Out(0)) { ic=i; pi=(*next)[i]->Mom(); }
    if ((*next)[i]==a->Split().p_c->Out(1)) { jc=i; pj=(*next)[i]->Mom(); }
    if ((*next)[i]==a->Split().p_s->Out(0)) { kc=i; pk=(*next)[i]->Mom(); }
  }
  a->Construct();
  Cluster_Amplitude *ampl(next->GetAmplitude());
  double ws, mu2;
  Splitting s=p_clus->KT2
    (ampl->Leg(ic),ampl->Leg(jc),ampl->Leg(kc),
     a->Split().p_c->Flav(),a->Split().m_kin,
     a->Split().m_type,1|(swap?2:0),ws,mu2);
  ampl->Delete();
  std::cout.precision(12);
  msg_Debugging()<<"New reco params: t = "<<s.m_t
		 <<", z = "<<s.m_z<<", phi = "<<s.m_phi<<"\n";
  msg_Debugging()<<"            vs.: t = "<<a->Split().m_t<<", z = "
		 <<a->Split().m_z<<", phi = "<<a->Split().m_phi
		 <<", kin = "<<a->Split().m_kin<<"\n";
  if (!IsEqual(s.m_t,a->Split().m_t,1.0e-6) || 
      !IsEqual(s.m_z,a->Split().m_z,1.0e-6) || 
      !IsEqual(s.m_phi,a->Split().m_phi,1.0e-6) ||
      !IsEqual(pi,ampl->Leg(ic)->Mom(),1.0e-6) || 
      !IsEqual(pj,ampl->Leg(jc)->Mom(),1.0e-6) || 
      !IsEqual(pk,ampl->Leg(kc)->Mom(),1.0e-6)) {
    msg_Error()<<"Faulty reco params: t = "<<s.m_t
	       <<", z = "<<s.m_z<<", phi = "<<s.m_phi<<"\n";
    msg_Error()<<"               vs.: t = "<<a->Split().m_t<<", z = "
	       <<a->Split().m_z<<", phi = "<<a->Split().m_phi
	       <<", kin = "<<a->Split().m_kin<<"\n\n";
    msg_Error()<<"  "<<pi<<" "<<pj<<" "<<pk<<"\n";
    msg_Error()<<"  "<<ampl->Leg(ic)->Mom()
	       <<" "<<ampl->Leg(jc)->Mom()
	       <<" "<<ampl->Leg(kc)->Mom()<<"\n";
    if (m_reco&2) abort();
  }
}

DECLARE_GETTER(Dire,"Dire",Shower_Base,Shower_Key);

Shower_Base *Getter<Shower_Base,Shower_Key,Dire>::
operator()(const Shower_Key &key) const
{
  return new Dire(key);
}

void Getter<Shower_Base,Shower_Key,Dire>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Dipole Parton Shower"; 
}
