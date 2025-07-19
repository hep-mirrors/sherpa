#include "MALARIC/Main/MCatNLO.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace MALARIC;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

MCatNLO::MCatNLO(const NLOMC_Key &key):
  NLOMC_Base("Alaric"), p_ms(NULL)
{
  m_subtype=subscheme::Alaric;
  p_mcatnlo = new Shower();
  p_gamma = new Gamma(this,p_mcatnlo);
  p_mcatnlo->SetGamma(p_gamma);
  p_mcatnlo->Init(key.p_model,key.p_isr);
  p_recoil = p_mcatnlo->Recoil();
  Settings& s = Settings::GetMainSettings();
  m_psmode=s["MC@NLO"]["PSMODE"].Get<int>();
  m_wcheck=s["MC@NLO"]["WEIGHT_CHECK"].Get<int>();
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
  m_weightsmap.Clear();
  m_weightsmap["Sudakov"] = Weights {Variations_Type::qcd};
  m_weightsmap["QCUT"] = Weights {Variations_Type::qcut};
  CleanUp();
  PrepareShower(ampl);
  if (p_rampl->NLO()&4) return 1;
  unsigned int nem=0;
  int stat(p_mcatnlo->Evolve(*m_ampls.back(),nem));
  m_weightsmap["Sudakov"] *= p_mcatnlo->GetWeightsMap().at("Sudakov");
  m_weightsmap["QCUT"] *= p_mcatnlo->GetWeightsMap().at("QCUT");
  if (m_wcheck && dabs(m_weightsmap.Nominal())>m_maxweight) {
    m_maxweight=dabs(m_weightsmap.Nominal());
    std::string rname="alaricnlo.random.dat";
    if (ATOOLS::msg->LogFile()!="")
      rname=ATOOLS::msg->LogFile()+"."+rname;
    ATOOLS::ran->WriteOutSavedStatus(rname.c_str());
    std::ofstream outstream(rname.c_str(),std::fstream::app);
    outstream<<std::endl;
    outstream<<"# Wrote status for weight="<<m_weightsmap.Nominal()<<" in event "
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
  ampl->SetFlag(p_mcatnlo->LastSplitting().p_sk->Id());
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
  int ntrip(0), noct(0);
  const IDip_Set &iinfo((*p_rampl->IInfo<StringIDipSet_Map>())[pname]);
  for (int nf(0);nf<2;++nf) {
  for (size_t i(0);i<m_ampls.back()->size();++i) {
    Parton *c((*m_ampls.back())[i]);
    msg_Debugging()<<"spectators for "
		   <<c->Flav()<<"("<<c->Id()<<") ... ";
    for (size_t j(0);j<m_ampls.back()->size();++j) {
      Parton *s((*m_ampls.back())[j]);
      if (iinfo.find(IDip_ID(c->Id()-1,s->Id()-1))!=iinfo.end()) {
	msg_Debugging()<<s->Flav()<<"("<<s->Id()<<") ";
	c->S().push_back(s);
	if (c->Flav().StrongCharge()==8) ++noct;
	else ++ntrip;
      }
    }
    msg_Debugging()<<"-> "<<c->S().size()<<" dipole(s)\n";
  }
  if (!(ntrip==0 && noct==2)) break;
  else msg_Debugging()<<"gg special, fill again\n";
  }
  if (ampl->NIn()+ampl->Leg(2)->NMax()==
      ampl->Legs().size()+1) m_ampls.back()->SetJF(NULL);
  Cluster_Amplitude *campl(ampl);
  while (campl->Next()) campl=campl->Next();
  return true;
}

double MCatNLO::KT2(const ATOOLS::NLO_subevt &sub,
		    const PHASIC::Kin_Args *ka,
		    const PHASIC::Ant_Args *aa)
{
  if (sub.m_ijt>=2) {
    if (sub.m_type==spt::splittingtype::soft) {
      double kt2;
      if ((p_mcatnlo->EvolScheme(0)&255)==0) {
	kt2=2.*(aa->m_pi*aa->m_pj)*
	  (aa->m_pj*aa->m_n)/(aa->m_pi*aa->m_n);
      }
      else if ((p_mcatnlo->EvolScheme(0)&255)==1) {
	kt2=2.*(aa->m_pi*aa->m_pj)*
	  (aa->m_pj*aa->m_K)/(aa->m_pi*aa->m_K);
      }
      else {
	THROW(not_implemented,"Unknown evolution scheme");
      }
      return kt2;
    }
    else {
      double Q2((ka->m_pi+ka->m_pk).Abs2());
      double yijk(ka->m_y), zi(ka->m_z), kt2;
      if ((p_mcatnlo->EvolScheme(0)>>8)==0) {
	kt2=Q2*yijk*zi*(1.0-zi);
      }
      else if ((p_mcatnlo->EvolScheme(0)>>8)==1) {
	kt2=Q2*yijk;
      }
      else {
	THROW(not_implemented,"Unknown evolution scheme");
      }
      return kt2;
    }
  }
  THROW(fatal_error,"Implement me");
}

DECLARE_GETTER(MALARIC::MCatNLO,"MC@NLO_Alaric",PDF::NLOMC_Base,PDF::NLOMC_Key);

NLOMC_Base *ATOOLS::Getter<PDF::NLOMC_Base,PDF::NLOMC_Key,MALARIC::MCatNLO>::
operator()(const NLOMC_Key &key) const
{
  return new MCatNLO(key);
}

void ATOOLS::Getter<PDF::NLOMC_Base,PDF::NLOMC_Key,MALARIC::MCatNLO>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Alaric MC@NLO generator"; 
}
