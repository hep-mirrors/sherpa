#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "PDF/Main/Jet_Criterion.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/ZAlign.H"
#include "ATOOLS/Org/Exception.H"

#define METS__reject_unordered
// #define CHECK__stepwise

namespace PHASIC {

  struct MCKey {
    size_t m_i, m_j, m_k;
    ATOOLS::Flavour m_fl;
    MCKey(const size_t &i,const size_t &j,const size_t &k,
	  const ATOOLS::Flavour &fl):
      m_i(i),m_j(j), m_k(k), m_fl(fl) {}
    bool operator<(const MCKey &ck) const
    { 
      if (m_i<ck.m_i) return true;
      if (m_i>ck.m_i) return false;
      if (m_j<ck.m_j) return true;
      if (m_j>ck.m_j) return false;
      if (m_k<ck.m_k) return true;
      if (m_k>ck.m_k) return false;
      return m_fl<ck.m_fl;
    }
  };// end of struct MCKey

  struct CS_Params {
    size_t m_i, m_j, m_k, m_oqcd;
    ATOOLS::Flavour m_fl;
    double m_kt2, m_op2, m_mu2, m_z, m_y;
    ATOOLS::Decay_Info *p_dec;
    ATOOLS::Vec4D m_pijt, m_pkt;
    ATOOLS::Poincare_Sequence m_lam;
    CS_Params(const size_t &i,const size_t &j,
	      const size_t &k,const ATOOLS::Flavour &fl):
      m_i(i),m_j(j), m_k(k), m_oqcd(0), m_fl(fl),
      m_kt2(-1.0), m_op2(-std::numeric_limits<double>::max()),
      m_mu2(-1.0), m_z(0.0), m_y(0.0), p_dec(NULL) {}
    bool operator<(const CS_Params &ck) const
    { 
      if (m_i<ck.m_i) return true;
      if (m_i>ck.m_i) return false;
      if (m_j<ck.m_j) return true;
      if (m_j>ck.m_j) return false;
      if (m_k<ck.m_k) return true;
      if (m_k>ck.m_k) return false;
      return m_fl<ck.m_fl;
    }
    void SetParams(const double &kt2,const double &z,const double &y,
		   const ATOOLS::Vec4D &pijt,const ATOOLS::Vec4D &pkt,
		   const ATOOLS::Poincare_Sequence &lam=
		   ATOOLS::Poincare_Sequence())
    { m_mu2=m_kt2=kt2, m_z=z; m_y=y; m_pijt=pijt; m_pkt=pkt; m_lam=lam; }
  };// end of struct CS_Params

  class METS_Scale_Setter: public Scale_Setter_Base {
  private:

    std::string m_muf2tag, m_mur2tag;

    ATOOLS::Algebra_Interpreter m_muf2calc, m_mur2calc;

    Tag_Setter m_tagset;

    ATOOLS::Vec4D_Vector   m_p;
    ATOOLS::Flavour_Vector m_f;

    SP(Color_Integrator) p_ci;

    size_t m_cnt, m_rej, m_mode, m_rproc, m_vmode, m_vproc;
    double m_lfrac, m_aqed, m_wthres;

    ATOOLS::DecayInfo_Vector m_decids;

    static double s_eps, s_kt2max;

    double CoreScale(ATOOLS::Cluster_Amplitude *const ampl);

    bool CheckColors(const ATOOLS::Cluster_Leg *li,
		     const ATOOLS::Cluster_Leg *lj,
		     const ATOOLS::Cluster_Leg *lk,
		     const ATOOLS::Flavour &mo) const;
    ATOOLS::ColorID CombineColors(const ATOOLS::Cluster_Leg *li,
				  const ATOOLS::Cluster_Leg *lj,
				  const ATOOLS::Cluster_Leg *lk,
				  const ATOOLS::Flavour &mo) const;

    double Lam(const double &s,
	       const double &sb,const double &sc) const;

    void KT2(const ATOOLS::Cluster_Leg *li,
	     const ATOOLS::Cluster_Leg *lj,
	     const ATOOLS::Cluster_Leg *lk,CS_Params &cs) const;

    bool Combine(ATOOLS::Cluster_Amplitude &ampl,
		 int i,int j,int k,const CS_Params &cs) const;

    double SetScales(const double &muf2,ATOOLS::Cluster_Amplitude *ampl);

    double CalculateStrict(const ATOOLS::Vec4D_Vector &momenta,const int mode);

  public:

    METS_Scale_Setter(const Scale_Setter_Arguments &args,
		      const int mode=1);

    ~METS_Scale_Setter();

    double CalculateScale(const ATOOLS::Vec4D_Vector &p,const int mode);
    double CalculateMyScale(const ATOOLS::Vec4D_Vector &p,const int mode);

    ATOOLS::Vec4D Momentum(const size_t &i) const;

    void SetScale(const std::string &mu2tag,Tag_Setter &mu2tagset,
		  ATOOLS::Algebra_Interpreter &mu2calc);

  };// end of class METS_Scale_Setter

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Loose_METS_Scale_Setter_Getter,"LOOSE_METS",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *Loose_METS_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new METS_Scale_Setter(args,0);
}

void Loose_METS_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"loose mets scale scheme";
}

DECLARE_GETTER(METS_Scale_Setter_Getter,"METS",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *METS_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new METS_Scale_Setter(args,1);
}

void METS_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"mets scale scheme";
}

DECLARE_GETTER(Strict_METS_Scale_Setter_Getter,"STRICT_METS",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *Strict_METS_Scale_Setter_Getter::
operator()(const Scale_Setter_Arguments &args) const
{
  return new METS_Scale_Setter(args,2);
}

void Strict_METS_Scale_Setter_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"strict mets scale scheme";
}

double METS_Scale_Setter::s_eps=1.0e-6;
double METS_Scale_Setter::s_kt2max=
       sqrt(std::numeric_limits<double>::max());

METS_Scale_Setter::METS_Scale_Setter
(const Scale_Setter_Arguments &args,const int mode):
  Scale_Setter_Base(args), m_tagset(this),
  m_cnt(0), m_rej(0), m_mode(mode), m_lfrac(0.0)
{
  m_p.resize(4);
  size_t pos(args.m_scale.find('{'));
  std::string mur2tag("MU_R2"), muf2tag("MU_F2");
  if (pos!=std::string::npos) {
    muf2tag=args.m_scale.substr(pos+1);
    pos=muf2tag.rfind('}');
    if (pos==std::string::npos)
      THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
    muf2tag=muf2tag.substr(0,pos);
    pos=muf2tag.find("}{");
    if (pos==std::string::npos) {
      mur2tag=muf2tag;
    }
    else {
      mur2tag=muf2tag.substr(pos+2);
      muf2tag=muf2tag.substr(0,pos);
    }
  }
  SetScale(muf2tag,m_tagset,m_muf2calc);
  SetScale(mur2tag,m_tagset,m_mur2calc);
  m_scale.resize(p_proc->NOut());
  SetCouplings();
  m_f=p_proc->Flavours();
  Subprocess_Info info(p_proc->Info().m_ii);
  info.Add(p_proc->Info().m_fi);
  m_decids=info.GetDecayInfos();
  m_aqed=(*MODEL::aqed)(sqr(Flavour(kf_Z).Mass()));
  for (size_t i(0);i<p_proc->NIn();++i) m_f[i]=m_f[i].Bar();
  m_rproc=p_proc->Info().Has(nlo_type::real);
  m_vproc=!p_proc->Parent()->Info().m_fi.NLOType()==nlo_type::lo;
  Data_Reader read(" ",";","!","=");
  if (!read.ReadFromFile(m_vmode,"METS_SCALE_VMODE")) m_vmode=8|2|4;
  else msg_Info()<<METHOD<<"(): Set NLO scale mode "<<m_vmode<<".\n";
  if (!read.ReadFromFile(m_wthres,"METS_WARNING_THRESHOLD")) m_wthres=0.1;
  if (m_vproc && (m_vmode&8)) m_mode=0;
}

METS_Scale_Setter::~METS_Scale_Setter()
{
  for (size_t i(0);i<m_ampls.size();++i) m_ampls[i]->Delete();
}

Vec4D METS_Scale_Setter::Momentum(const size_t &i) const
{
  if (i>m_p.size()) THROW(fatal_error,"Momentum index too large");
  return m_p[i];
}

double METS_Scale_Setter::CalculateStrict
(const Vec4D_Vector &momenta,const int mode)
{
  if (p_caller->Shower()==NULL) THROW(fatal_error,"No shower");
  DEBUG_FUNC(p_caller->Name());
  p_caller->Integrator()->SetMomenta(momenta);
  p_caller->Generator()->SetClusterDefinitions
    (p_caller->Shower()->GetClusterDefinitions());
  Cluster_Amplitude *ampl
    (p_caller->Generator()->ClusterConfiguration(p_caller,m_vproc));
  if (ampl==NULL) {
    msg_Debugging()<<METHOD<<"(): No CSS history for '"
		   <<p_caller->Name()<<"'. Set \\hat{s}.\n";
    ++m_rej;
    double frac(m_rej/(double)m_cnt);
    if (frac>1.25*m_lfrac && m_cnt>5000) {
      m_lfrac=frac;
      if (m_lfrac>m_wthres)
      msg_Error()<<METHOD<<"(): No CSS history for '"
		 <<p_caller->Name()<<"' in >"
		 <<(int(m_lfrac*10000)/100.0)
		 <<"% of calls. Set \\hat{s}."<<std::endl;
    }
    return SetScales((m_p[0]+m_p[1]).Abs2(),NULL);
  }
  Cluster_Amplitude *rampl(ampl);
  while (rampl->Next()) rampl=rampl->Next();
  double muf2(SetScales(rampl->Mu2(),ampl));
  if (p_caller->LookUp()) ampl->Delete();
  else m_ampls.push_back(ampl);
  return muf2;
}

double METS_Scale_Setter::CalculateScale
(const Vec4D_Vector &momenta,const int mode) 
{
  if (mode==0) return CalculateMyScale(momenta,mode);
  if (m_mode==2 || (m_mode==1 && !p_caller->LookUp())) {
    p_caller->Integrator()->SetMomenta(momenta);
    p_caller->Integrator()->SwapInOrder();
    double muf2(CalculateMyScale(p_caller->Integrator()->Momenta(),1));
    p_caller->Integrator()->RestoreInOrder();
    p_cpls->Calculate();
    return muf2;
  }
  return m_scale[stp::fac];
}

double METS_Scale_Setter::CalculateMyScale
(const Vec4D_Vector &momenta,const int mode) 
{
  if (m_escale.size()) {
    m_scale[stp::fac]=m_escale[stp::fac];
    m_scale[stp::ren]=m_escale[stp::ren];
    p_cpls->Calculate();
    return m_scale[stp::fac];    
  }
  ++m_cnt;
  m_p=momenta;
  p_ci=p_proc->Integrator()->ColorIntegrator();
  for (size_t i(0);i<p_proc->NIn();++i) m_p[i]=-m_p[i];
  if (mode==0) {
    while (m_ampls.size()) {
      m_ampls.back()->Delete();
      m_ampls.pop_back();
    }
  }
  if (m_mode==2 || (m_mode==1 && !p_caller->LookUp())) {
    m_scale2=p_proc->Integrator()->ISR()->On() && m_f[0]!=m_f[1];
    return CalculateStrict(momenta,mode);
  }
  DEBUG_FUNC(p_proc->Name());
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  ampl->SetNIn(p_proc->NIn());
  if (p_ci==NULL) {
    for (size_t i(0);i<m_p.size();++i) ampl->CreateLeg(m_p[i],m_f[i]);
  }
  else {
    Int_Vector ci(p_ci->I()), cj(p_ci->J());
    for (size_t i(0);i<m_p.size();++i) 
      ampl->CreateLeg(m_p[i],m_f[i],ColorID(ci[i],cj[i]));
  }
  Single_Process *proc(p_proc->Get<Single_Process>());
  std::vector<std::set<CS_Params> > alltrials(ampl->Legs().size()-4);
  std::vector<bool> ords(ampl->Legs().size()-3,true);
  std::vector<double> ops(ampl->Legs().size()-3,0.0);
  double kt2core(ampl->Legs().size()>4?0.0:CoreScale(ampl));
  ampl->SetOrderQCD(proc->OrderQCD());
  while (ampl->Legs().size()>4) {
    msg_Debugging()<<"Actual = "<<*ampl<<"\n";
    std::set<CS_Params> &trials(alltrials[ampl->Legs().size()-5]);
    size_t iw(0), jw(0), kw(0);
    CS_Params ckw(0,0,0,kf_none);
    msg_Debugging()<<"Weights: {\n";
    for (size_t i(0);i<ampl->Legs().size();++i) {
      msg_Indent();
      Cluster_Leg *li(ampl->Leg(i));
      for (size_t j(Max((size_t)2,i+1));j<ampl->Legs().size();++j) {
	Cluster_Leg *lj(ampl->Leg(j));
	if (!proc->Combinable(li->Id(),lj->Id())) continue;
	Decay_Info *dec(NULL);
	for (size_t l(0);l<m_decids.size();++l)
	  if (m_decids[l].m_id==li->Id()+lj->Id()) {
	  msg_Debugging()<<"cut propagator "<<ID(li->Id()+lj->Id())<<"\n";
	  dec=&m_decids[l];
	  break;
	}
	const Flavour_Vector &cf(proc->CombinedFlavour(li->Id()+lj->Id()));
	for (size_t k(0);k<ampl->Legs().size();++k) {
	  Cluster_Leg *lk(ampl->Leg(k));
	  if (k!=i && k!=j) {
	    for (size_t f(0);f<cf.size();++f) {
	      CS_Params cs(li->Id(),lj->Id(),lk->Id(),cf[f]);
	      cs.p_dec=dec;
	      if (trials.find(cs)!=trials.end()) continue;
	      if (!CheckColors(li,lj,lk,cf[f])) {
		msg_Debugging()<<"Veto colours: "<<cf[f]<<" = "
			       <<ID(cs.m_i)<<" & "<<ID(cs.m_j)
			       <<" <-> "<<ID(cs.m_k)<<"\n";
		trials.insert(cs);
		continue;
	      }
	      KT2(li,lj,lk,cs);
	      if (cs.m_kt2==-1.0) {
		msg_Debugging()<<"Veto kinematics: "<<cf[f]<<" = "
			       <<ID(cs.m_i)<<" & "<<ID(cs.m_j)
			       <<" <-> "<<ID(cs.m_k)<<"\n";
		trials.insert(cs);
		continue;
	      }
	      if (m_vproc) {
		if (m_vmode&2) cs.m_mu2=cs.m_kt2;
		if (m_vmode&4) cs.m_mu2*=4.0;
	      }
	      if (m_rproc && ampl->Prev()==NULL) cs.m_op2=
		1.0/PDF::Qij2(li->Mom(),lj->Mom(),lk->Mom(),
			      kf_gluon,kf_gluon);
	      msg_Debugging()<<ID(cs.m_i)<<" & "<<ID(cs.m_j)<<" <-> "
			     <<ID(cs.m_k)<<" ["<<cf[f]
			     <<"]: "<<cs.m_op2<<" -> ";
	      if (cf[f].Strong() &&
		  li->Flav().Strong() &&
		  lj->Flav().Strong()) {
		// strong clustering, reweight with as
		cs.m_op2*=MODEL::as->BoundedAlphaS(cs.m_kt2);
	      }
	      else {
		// crude: reweight with em coupling
		cs.m_op2*=m_aqed;
		if (cf[f].IsPhoton() ||
		    li->Flav().IsPhoton() ||
		    lj->Flav().IsPhoton()) {
		  // if photon, reweight with charge
		  if (li->Flav().IsPhoton() ||
		      lj->Flav().IsPhoton()) {
		    cs.m_op2*=dabs(cf[f].Charge());}
		  else {
		    cs.m_op2*=dabs(li->Flav().Charge());
		  }
		}
		else if (cf[f].Mass()) {
		  if (cf[f].Width()) {
		    // if resonance, reweight with breit-wigner
		    double s((li->Mom()+lj->Mom()).Abs2());
		    double m2(sqr(cf[f].Mass()));
		    cs.m_op2*=cs.m_kt2/
		      sqrt(sqr(s-m2)+m2*sqr(cf[f].Width()));
		  }
		  else {
		    // if non-resonant, reweight with massive prop
		    cs.m_op2*=cs.m_kt2/(cs.m_kt2+sqr(cf[f].Mass()));
		  }
		}
	      }
	      msg_Debugging()<<cs.m_op2<<"\n";
	      if (cs.m_op2>ckw.m_op2) {
		ckw=cs;
		iw=i;
		jw=j;
		kw=k;
	      }
	    }
	  }
	}
      }
    }
    msg_Debugging()<<"}\n";
    trials.insert(ckw);
    if (iw==0 && jw==0 && kw==0) {
      if (ords[ampl->Legs().size()-5]) {
	msg_Debugging()<<"trying unordered mode at level "
		       <<ampl->Legs().size()-5<<"\n";
	ords[ampl->Legs().size()-5]=false;
	trials.clear();
	continue;
      }
      if (ampl->Prev()==NULL) {
	msg_Debugging()<<METHOD<<"(): No CSS history for '"
		       <<p_proc->Name()<<"'. Set \\hat{s}.\n";
	++m_rej;
	double frac(m_rej/(double)m_cnt);
	if (frac>1.25*m_lfrac && m_cnt>5000) {
	  m_lfrac=frac;
	  if (m_lfrac>m_wthres)
	  msg_Error()<<METHOD<<"(): No CSS history for '"
		     <<p_proc->Name()<<"' in >"
		     <<(int(m_lfrac*10000)/100.0)
		     <<"% of calls. Set \\hat{s}."<<std::endl;
	}
	ampl->Delete();
	return SetScales((m_p[0]+m_p[1]).Abs2(),NULL);
      }
      ampl=ampl->Prev();
      ampl->DeleteNext();
      trials.clear();
      continue;
    }
    msg_Debugging()<<"Cluster "<<ckw.m_fl<<" "
		   <<ID(ckw.m_i)<<" & "<<ID(ckw.m_j)
		   <<" <-> "<<ID(ckw.m_k)
		   <<" => "<<sqrt(ckw.m_kt2)
		   <<" ("<<sqrt(ckw.m_op2)<<") <-> "
		   <<sqrt(ops[ampl->Legs().size()-4])<<"\n";
    if (ops[ampl->Legs().size()-4]>ckw.m_kt2) {
      msg_Debugging()<<"unordered configuration [ ord = "
		     <<ords[ampl->Legs().size()-4]<<" ]\n";
#ifdef METS__reject_unordered
      if (ords[ampl->Legs().size()-4]) continue;
#endif
    }
    ampl->SetKT2(ckw.m_kt2);
    ampl->SetMu2(ckw.m_mu2>0.0?ckw.m_mu2:ckw.m_kt2);
    ampl=ampl->InitNext();
    ampl->CopyFrom(ampl->Prev());
    if (!Combine(*ampl,iw,jw,kw,ckw) ||
	ampl->OrderQCD()<ckw.m_oqcd) {
      msg_Debugging()<<"combine failed\n";
      ampl=ampl->Prev();
      ampl->DeleteNext();
      continue;
    }
    ampl->SetOrderQCD(ampl->OrderQCD()-ckw.m_oqcd);
    if (ckw.p_dec) ampl->Decays().push_back(*ckw.p_dec);
    ops[ampl->Legs().size()-4]=ckw.m_kt2;
#ifdef CHECK__stepwise
    Vec4D psum;
    for (size_t i(0);i<ampl->Legs().size();++i) {
      psum+=ampl->Leg(i)->Mom();
    }
    if (!IsEqual(psum,Vec4D(),1.0e-3)) {
      msg_Error()<<METHOD<<"(): Momentum not conserved\n  in process '"
		 <<p_proc->Name()<<"'\n  \\sum p = "<<psum
		 <<" in step \n"<<*ampl->Prev()<<"\n"<<*ampl<<std::endl;
    }
#endif
    if (ampl->Legs().size()==4) {
      if (!proc->Combinable(ampl->Leg(0)->Id(),ampl->Leg(1)->Id()) &&
	  !proc->Combinable(ampl->Leg(0)->Id(),ampl->Leg(2)->Id()) &&
	  !proc->Combinable(ampl->Leg(0)->Id(),ampl->Leg(3)->Id())) {
	ampl=ampl->Prev();
	ampl->DeleteNext();
	continue;
      }
      if (ampl->Decays().size()!=m_decids.size()) {
	msg_Debugging()<<"unclustered decay\n";
	ampl=ampl->Prev();
	ampl->DeleteNext();
	continue;
      }
      kt2core=CoreScale(ampl);
      if (ops[ampl->Legs().size()-4]>kt2core) {
	msg_Debugging()<<"unordered configuration (core) [ ord = "
		       <<ords[ampl->Legs().size()-4]<<" ]\n";
#ifdef METS__reject_unordered
	if (ords[ampl->Legs().size()-4]) {
	  ampl=ampl->Prev();
	  ampl->DeleteNext();
	  continue;
	}
#endif
      }
    }
  }
  while (ampl->Prev()) ampl=ampl->Prev();
  double muf2(SetScales(kt2core,ampl));
  ampl->Delete();
  return muf2;
}

double METS_Scale_Setter::CoreScale(Cluster_Amplitude *const ampl)
{
  ampl->SetProcs(p_proc);
  double kt2cmin(p_proc->Shower()->GetClusterDefinitions()->CoreScale(ampl));
  ampl->SetKT2(kt2cmin);
  ampl->SetMu2(kt2cmin);
  msg_Debugging()<<"Core = "<<*ampl<<" -> "<<sqrt(kt2cmin)<<"\n";
  return kt2cmin;
}

double METS_Scale_Setter::SetScales(const double &muf2,Cluster_Amplitude *ampl)
{
  double mur2(muf2);
  if (ampl) {
    msg_Debugging()<<"Setting scales {\n";
    mur2=1.0;
    double as(1.0), oqcd(0.0), mum2(1.0);
    for (size_t idx(2);ampl->Next();++idx,ampl=ampl->Next()) {
      m_scale[idx]=Max(ampl->Mu2(),MODEL::as->CutQ2());
      m_scale[idx]=Min(m_scale[idx],sqr(rpa.gen.Ecms()));
      mum2=Min(mum2,m_scale[idx]);
      if (m_rproc && ampl->Prev()==NULL) continue;
      double coqcd(ampl->OrderQCD()-ampl->Next()->OrderQCD());
      if (coqcd>0.0) {
	double cas(MODEL::as->BoundedAlphaS(m_scale[idx]));
	msg_Debugging()<<"  \\mu_{"<<idx<<"} = "<<sqrt(m_scale[idx])
		       <<", as = "<<cas<<", O(QCD) = "<<coqcd<<"\n";
	mur2*=pow(m_scale[idx],coqcd);
	as*=pow(cas,coqcd);
	oqcd+=coqcd;
      }
    }
    if (ampl->OrderQCD()) {
      double mu2(Max(ampl->Mu2(),MODEL::as->CutQ2()));
      mum2=Min(mum2,mu2);
      double cas(MODEL::as->BoundedAlphaS(mu2));
      msg_Debugging()<<"  \\mu_{0} = "<<sqrt(mu2)<<", as = "<<cas
		     <<", O(QCD) = "<<ampl->OrderQCD()<<"\n";
      mur2*=pow(mu2,ampl->OrderQCD());
      as*=pow(cas,ampl->OrderQCD());
      oqcd+=ampl->OrderQCD();
    }
    if (m_vproc && (m_vmode&1)) {
      double cas(MODEL::as->BoundedAlphaS(muf2));
      msg_Debugging()<<"  \\mu_{"<<0<<"} = "<<sqrt(muf2)
		     <<", as = "<<cas<<", O(QCD) = "<<1<<"\n";
      mur2*=pow(muf2,1);
      as*=pow(cas,1);
      oqcd+=1;
    }
    if (oqcd==0.0) mur2=muf2;
    else {
      mur2=pow(mur2,1.0/oqcd);
      as=pow(as,1.0/oqcd);
      mur2=MODEL::as->WDBSolve(as,mum2,rpa.gen.CplScale());
      if (!IsEqual((*MODEL::as)(mur2),as))
	msg_Error()<<METHOD<<"(): Failed to determine \\mu."<<std::endl; 
    }
    msg_Debugging()<<"} -> as = "<<as<<" -> "<<sqrt(mur2)<<"\n";
  }
  m_scale[stp::fac]=muf2;
  m_scale[stp::ren]=mur2;
  msg_Debugging()<<"Core / QCD scale = "<<sqrt(m_scale[stp::fac])
		 <<" / "<<sqrt(m_scale[stp::ren])<<"\n";
  m_scale[stp::fac]=m_muf2calc.Calculate()->Get<double>();
  m_scale[stp::ren]=m_mur2calc.Calculate()->Get<double>();
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(m_scale[stp::fac])<<"\n"
		 <<"  \\mu_r = "<<sqrt(m_scale[stp::ren])<<"\n"
		 <<"} <- "<<p_proc->Name()<<"\n";
  p_cpls->Calculate();
  if (ampl) {
    ampl->SetMuF2(m_scale[stp::fac]);
    ampl->SetMuR2(m_scale[stp::ren]);
    while (ampl->Prev()) {
      ampl=ampl->Prev();
      ampl->SetMuF2(m_scale[stp::fac]);
      ampl->SetMuR2(m_scale[stp::ren]);
    }
  }
  return m_scale[stp::fac];
}

void METS_Scale_Setter::SetScale
(const std::string &mu2tag,Tag_Setter &mu2tagset,Algebra_Interpreter &mu2calc)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<p_caller->Name()<<"' {\n";
  msg_Indent();
  mu2tagset.SetCalculator(&mu2calc);
  mu2calc.SetTagReplacer(&mu2tagset);
  mu2tagset.SetTags(&mu2calc);
  mu2calc.Interprete(mu2tag);
  if (msg_LevelIsDebugging()) mu2calc.PrintEquation();
  msg_Debugging()<<"}\n";
}

double METS_Scale_Setter::Lam
(const double &s,const double &sb,const double &sc) const
{
  return sqr(s-sb-sc)-4.0*sb*sc;
}

void METS_Scale_Setter::KT2
(const Cluster_Leg *li,const Cluster_Leg *lj,
 const Cluster_Leg *lk,CS_Params &cs) const
{
  if ((li->Id()&3)<(lj->Id()&3)) std::swap<const Cluster_Leg*>(li,lj);
  Vec4D pi(li->Mom()), pj(lj->Mom()), pk(lk->Mom());
  double mi2=sqr(li->Flav().Mass()), mj2=sqr(lj->Flav().Mass());
  double mij2=sqr(cs.m_fl.Mass()), mk2=sqr(lk->Flav().Mass());
  if (li->Flav().Strong() && lj->Flav().Strong() &&
      cs.m_fl.Strong()) cs.m_oqcd=1;
  cs.m_op2=1.0/PDF::Qij2(pi,pj,pk,li->Flav(),lj->Flav());
  if ((li->Id()&3)==0) {
    if ((lj->Id()&3)==0) {
      if ((lk->Id()&3)==0) {
	Kin_Args ffp(ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,3));
	double kt2=2.0*(pi*pj)*ffp.m_z*(1.0-ffp.m_z)
	  -sqr(1.0-ffp.m_z)*mi2-sqr(ffp.m_z)*mj2;
	if (ffp.m_stat<0) kt2=-1.0;
 	cs.SetParams(kt2,ffp.m_z,ffp.m_y,ffp.m_pi,ffp.m_pk,ffp.m_lam);
	cs.m_mu2*=p_proc->Shower()->CplFac
	  (li->Flav(),lj->Flav(),lk->Flav(),0,cs.m_oqcd?1:2,kt2);
      }
      else {
	Kin_Args fip(ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,3));
	double kt2=2.0*(pi*pj)*fip.m_z*(1.0-fip.m_z)
	  -sqr(1.0-fip.m_z)*mi2-sqr(fip.m_z)*mj2;
	Vec4D sum(rpa.gen.PBeam(0)+rpa.gen.PBeam(1));
	if (fip.m_pk.PPlus()>sum.PPlus() ||
	    fip.m_pk.PMinus()>sum.PMinus() || fip.m_stat<0) kt2=-1.0;
 	cs.SetParams(kt2,fip.m_z,fip.m_y,fip.m_pi,-fip.m_pk);
	cs.m_mu2*=p_proc->Shower()->CplFac
	  (li->Flav(),lj->Flav(),lk->Flav().Bar(),2,cs.m_oqcd?1:2,kt2);
      }
    }
  }
  else {
    if ((lj->Id()&3)==0) {
      if ((lk->Id()&3)==0) {
	Kin_Args ifp(ClusterIFDipole(mi2,mj2,mij2,mk2,0.0,-pi,pj,pk,pk,3|4));
	double kt2=-2.0*(pi*pj)*(1.0-ifp.m_z)-mj2-sqr(1.0-ifp.m_z)*mi2;
	Vec4D sum(rpa.gen.PBeam(0)+rpa.gen.PBeam(1));
	if (ifp.m_pi.PPlus()>sum.PPlus() ||
	    ifp.m_pi.PMinus()>sum.PMinus() || ifp.m_stat<0) kt2=-1.0;
 	cs.SetParams(kt2,ifp.m_z,ifp.m_y,-ifp.m_pi,ifp.m_pk);
	cs.m_mu2*=p_proc->Shower()->CplFac
	  (li->Flav().Bar(),lj->Flav(),lk->Flav(),1,cs.m_oqcd?1:2,kt2);
      }
      else {
	Kin_Args iip(ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,3));
	double kt2=-2.0*(pi*pj)*(1.0-iip.m_z)-mj2-sqr(1.0-iip.m_z)*mi2;
	Vec4D sum(rpa.gen.PBeam(0)+rpa.gen.PBeam(1));
	if (iip.m_pi.PPlus()>sum.PPlus() ||
	    iip.m_pi.PMinus()>sum.PMinus() || iip.m_stat<0) kt2=-1.0;
 	cs.SetParams(kt2,iip.m_z,iip.m_y,-iip.m_pi,-iip.m_pk,iip.m_lam);
	cs.m_mu2*=p_proc->Shower()->CplFac
	  (li->Flav().Bar(),lj->Flav(),lk->Flav().Bar(),3,cs.m_oqcd?1:2,kt2);
      }
    }
  }
}
  
bool METS_Scale_Setter::Combine(Cluster_Amplitude &ampl,int i,int j,int k,
				const CS_Params &cs) const
{
  if (i>j) std::swap<int>(i,j);
  Cluster_Leg *li(ampl.Leg(i)), *lj(ampl.Leg(j)), *lk(ampl.Leg(k));
  li->SetCol(CombineColors(li,lj,lk,cs.m_fl));
  li->SetFlav(cs.m_fl);
  li->SetMom(cs.m_pijt);
  lk->SetMom(cs.m_pkt);
  if (i<2) {
    for (size_t m(0);m<ampl.Legs().size();++m) {
      if ((int)m==i || (int)m==j || (int)m==k) continue;
      ampl.Leg(m)->SetMom(cs.m_lam*ampl.Leg(m)->Mom());
      ampl.Leg(m)->SetK(0);
    }
  }
  li->SetId(li->Id()+lj->Id());
  li->SetK(lk->Id());
  std::vector<Cluster_Leg*>::iterator lit(ampl.Legs().begin());
  for (int l(0);l<j;++l) ++lit;
  (*lit)->Delete();
  ampl.Legs().erase(lit);
  return true;
}

bool METS_Scale_Setter::CheckColors
(const ATOOLS::Cluster_Leg *li,const ATOOLS::Cluster_Leg *lj,
 const ATOOLS::Cluster_Leg *lk,const ATOOLS::Flavour &mo) const
{
  if (mo.StrongCharge()==8) {
    if (!lk->Flav().Strong()) return false;
  }
  else if (mo.Strong()) {
    if (!(lk->Flav().StrongCharge()==8 ||
	  lk->Flav().StrongCharge()==-mo.StrongCharge())) return false;
  }
  else {
    if (lk->Flav().StrongCharge()==8) return false;
    if (li->Col().m_i==-1 && lj->Col().m_i==-1 &&
	lk->Col().m_i==-1) return true;
    ColorID ci(li->Col()), cj(lj->Col());
    if (ci.m_i==cj.m_j && ci.m_j==0 && cj.m_i==0) return true;
    if (ci.m_j==cj.m_i && ci.m_i==0 && cj.m_j==0) return true;
    return false;
  }
  if (li->Col().m_i==-1 && lj->Col().m_i==-1 &&
      lk->Col().m_i==-1) return true;
  ColorID ci(li->Col()), cj(lj->Col()), ck(lk->Col());
  if (ci.m_i<0 && cj.m_i<0 && ck.m_i<0) return true;
  if (li->Flav().StrongCharge()==3) {
    if (lj->Flav().StrongCharge()==-3) {
      if (lk->Flav().StrongCharge()==0) return true;
      if (ci.m_i==ck.m_j || cj.m_j==ck.m_i ||
	  (ci.m_i==cj.m_j && (ck.m_i>0 || ck.m_j>0))) return true;
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (lk->Flav().StrongCharge()==0) return false;
      if (ci.m_i==cj.m_j && 
	  (cj.m_i==ck.m_j || ck.Singlet())) return true;
      if ((ci.m_i==ck.m_j || ck.Singlet()) && 
	  cj.Singlet()) return true;
    }
    else {
      if (lk->Flav().StrongCharge()==8) return false;
      return true;
    }
  }
  else if (li->Flav().StrongCharge()==-3) {
    if (lj->Flav().StrongCharge()==3) {
      if (lk->Flav().StrongCharge()==0) return true;
      if (ci.m_j==ck.m_i || cj.m_i==ck.m_j ||
	  (ci.m_j==cj.m_i && (ck.m_i>0 || ck.m_j>0))) return true;
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (lk->Flav().StrongCharge()==0) return false;
      if (ci.m_j==cj.m_i && 
	  (cj.m_j==ck.m_i || ck.Singlet())) return true;
      if ((ci.m_j==ck.m_i || ck.Singlet()) && 
	  cj.Singlet()) return true;
    }
    else {
      if (lk->Flav().StrongCharge()==8) return false;
      return true;
    }
  }
  else if (li->Flav().StrongCharge()==8) {
    if (lk->Flav().StrongCharge()==0) return false;
    if (lj->Flav().StrongCharge()==8) {
      if (ci.m_i==cj.m_j && 
	  (ci.m_j==ck.m_i || cj.m_i==ck.m_j ||
	   (ci.m_j==cj.m_i && lk->Flav().StrongCharge()!=8))) 
	return true;
      if (ci.m_j==cj.m_i && 
	  (ci.m_i==ck.m_j || cj.m_j==ck.m_i ||
	   (ci.m_i==cj.m_j && lk->Flav().StrongCharge()!=8)))
	return true;
    }
    else if (lj->Flav().StrongCharge()==3) {
      if (ci.m_j==cj.m_i &&
	  (ci.m_i==ck.m_j || ck.Singlet())) return true;
      if ((cj.m_i==ck.m_j || ck.Singlet()) &&
	  ci.Singlet()) return true;
    }
    else if (lj->Flav().StrongCharge()==-3) {
      if (ci.m_i==cj.m_j &&
	  (ci.m_j==ck.m_i || ck.Singlet())) return true;
      if ((cj.m_j==ck.m_i || ck.Singlet()) &&
	  ci.Singlet()) return true;
    }
    else {
      return false;
    }
  }
  else {
    if (lj->Flav().StrongCharge()==8 ||
	lk->Flav().StrongCharge()==8) {
      return false;
    }
    return true;
  }
  return false;
}

ColorID METS_Scale_Setter::CombineColors
(const Cluster_Leg *li,const Cluster_Leg *lj,const Cluster_Leg *lk,
 const ATOOLS::Flavour &mo) const
{
  ColorID ci(li->Col()), cj(lj->Col()), ck(lk->Col());
  if (ci.m_i==-1 && cj.m_i==-1 && ck.m_i==-1) return ColorID();
  if (!mo.Strong()) return ColorID(0,0);
  if (li->Flav().StrongCharge()==3) {
    if (lj->Flav().StrongCharge()==-3) {
      return ColorID(ci.m_i,cj.m_j);
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (cj.Singlet()) return ColorID(ci.m_i,0);
      return ColorID(cj.m_i,0);
    }
    else {
      return ColorID(ci.m_i,0);
    }
  }
  else if (li->Flav().StrongCharge()==-3) {
    if (lj->Flav().StrongCharge()==3) {
      return ColorID(cj.m_i,ci.m_j);
    }
    else if (lj->Flav().StrongCharge()==8) {
      if (cj.Singlet()) return ColorID(0,ci.m_j);
      return ColorID(0,cj.m_j);
    }
    else {
      return ColorID(0,ci.m_j);
    }
  }
  else if (li->Flav().StrongCharge()==8) {
    if (lj->Flav().StrongCharge()==8) {
      if (ci.m_i==cj.m_j && 
	  (ci.m_j==ck.m_i || cj.m_i==ck.m_j ||
	   (ci.m_j==cj.m_i && lk->Flav().StrongCharge()!=8))) 
	return ColorID(cj.m_i,ci.m_j);
      if (ci.m_j==cj.m_i && 
	  (ci.m_i==ck.m_j || cj.m_j==ck.m_i ||
	   (ci.m_i==cj.m_j && lk->Flav().StrongCharge()!=8)))
	return ColorID(ci.m_i,cj.m_j);
      THROW(fatal_error,"Invalid clustering");
    }
    else if (lj->Flav().StrongCharge()==3) {
      if (ci.Singlet()) return ColorID(cj.m_i,0);
      return ColorID(ci.m_i,0);
    }
    else if (lj->Flav().StrongCharge()==-3) {
      if (ci.Singlet()) return ColorID(0,cj.m_j);
      return ColorID(0,ci.m_j);
    }
    else {
      THROW(fatal_error,"Invalid combination");
    }
  }
  else {
    if (lj->Flav().StrongCharge()==3) {
      return ColorID(cj.m_i,0);
    }
    else if (lj->Flav().StrongCharge()==-3) {
      return ColorID(0,cj.m_j);
    }
    else {
      return ColorID(0,0);
    }
  }
  return ColorID();
}
