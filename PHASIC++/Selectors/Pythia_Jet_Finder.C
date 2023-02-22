#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

namespace PHASIC {

  class Pythia_Jet_Finder: public Selector_Base,
		    public ATOOLS::Tag_Replacer {

    double m_qcut;
    std::string m_cuttag;
    bool m_on;

    ATOOLS::Algebra_Interpreter *p_yccalc;
    ATOOLS::Cluster_Amplitude *p_ampl;

    double Value(Cluster_Amplitude *ampl,int mode)
    {
      DEBUG_FUNC("mode = "<<mode);
      msg_Debugging()<<*ampl<<"\n";
      double q2min(std::numeric_limits<double>::max());
      size_t imin(0), jmin(0), kmin(0);
      Flavour mofl;
      for (size_t i(0);i<ampl->Legs().size();++i) {
	Cluster_Leg *li(ampl->Leg(i));
	for (size_t j(Max(i+1,ampl->NIn()));j<ampl->Legs().size();++j) {
	  Cluster_Leg *lj(ampl->Leg(j));
	  if (j<ampl->NIn()) continue;
	  for (size_t k(0);k<ampl->Legs().size();++k) {
	    if (k==i || k==j) continue;
	    Cluster_Leg *lk(ampl->Leg(k));
	    if (i<ampl->NIn() && k>=ampl->NIn()) continue;
	    if (lk->Flav().Strong() &&
		li->Flav().Strong() && lj->Flav().Strong()) {
	      if (i<ampl->NIn()) li->SetMom(-li->Mom());
	      if (k<ampl->NIn()) lk->SetMom(-lk->Mom());
	      Flavour mofl2(li->Flav().Bar());
	      if (mofl2.IsGluon()) mofl2=lj->Flav().Bar();
	      if (mofl2==lj->Flav()) mofl2=Flavour(kf_gluon);
	      if (i<ampl->NIn() && !ampl->Proc<Process_Base>()->
		  Integrator()->ISR()->PDF(i)->Contains(mofl2)) continue;
	      double q2ijk(pT2pythia(ampl,*li,*lj,*lk,i<ampl->NIn()?-1:1));
	      msg_Debugging()<<"Q_{"<<ID(li->Id())<<ID(lj->Id())
			     <<","<<ID(lk->Id())<<"} = "<<sqrt(q2ijk)<<"\n";
	      if (i<ampl->NIn()) li->SetMom(-li->Mom());
	      if (k<ampl->NIn()) lk->SetMom(-lk->Mom());
	      if (q2ijk<q2min) {
		q2min=q2ijk;
		mofl=Flavour(kf_gluon);
		if (li->Flav().IsGluon()) mofl=lj->Flav();
		if (lj->Flav().IsGluon()) mofl=li->Flav();
		imin=i;
		jmin=j;
		kmin=k;
	      }
	    }
	  }
	}
      }
      if (mode==0 || imin==jmin) return q2min;
      Vec4D_Vector p=Combine(*ampl,imin,jmin,kmin,mofl);
      if (p.empty()) {
	msg_Error()<<METHOD<<"(): Combine failed. Use R configuration."<<std::endl;
	return Value(ampl,0);
      }
      Cluster_Amplitude *bampl(Cluster_Amplitude::New());
      bampl->SetProc(ampl->Proc<void>());
      bampl->SetMS(ampl->MS());
      bampl->SetNIn(ampl->NIn());
      bampl->SetJF(ampl->JF<void>());
      for (int i(0), j(0);i<ampl->Legs().size();++i) {
	if (i==jmin) continue;
	if (i==imin) {
	  bampl->CreateLeg(p[j],mofl,ampl->Leg(i)->Col());
	  bampl->Legs().back()->SetId(ampl->Leg(imin)->Id()|ampl->Leg(jmin)->Id());
	  bampl->Legs().back()->SetK(ampl->Leg(kmin)->Id());	
	}
	else {
	  bampl->CreateLeg(p[j],ampl->Leg(i)->Flav(),ampl->Leg(i)->Col());
	}
	++j;
      }
      q2min=Value(bampl,0);
      bampl->Delete();
      return q2min;
    }

    double pT2pythia(Cluster_Amplitude* ampl, const Cluster_Leg& RadAfterBranch,
                const Cluster_Leg& EmtAfterBranch,
                const Cluster_Leg& RecAfterBranch, int ShowerType){
      // Save type: 1 = FSR pT definition, else ISR definition
      int Type   = ShowerType;
      // Calculate virtuality of splitting
      int sign = (Type==1) ? 1 : -1;
      Vec4D Q(RadAfterBranch.Mom() + sign*EmtAfterBranch.Mom());
      double Qsq = sign * Q.Abs2();
      // Mass term of radiator
      DEBUG_VAR(ampl->MS());
      double m2Rad = ( RadAfterBranch.Flav().Kfcode() >= 4
                   && RadAfterBranch.Flav().Kfcode() < 7)
                   ? ampl->MS()->Mass2(RadAfterBranch.Flav())
                   : 0.;
      // Construct 2->3 variables for FSR
      Vec4D  sum     = RadAfterBranch.Mom() + RecAfterBranch.Mom()
                     + EmtAfterBranch.Mom();
      double m2Dip = sum.Abs2();
      double x1 = 2. * (sum * RadAfterBranch.Mom()) / m2Dip;
      double x3 = 2. * (sum * EmtAfterBranch.Mom()) / m2Dip;
      // Construct momenta of dipole before/after splitting for ISR 
      Vec4D qBR(RadAfterBranch.Mom() - EmtAfterBranch.Mom() + RecAfterBranch.Mom());
      Vec4D qAR(RadAfterBranch.Mom() + RecAfterBranch.Mom()); 
      // Calculate z of splitting, different for FSR and ISR
      double z = (Type==1) ? x1/(x1+x3)
                         : (qBR.Abs2())/( qAR.Abs2());
      // Separation of splitting, different for FSR and ISR
      double pTpyth = (Type==1) ? z*(1.-z) : (1.-z);
      // pT^2 = separation*virtuality
      pTpyth *= (Qsq - sign*m2Rad);
      if(pTpyth < 0.) pTpyth = 0.;
      // Return pT
      return pTpyth;
    }

    ATOOLS::Vec4D_Vector  Combine
    (const Cluster_Amplitude &ampl,int i,int j,int k,const ATOOLS::Flavour &mo)
    {
      Mass_Selector *p_ms=ampl.MS();
      if (i>j) std::swap<int>(i,j);
      Vec4D_Vector after(ampl.Legs().size()-1);
      double mb2(0.0);
      if (i<2) {
	mb2=ampl.Leg(1-i)->Mom().Abs2();
	double mfb2(p_ms->Mass2(ampl.Leg(1-i)->Flav()));
	if ((mfb2==0.0 && IsZero(mb2,1.0e-6)) || IsEqual(mb2,mfb2,1.0e-6)) mb2=mfb2;
      }
      Vec4D pi(ampl.Leg(i)->Mom()), pj(ampl.Leg(j)->Mom());
      Vec4D pk(ampl.Leg(k)->Mom()), pb(i<2?ampl.Leg(1-i)->Mom():Vec4D());
      double mi2=pi.Abs2(), mfi2=p_ms->Mass2(ampl.Leg(i)->Flav());
      double mj2=pj.Abs2(), mfj2=p_ms->Mass2(ampl.Leg(j)->Flav());
      double mk2=pk.Abs2(), mfk2=p_ms->Mass2(ampl.Leg(k)->Flav());
      if ((mfi2==0.0 && IsZero(mi2,1.0e-6)) || IsEqual(mi2,mfi2,1.0e-6)) mi2=mfi2;
      if ((mfj2==0.0 && IsZero(mj2,1.0e-6)) || IsEqual(mj2,mfj2,1.0e-6)) mj2=mfj2;
      if ((mfk2==0.0 && IsZero(mk2,1.0e-6)) || IsEqual(mk2,mfk2,1.0e-6)) mk2=mfk2;
      double mij2=p_ms->Mass2(mo);
      Kin_Args lt;
      if (i>1) {
	if (k>1) lt=ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,2);
	else lt=ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,2);
	if ((k==0 && lt.m_pk[3]<0.0) ||
	    (k==1 && lt.m_pk[3]>0.0) || lt.m_pk[0]<0.0) return Vec4D_Vector();
      }
      else {
	if (k>1) {
	  lt=ClusterIFDipole(mi2,mj2,mij2,mk2,mb2,-pi,pj,pk,-pb,2);
	}
	else lt=ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,2);
	if ((i==0 && lt.m_pi[3]<0.0) ||
	    (i==1 && lt.m_pi[3]>0.0) || lt.m_pi[0]<0.0) return Vec4D_Vector();
      }
      if (lt.m_stat<0) return Vec4D_Vector();
      for (size_t l(0), m(0);m<ampl.Legs().size();++m) {
	if (m==(size_t)j) continue;
	if (m==(size_t)i) after[l]=i>1?lt.m_pi:-lt.m_pi;
	else if (m==(size_t)k) after[l]=k>1?lt.m_pk:-lt.m_pk;
	else after[l]=lt.m_lam*ampl.Leg(m)->Mom();
	++l;
      }
      return after;
    }

  public:

    Pythia_Jet_Finder(): Selector_Base("Pythia_Jet_Finder"),
			 m_qcut(1.), p_yccalc(NULL), p_ampl(NULL)
    {
    }

    Pythia_Jet_Finder(Process_Base *const proc,const std::string &ycut):
      Selector_Base("PythiaJetFinder",proc), m_cuttag(ycut),
      m_on(true), p_yccalc(NULL)
    {
      p_ampl = Cluster_Amplitude::New();
      p_ampl->SetNIn(m_nin);
      for (int i(0);i<m_nin+m_nout;++i)
	p_ampl->CreateLeg(Vec4D(),i<m_nin?p_fl[i].Bar():p_fl[i],ColorID());
      p_ampl->SetJF(this);
      p_ampl->SetMS(p_proc->Generator());
      p_yccalc = new Algebra_Interpreter();
      p_yccalc->SetTagReplacer(this);
      for (int i=0;i<m_n;++i)
	p_yccalc->AddTag("p["+ToString(i)+"]",ToString(Vec4D()));
      p_yccalc->Interprete(m_cuttag);
    }

    ~Pythia_Jet_Finder() 
    {
      if (p_ampl) p_ampl->Delete();
      if (p_yccalc) delete p_yccalc;
    }

    bool Trigger(Selector_List &sl)
    {
      m_pass=true;
      p_ampl->SetProc(p_proc);
      for (size_t i(0);i<sl.size();++i)
	p_ampl->Leg(i)->SetMom((int)i<m_nin?-sl[i].Momentum():sl[i].Momentum());
      m_qcut=p_yccalc->Calculate()->Get<double>();
      if (!m_on) return true;
      msg_Debugging()<<METHOD<<"("<<this<<"): '"<<p_proc->Name()
		     <<"' Q_cut = "<<m_qcut<<(m_on?" {":", off")<<"\n";
      p_ampl->Decays()=p_proc->Info().m_fi.GetDecayInfos();
      double jcv=Value(p_ampl,0);
      bool res=m_pass=jcv>sqr(m_qcut);
      msg_Debugging()<<"} -> "<<res<<"\n";
      return 1-m_sel_log->Hit(!res);
    }

    bool RSTrigger(NLO_subevtlist *const subs)
    {
      for (size_t i(0);i<m_nin+m_nout;++i)
	p_ampl->Leg(i)->SetMom(i<m_nin && subs->back()->p_mom[i][0]>0.0?
			       -subs->back()->p_mom[i]:subs->back()->p_mom[i]);
      m_qcut=p_yccalc->Calculate()->Get<double>();
      if (!m_on) return true;
      int res(0);
      m_pass=0;
      for (size_t n(0);n<subs->size();++n) {
	msg_Debugging()<<METHOD<<"("<<n<<"): '"<<p_proc->Name()
		       <<"' Q_cut = "<<m_qcut<<(m_on?" {":", off")<<"\n";
	{
	  msg_Indent();
	  p_ampl->SetProc(p_proc);
	  if (p_ampl->Legs().size()<(*subs)[n]->m_n)
	    p_ampl->CreateLeg(Vec4D(),Flavour(kf_jet),ColorID());
	  else if (p_ampl->Legs().size()>(*subs)[n]->m_n) {
	    p_ampl->Legs().back()->Delete();
	    p_ampl->Legs().pop_back();
	  }
	  size_t idij((1<<(*subs)[n]->m_i)|(1<<(*subs)[n]->m_j));
	  if ((*subs)[n]->m_i==(*subs)[n]->m_j) idij=0;
	  for (size_t i(0);i<(*subs)[n]->m_n;++i) {
	    p_ampl->Leg(i)->SetFlav
	      ((int)i<m_nin?(*subs)[n]->p_fl[i].Bar():(*subs)[n]->p_fl[i]);
	    p_ampl->Leg(i)->SetMom(i<m_nin && subs->back()->p_mom[i][0]>0.0?
				   -(*subs)[n]->p_mom[i]:(*subs)[n]->p_mom[i]);
	    p_ampl->Leg(i)->SetId((*subs)[n]->p_id[i]);
	    p_ampl->Leg(i)->SetK((*subs)[n]->p_id[i]==idij?
				 (1<<(*subs)[n]->m_k):0);
	  }
	  p_ampl->Decays()=p_proc->Info().m_fi.GetDecayInfos();
	  double jcv=Value(p_ampl,idij?0:1);
	  (*subs)[n]->m_trig=jcv>sqr(m_qcut);
	  if ((*subs)[n]->m_trig) res=m_pass=1;
	}
	msg_Debugging()<<"} -> "<<(*subs)[n]->m_trig<<"\n";
      }
      return 1-m_sel_log->Hit(!res);
    }

    void BuildCuts(Cut_Data *cuts) 
    {
    }

    std::string ReplaceTags(std::string &expr) const
    {
      return p_yccalc->ReplaceTags(expr);
    }

    Term *ReplaceTags(Term *term) const
    {
      term->Set(p_ampl->Leg(term->Id())->Mom());
      return term;
    }

    void AssignId(Term *term)
    {
      term->SetId(ToType<int>
		  (term->Tag().substr
		   (2,term->Tag().length()-3)));
    }

    inline double Qcut() const { return m_qcut; }

    inline void SetOn(const bool on) { m_on=on; }

  };// end of class Pythia_Jet_Finder

}

using namespace PHASIC;

DECLARE_ND_GETTER(Pythia_Jet_Finder,"PythiaJets",Selector_Base,Selector_Key,false);
  
Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Pythia_Jet_Finder>::
operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<1) THROW(critical_error,"Invalid syntax");
  Pythia_Jet_Finder *jf(new Pythia_Jet_Finder(key.p_proc,key[0][0]));
  if (key.front().size()>1 && key[0][1]=="LO" && 
      !(key.front().size()>2 && key[0][2]=="CUT")) 
    jf->SetOn(false);
  return jf;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Pythia_Jet_Finder>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Pythia jet finder"; 
}
