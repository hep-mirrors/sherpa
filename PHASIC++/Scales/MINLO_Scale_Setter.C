#include "PHASIC++/Scales/MINLO_Scale_Setter.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(MINLO_Scale_Setter,"MINLO",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,MINLO_Scale_Setter>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new MINLO_Scale_Setter(args);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
		    MINLO_Scale_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"minlo scale scheme";
}

MINLO_Scale_Setter::MINLO_Scale_Setter
(const Scale_Setter_Arguments &args):
  Scale_Setter_Base(args), m_tagset(this), p_ampl(NULL)
{
  std::string tag(args.m_scale), core;
  m_nproc=p_proc->Info().Has(nlo_type::real) ||
    p_proc->Info().Has(nlo_type::rsub) ||
    p_proc->Info().Has(nlo_type::loop) ||
    p_proc->Info().Has(nlo_type::vsub);
  p_isr=p_proc->Integrator()->ISR();
  size_t pos(tag.find('['));
  if (pos!=std::string::npos) {
    tag=tag.substr(pos+1);
    pos=tag.find(']');
    if (pos==std::string::npos) 
      THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
    core=tag.substr(0,pos);
    tag=tag.substr(pos+1);
  }
  if (tag.find('{')==std::string::npos &&
      tag.find('}')==std::string::npos) tag+="{MU_F2}{MU_R2}{MU_Q2}";
  while (true) {
    size_t pos(tag.find('{'));
    if (pos==std::string::npos) {
      if (!m_calcs.empty()) break;
      else { THROW(fatal_error,"Invalid scale '"+args.m_scale+"'"); }
    }
    tag=tag.substr(pos+1);
    pos=tag.find('}');
    if (pos==std::string::npos) 
      THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
    std::string ctag(tag.substr(0,pos));
    tag=tag.substr(pos+1);
    m_calcs.push_back(new Algebra_Interpreter());
    m_calcs.back()->AddFunction(MODEL::as->GetAIGMeanFunction());
    m_calcs.back()->SetTagReplacer(&m_tagset);
    if (m_calcs.size()==1) m_tagset.SetCalculator(m_calcs.back());
    SetScale(ctag,*m_calcs.back());
  }
  m_scale.resize(Max(m_scale.size(),m_calcs.size()));
  SetCouplings();
  Data_Reader read(" ",";","!","=");
  if (!read.ReadFromFile(m_noutmin,"MINLO_NOUT_MIN")) m_noutmin=2;
  if (!read.ReadFromFile(m_cmode,"MINLO_CLUSTER_MODE")) m_cmode=1;
  if (!read.ReadFromFile(m_hqmode,"MINLO_HQ_MODE")) m_hqmode=1;
  if (!read.ReadFromFile(m_order,"MINLO_FORCE_ORDER")) m_order=0;
  if (!read.ReadFromFile(m_orderrs,"MINLO_ORDER_RS")) m_orderrs=1;
  if (!read.ReadFromFile(m_usecomb,"MINLO_USE_COMBINABLE")) m_usecomb=0;
  if (!read.ReadFromFile(m_usepdfinfo,"MINLO_USE_PDFINFO")) m_usepdfinfo=1;
  if (!read.ReadFromFile(m_nlocpl,"MINLO_NLO_COUPLING_MODE")) m_nlocpl=1;
  if (!read.ReadFromFile(m_mufmode,"MINLO_MUF_VARIATION_MODE")) m_mufmode=1;
  if (!read.ReadFromFile(m_dr,"MINLO_DELTA_R")) m_dr=0.4;
  if (!read.ReadFromFile(m_muf2min,"MINLO_MUF2_MIN"))
    m_muf2min=p_isr->PDF(0)->Q2Min();
  if (core=="" && !read.ReadFromFile(core,"CORE_SCALE")) core="VAR{H_TM2/4}";
  p_core=Core_Scale_Getter::GetObject(core,Core_Scale_Arguments(p_proc,core));
  if (p_core==NULL) THROW(fatal_error,"Invalid core scale '"+core+"'");
  if (!read.ReadFromFile(m_nfgsplit,"DIPOLE_NF_GSPLIT"))
    m_nfgsplit=Flavour(kf_jet).Size()/2;
  else msg_Tracking()<<METHOD<<"(): Set dipole N_f="<<m_nfgsplit<<"\n.";
  m_rsf=ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
  if (m_rsf!=1.0) msg_Debugging()<<METHOD<<
		    "(): Renormalization scale factor "<<sqrt(m_rsf)<<"\n";
}

MINLO_Scale_Setter::~MINLO_Scale_Setter()
{
  for (size_t i(0);i<m_calcs.size();++i) delete m_calcs[i];
  delete p_core;
}

double MINLO_Scale_Setter::Calculate(const Vec4D_Vector &momenta,const size_t &mode) 
{
  m_p=momenta;
  for (size_t i(0);i<p_proc->NIn();++i) m_p[i]=-m_p[i];
  if (p_ampl) p_ampl->Delete();
  m_rproc=p_caller->Info().Has(nlo_type::real);
  m_vproc=p_caller->Info().Has(nlo_type::loop) || p_caller->Info().Has(nlo_type::vsub);
  m_f=p_caller->Flavours();
  for (size_t i(0);i<p_caller->NIn();++i) m_f[i]=m_f[i].Bar();
  DEBUG_FUNC(p_proc->Name()<<" from "<<p_caller->Name()<<", R="<<m_rproc<<", V="<<m_vproc);
  Cluster_Amplitude *ampl(p_ampl = Cluster_Amplitude::New());
  ampl->SetNIn(p_proc->NIn());
  for (size_t i(0);i<m_p.size();++i) ampl->CreateLeg(m_p[i],m_f[i]);
  Single_Process *proc(p_proc->Get<Single_Process>());
  std::vector<std::set<MCS_Params> > alltrials(ampl->Legs().size()-(m_nin+m_noutmin-1));
  std::vector<std::vector<std::pair<size_t,double> > >
    ops(ampl->Legs().size()-(m_nin+m_noutmin-1));
  ops[ampl->Legs().size()-(m_nin+m_noutmin)].push_back
    (std::pair<size_t,double>((1<<ampl->Legs().size())-1,0.0));
  if (ampl->Legs().size()==m_nin+m_noutmin) CoreScale(ampl);
  ampl->SetOrderQCD(p_caller->MaxOrder(0));
  while (ampl->Legs().size()>=m_nin+m_noutmin) {
    msg_Debugging()<<"Actual = "<<*ampl<<"\n";
    std::set<MCS_Params> &trials(alltrials[ampl->Legs().size()-(m_nin+m_noutmin)]);
    MCS_Params ckw(0,0,0,kf_none);
    msg_Debugging()<<"Weights: {\n";
    for (size_t i(0);i<ampl->Legs().size();++i) {
      msg_Indent();
      Cluster_Leg *li(ampl->Leg(i));
      for (size_t j(Max((size_t)2,i+1));j<ampl->Legs().size();++j) {
	Cluster_Leg *lj(ampl->Leg(j));
	Flavour_Vector cf;
	if (m_usecomb) {
	  if (!proc->Combinable(li->Id(),lj->Id())) continue;
	  cf=proc->CombinedFlavour(li->Id()+lj->Id());
	}
	else {
	  if (!li->Flav().Strong() || !lj->Flav().Strong()) continue;
	  if (m_hqmode==1 && (li->Flav().Mass() || lj->Flav().Mass())) continue;
	  if (li->Flav().IsGluon()) cf.push_back(lj->Flav());
	  else if (lj->Flav().IsGluon()) cf.push_back(li->Flav());
	  else if (li->Flav()==lj->Flav().Bar()) cf.push_back(kf_gluon);
	}
	for (size_t k(0);k<ampl->Legs().size();++k) {
	  Cluster_Leg *lk(ampl->Leg(k));
	  if (k!=i && k!=j) {
	    for (size_t f(0);f<cf.size();++f) {
	      if (!cf[f].Strong() || !lk->Flav().Strong() ||
		  !li->Flav().Strong() || !lj->Flav().Strong()) continue;
	      if (i<ampl->NIn() && m_usepdfinfo) {
		if (p_isr->PDF(i)==NULL) THROW(fatal_error,"No PDF");
		if (!p_isr->PDF(i)->Contains(cf[f])) continue;
	      }
	      MCS_Params cs(i,j,k,cf[f]);
	      if (trials.find(cs)!=trials.end()) continue;
	      if (cf[f].IsGluon() &&
		  !li->Flav().IsGluon() &&
		  li->Flav().Kfcode()>m_nfgsplit &&
		  !lj->Flav().IsGluon() &&
		  lj->Flav().Kfcode()>m_nfgsplit) {
		msg_Debugging()<<"Veto flavour: "<<cf[f]<<" = "
			       <<ID(li->Id())<<" & "<<ID(lj->Id())
			       <<" <-> "<<ID(lk->Id())<<"\n";
		trials.insert(cs);
		continue;
	      }
	      KT2(li,lj,lk,cs);
	      if (cs.m_kt2==-1.0) {
		msg_Debugging()<<"Veto kinematics: "<<cf[f]<<" = "
			       <<ID(li->Id())<<" & "<<ID(lj->Id())
			       <<" <-> "<<ID(lk->Id())<<"\n";
		trials.insert(cs);
		continue;
	      }
	      msg_Debugging()<<ID(li->Id())<<" & "<<ID(lj->Id())
			     <<" <-> "<<ID(lk->Id())<<" ["<<cf[f]
			     <<"]: "<<sqrt(cs.m_kt2)
			     <<" ("<<sqrt(cs.m_op2)<<")\n";
	      if (cs.m_op2>ckw.m_op2) ckw=cs;
	    }
	  }
	}
      }
    }
    msg_Debugging()<<"}\n";
    trials.insert(ckw);
    if (ckw.m_i==0 && ckw.m_j==0 && ckw.m_k==0) {
      double kt2core(CoreScale(ampl).m_op2);
      bool ord(true);
      std::vector<std::pair<size_t,double> > 
    	&pops(ops[ampl->Legs().size()-(m_nin+m_noutmin)]);
      if ((m_order || (m_rproc && ampl->Prev() &&
		      ampl->Prev()->Prev()==NULL)) &&
	  kt2core<pops.front().second) {
    	msg_Debugging()<<"unordered configuration (core): "
    		       <<sqrt(kt2core)<<" vs. "
    		       <<sqrt(pops.front().second)<<" "
    		       <<ID(pops.front().first)<<"\n";
    	ampl=ampl->Prev();
    	ampl->DeleteNext();
    	continue;
      }
      else {
	msg_Debugging()<<"Final configuration:\n";
	msg_Debugging()<<*ampl<<"\n";
	while (ampl->Prev()) {
	  ampl=ampl->Prev();
	  msg_Debugging()<<*ampl<<"\n";
	}
	double muf2(SetScales(ampl,mode));
	return muf2;
      }
    }
    msg_Debugging()<<"Cluster "<<ckw.m_fl
		   <<" "<<ID(ampl->Leg(ckw.m_i)->Id())
		   <<" & "<<ID(ampl->Leg(ckw.m_j)->Id())
		   <<" <-> "<<ID(ampl->Leg(ckw.m_k)->Id())
		   <<" => "<<sqrt(ckw.m_kt2)<<"\n";
    std::vector<std::pair<size_t,double> > 
      &cops(ops[ampl->Legs().size()-(m_nin+m_noutmin+1)]),
      &pops(ops[ampl->Legs().size()-(m_nin+m_noutmin)]);
    cops=pops;
    size_t sid(ampl->Leg(ckw.m_i)->Id()|ampl->Leg(ckw.m_j)->Id());
    size_t lmin(100), li(0);
    for (size_t i(0);i<cops.size();++i) {
      if ((cops[i].first&sid)==sid &&
	  IdCount(cops[i].first)<lmin) {
	lmin=IdCount(cops[i].first);
	li=i;
      }
    }
    cops[li].second=ckw.m_kt2;
    if (!m_orderrs && m_rproc &&
	ampl->Prev()==NULL) cops[li].second=0.0;
    msg_Debugging()<<"set last k_T = "<<sqrt(cops[li].second)
		   <<" "<<ID(cops[li].first)<<"\n";
    if ((m_order || (m_rproc && ampl->Prev() &&
		     ampl->Prev()->Prev()==NULL)) &&
	cops[li].second<pops[li].second) {
      msg_Debugging()<<"unordered configuration: "
		     <<sqrt(cops[li].second)<<" vs. "
		     <<sqrt(pops[li].second)<<" "
		     <<ID(pops[li].first)<<"\n";
      continue;
    }
    ampl->SetKT2(ckw.m_kt2);
    ampl->SetMu2(ckw.m_kt2);
    ampl=ampl->InitNext();
    ampl->CopyFrom(ampl->Prev());
    if (!Combine(*ampl,ckw)) {
      msg_Debugging()<<"combine failed\n";
      ampl=ampl->Prev();
      ampl->DeleteNext();
      continue;
    }
    ampl->SetOrderQCD(ampl->OrderQCD()-1);
  }
  THROW(fatal_error,"Invalid amplitude");
  return 0.0;
}

PDF::CParam MINLO_Scale_Setter::CoreScale(Cluster_Amplitude *const ampl) const
{
  ampl->SetProc(p_proc);
  PDF::CParam kt2(p_core->Calculate(ampl));
  ampl->SetKT2(kt2.m_kt2);
  ampl->SetMu2(kt2.m_mu2);
  for (Cluster_Amplitude *campl(ampl);
       campl;campl=campl->Prev()) campl->SetMuQ2(kt2.m_op2);
  return kt2;
}

double MINLO_Scale_Setter::SetScales(Cluster_Amplitude *ampl,const size_t &mode)
{
  m_scale[stp::res]=ampl->MuQ2();
  m_scale[stp::fac]=m_q02[0]=m_q02[1]=ampl->KT2();
  std::vector<double> scale(p_proc->NOut()+1);
  msg_Debugging()<<"Setting scales {\n";
  double mur2(1.0), as(1.0), ass(0.0), oqcd(0.0);
  for (size_t idx(2);ampl->Next();++idx,ampl=ampl->Next()) {
    scale[idx]=Max(ampl->Mu2(),m_rsf*MODEL::as->CutQ2());
    scale[idx]=Min(scale[idx],sqr(rpa->gen.Ecms()));
    if (m_rproc && ampl->Prev()==NULL) {
      m_scale[stp::fac]=m_q02[0]=m_q02[1]=ampl->Next()->KT2();
      continue;
    }
    if (ampl->KT2()>m_scale[stp::res]) m_scale[stp::res]=ampl->KT2();
    double coqcd(ampl->OrderQCD()-ampl->Next()->OrderQCD());
    if (coqcd!=1.0) THROW(fatal_error,"Non-QCD clustering");
    double cas(MODEL::as->BoundedAlphaS(m_rsf*scale[idx]));
    msg_Debugging()<<"  \\mu_{"<<idx<<"} = "
		   <<sqrt(m_rsf)<<" * "<<sqrt(scale[idx])
		   <<", as = "<<cas<<", O(QCD) = "<<coqcd<<"\n";
    mur2*=pow(m_rsf*scale[idx],coqcd);
    as*=pow(cas,coqcd);
    ass+=cas*coqcd;
    oqcd+=coqcd;
  }
  if (ampl->OrderQCD()-(m_vproc?1:0)) {
    double mu2(Max(ampl->Mu2(),m_rsf*MODEL::as->CutQ2()));
    int coqcd(ampl->OrderQCD()-(m_vproc?1:0));
    double cas(MODEL::as->BoundedAlphaS(m_rsf*mu2));
    msg_Debugging()<<"  \\mu_{0} = "<<sqrt(m_rsf)<<" * "<<sqrt(mu2)
		   <<", as = "<<cas<<", O(QCD) = "<<coqcd<<"\n";
    mur2*=pow(m_rsf*mu2,coqcd);
    as*=pow(cas,coqcd);
    ass+=cas*coqcd;
    oqcd+=coqcd;
  }
  if (oqcd==0.0) mur2=m_rsf*ampl->Mu2();
  else {
    mur2=pow(mur2,1.0/oqcd);
    ass/=oqcd;
    if (m_nproc && m_nlocpl==1) {
      msg_Debugging()<<"  as_{NLO} = "<<ass<<"\n";
      as=pow(as*ass,1.0/(oqcd+1.0));
    }
    else {
      as=pow(as,1.0/oqcd);
    }
    mur2=MODEL::as->WDBSolve(as,m_rsf*MODEL::as->CutQ2(),
			     m_rsf*1.01*sqr(rpa->gen.Ecms()));
    if (!IsEqual((*MODEL::as)(mur2),as))
      msg_Error()<<METHOD<<"(): Failed to determine \\mu."<<std::endl; 
  }
  msg_Debugging()<<"} -> as = "<<as<<" -> "<<sqrt(mur2)<<"\n";
  m_scale[stp::ren]=mur2;
  msg_Debugging()<<"Core / QCD scale = "<<sqrt(m_scale[stp::fac])
		 <<" / "<<sqrt(m_scale[stp::ren])<<"\n";
  for (size_t i(0);i<m_calcs.size();++i)
    m_scale[i]=m_calcs[i]->Calculate()->Get<double>();
  m_scale[stp::fac]=Max(m_scale[stp::fac],m_muf2min);
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(m_scale[stp::fac])<<"\n"
		 <<"  \\mu_r = "<<sqrt(m_scale[stp::ren])<<"\n"
		 <<"  \\mu_q = "<<sqrt(m_scale[stp::res])<<"\n";
  for (size_t i(stp::size);i<m_scale.size();++i)
    msg_Debugging()<<"  \\mu_"<<i<<" = "<<sqrt(m_scale[i])<<"\n";
  msg_Debugging()<<"} <- "<<(p_caller?p_caller->Name():"")<<"\n";
  ampl->SetMuF2(m_scale[stp::fac]);
  ampl->SetMuR2(m_scale[stp::ren]);
  ampl->SetMuQ2(m_scale[stp::res]);
  while (ampl->Prev()) {
    ampl=ampl->Prev();
    ampl->SetMuF2(m_scale[stp::fac]);
    ampl->SetMuR2(m_scale[stp::ren]);
    ampl->SetMuQ2(m_scale[stp::res]);
  }
  if (m_mufmode==1) m_q02[1]=m_scale[stp::fac];
  return m_scale[stp::fac];
}

void MINLO_Scale_Setter::SetScale
(const std::string &mu2tag,Algebra_Interpreter &mu2calc)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<p_caller->Name()<<"' {\n";
  msg_Indent();
  m_tagset.SetTags(&mu2calc);
  mu2calc.Interprete(mu2tag);
  if (msg_LevelIsDebugging()) mu2calc.PrintEquation();
  msg_Debugging()<<"}\n";
}
  
bool MINLO_Scale_Setter::Combine
(Cluster_Amplitude &ampl,const MCS_Params &cs) const
{
  int i(cs.m_i), j(cs.m_j), k(cs.m_k);
  if (i>j) std::swap<int>(i,j);
  Cluster_Leg *li(ampl.Leg(i)), *lj(ampl.Leg(j)), *lk(ampl.Leg(k));
  li->SetFlav(cs.m_fl);
  li->SetMom(cs.m_pijt);
  lk->SetMom(cs.m_pkt);
  for (size_t m(0);m<ampl.Legs().size();++m) ampl.Leg(m)->SetK(0);
  if (i<2) {
    for (size_t m(0);m<ampl.Legs().size();++m) {
      if ((int)m==i || (int)m==j || (int)m==k) continue;
      ampl.Leg(m)->SetMom(cs.m_lam*ampl.Leg(m)->Mom());
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

void MINLO_Scale_Setter::KT2
(const Cluster_Leg *li,const Cluster_Leg *lj,
 const Cluster_Leg *lk,MCS_Params &cs) const
{
  if ((li->Id()&3)<(lj->Id()&3)) std::swap<const Cluster_Leg*>(li,lj);
  Vec4D pi(li->Mom()), pj(lj->Mom()), pk(lk->Mom());
  double mi2=sqr(li->Flav().Mass()), mj2=sqr(lj->Flav().Mass());
  double mij2=sqr(cs.m_fl.Mass()), mk2=sqr(lk->Flav().Mass());
  if (li->Stat()==3) mi2=pi.Abs2();
  if (lj->Stat()==3) mj2=pj.Abs2();
  if (lk->Stat()==3) mk2=pk.Abs2();
  if (m_cmode==0) {// CSS algorithm
    if ((li->Id()&3)==0) {
      if ((lk->Id()&3)==0) {
	Kin_Args ffp(ClusterFFDipole(mi2,mj2,mij2,mk2,pi,pj,pk,3));
	double kt2=cs.m_kt2=2.0*(pi*pj)*ffp.m_z*(1.0-ffp.m_z)
	  -sqr(1.0-ffp.m_z)*mi2-sqr(ffp.m_z)*mj2;
	if (ffp.m_stat<0) kt2=-1.0;
	cs.SetParams(kt2,ffp.m_pi,ffp.m_pk,ffp.m_lam);
      }
      else {
	Kin_Args fip(ClusterFIDipole(mi2,mj2,mij2,mk2,pi,pj,-pk,3));
	double kt2=cs.m_kt2=2.0*(pi*pj)*fip.m_z*(1.0-fip.m_z)
	  -sqr(1.0-fip.m_z)*mi2-sqr(fip.m_z)*mj2;
	Vec4D sum(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
	if ((cs.m_k==0 && fip.m_pk[3]<0.0) ||
	    (cs.m_k==1 && fip.m_pk[3]>0.0) ||
	    fip.m_pk[0]<0.0 || fip.m_stat<0 ||
	    fip.m_pk[0]>rpa->gen.PBeam(cs.m_k)[0]) kt2=-1.0;
	cs.SetParams(kt2,fip.m_pi,-fip.m_pk,fip.m_lam);
      }
    }
    else {
      if ((lk->Id()&3)==0) {
	Kin_Args ifp(ClusterIFDipole(mi2,mj2,mij2,mk2,0.0,-pi,pj,pk,pk,3|4));
	double kt2=cs.m_kt2=-2.0*(pi*pj)*(1.0-ifp.m_y)*(1.0-ifp.m_z);
	Vec4D sum(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
	if ((cs.m_i==0 && ifp.m_pi[3]<0.0) ||
	    (cs.m_i==1 && ifp.m_pi[3]>0.0) ||
	    ifp.m_pi[0]<0.0 || ifp.m_stat<0 ||
	    ifp.m_pi[0]>rpa->gen.PBeam(cs.m_i)[0]) kt2=-1.0;
	cs.SetParams(kt2,-ifp.m_pi,ifp.m_pk,ifp.m_lam);
      }
      else {
	Kin_Args iip(ClusterIIDipole(mi2,mj2,mij2,mk2,-pi,pj,-pk,3));
	double kt2=cs.m_kt2=-2.0*(pi*pj)*(1.0-iip.m_z-iip.m_y);
	Vec4D sum(rpa->gen.PBeam(0)+rpa->gen.PBeam(1));
	if ((cs.m_i==0 && iip.m_pi[3]<0.0) ||
	    (cs.m_i==1 && iip.m_pi[3]>0.0) ||
	    iip.m_pi[0]<0.0 || iip.m_stat<0 ||
	    iip.m_pi[0]>rpa->gen.PBeam(cs.m_i)[0]) kt2=-1.0;
	cs.SetParams(kt2,-iip.m_pi,-iip.m_pk,iip.m_lam);
	double pp(pj.PPlus()), pm(pj.PMinus());
	if (pi[3]>0) std::swap<double>(pp,pm);
	if (pp>pm) cs.m_op2*=1.000001;
	else cs.m_op2/=1.000001;
      }
    }
  }
  else if (m_cmode==1) {// KT algorithm, E-scheme
    if ((li->Id()&3)==0) {
      if (li->Flav().Mass() || lj->Flav().Mass()) {
	// HQ: massive Durham algorithm
	double kt2=2.0*Min(pi.PSpat2(),pj.PSpat2())*
	  (1.0-pi.CosTheta(pj));
	cs.SetParams(kt2,pi+pj,pk);
      }
      else {
	// LQ: kT algorithm
	double kt2=Min(pi.PPerp2(),pj.PPerp2())*
	  (cosh(pi.DY(pj))-cos(pi.DPhi(pj)))/sqr(m_dr);
	cs.SetParams(kt2,pi+pj,pk);
      }
    }
    else {
      double kt2=pj.PPerp2();
      cs.SetParams(kt2,pi+pj,pk);
      double pp(pj.PPlus()), pm(pj.PMinus());
      if (pi[3]>0) std::swap<double>(pp,pm);
      if (pp>pm) cs.m_op2*=1.000001;
      else cs.m_op2/=1.000001;
    }
  }
  else {
    THROW(fatal_error,"Invalid clustering scheme");
  }
}
