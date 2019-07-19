#include "DIM/Shower/Shower.H"

#include "DIM/Main/Gamma.H"
#include "DIM/Tools/Amplitude.H"
#include "DIM/Tools/Weight.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PDF/Main/Jet_Criterion.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Variations.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace DIM;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

Shower::Shower():
  p_model(NULL), p_as(NULL),
  p_gamma(NULL)
{
  p_pdf[1]=p_pdf[0]=NULL;

  Settings& s = Settings::GetMainSettings();
  std::string dipole_string = s["DIPOLES"]["CASE"].Get<std::string>();
  if      (dipole_string == "CS")    m_dipole_case = EXTAMP::DipoleCase::CS;
  else if (dipole_string == "IDa")   m_dipole_case = EXTAMP::DipoleCase::IDa;  // ee > bbWW with mapping a
  else if (dipole_string == "IDb")   m_dipole_case = EXTAMP::DipoleCase::IDb;  // ee > bbWW with mapping b
  else if (dipole_string == "IDin")  m_dipole_case = EXTAMP::DipoleCase::IDin; // pp > bbWW
  else if (dipole_string == "RES")   m_dipole_case = EXTAMP::DipoleCase::RES;  // ee > guu
  else if (dipole_string == "ID")    m_dipole_case = EXTAMP::DipoleCase::ID;   // ee > guu
  else                               m_dipole_case = EXTAMP::DipoleCase::CS;
}

Shower::~Shower()
{
  for (Kernel_Vector::const_iterator it(m_cks.begin());
       it!=m_cks.end();++it) delete *it;
}

struct FTrip {
public:
  Flavour m_a, m_b, m_c;
public:
  inline FTrip(const Flavour &a,const Flavour &b,const Flavour &c):
    m_a(a), m_b(b), m_c(c) {}
  bool operator<(const FTrip &f) const
  {
    if (m_a<f.m_a) return true;
    if (m_a==f.m_a) {
      if (m_b<f.m_b) return true;
      if (m_b==f.m_b) {
	return m_c<f.m_c;
      }
    }
    return false;
  }
};

void Shower::Init(MODEL::Model_Base *const model,
		  PDF::ISR_Handler *const isr)
{
  DEBUG_FUNC(this);
  Settings& s = Settings::GetMainSettings();
  p_model=model;
  p_as=(MODEL::Running_AlphaS*)p_model->GetScalarFunction("alpha_S");
  for (int i=0;i<2;++i) p_pdf[i]=isr->PDF(i);
  m_tmin[0] = s["CSS_FS_PT2MIN"].Get<double>();
  m_tmin[1] = s["CSS_IS_PT2MIN"].Get<double>();
  m_cplfac[0] = s["CSS_FS_AS_FAC"].Get<double>();
  m_cplfac[1] = s["CSS_IS_AS_FAC"].Get<double>();
  m_rsf=ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
  m_fsf=ToType<double>(rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR"));
  m_rcf = s["NLO_CSS_RECALC_FACTOR"].Get<double>();
  m_kfac = s["CSS_KFACTOR_SCHEME"].Get<int>();
  m_cpl = s["CSS_COUPLING_SCHEME"].Get<int>();
  m_pdfmin[0]=s["CSS_PDF_MIN"].Get<double>();
  m_pdfmin[1]=s["CSS_PDF_MIN_X"].Get<double>();
  m_maxem=s["NLO_CSS_MAXEM"].Get<unsigned int>();
  m_maxrewem=s["REWEIGHT_MCATNLO_EM"].Get<unsigned int>();
  m_rewtmin=s["CSS_REWEIGHT_SCALE_CUTOFF"].Get<double>();
  m_oef=s["CSS_OEF"].Get<double>();
  if (msg_LevelIsDebugging()) {
    msg_Out()<<METHOD<<"(): {\n\n"
	     <<"   // available gauge calculators\n\n";
    Gauge_Getter::PrintGetterInfo(msg->Out(),25);
    msg_Out()<<"\n   // available lorentz calculators\n\n";
    Lorentz_Getter::PrintGetterInfo(msg->Out(),25);
    msg_Out()<<"\n}"<<std::endl;
  }
  int types(s["CSS_KERNEL_TYPE"].Get<int>());
  std::set<FTrip> sfs;
  const Vertex_Table *vtab(model->VertexTable());
  for (Vertex_Table::const_iterator
	 vlit=vtab->begin();vlit!=vtab->end();++vlit) {
    for (Vertex_List::const_iterator 
	   vit=vlit->second.begin();vit!=vlit->second.end();++vit) {
      Single_Vertex *v(*vit);
      if (v->NLegs()>3) continue;
      if (sfs.find(FTrip(v->in[0],v->in[1],v->in[2]))
	  !=sfs.end()) continue;
      msg_Indent();
      sfs.insert(FTrip(v->in[0],v->in[1],v->in[2]));
      sfs.insert(FTrip(v->in[0],v->in[2],v->in[1]));
      msg_IODebugging()<<"Add "<<v->in[0].Bar()<<" -> "
		       <<v->in[1]<<" "<<v->in[2]<<" {\n";
      {
	msg_Indent();
	for (int type(0);type<4;++type)
	  if (types&(1<<type))
	    for (int mode(0);mode<2;++mode)
	      AddKernel(new Kernel(this,Kernel_Key(v,mode,type)));
      }
      msg_IODebugging()<<"}\n";
    }
  }
}

void Shower::AddKernel(Kernel *const k)
{
  if (k->On()<0) {
    delete k;
    return;
  }
  k->GF()->SetLimits();
  if (k->On()) m_sks[k->LF()->Flav(0)].push_back(k);
  m_cks.push_back(k);
  m_kmap[k->Type()|(k->Type()&1?(k->Mode()?4:0):0)]
    [k->LF()->Flav(1)][k->LF()->Flav(2)]=k;
}


void Shower::SetMS(ATOOLS::Mass_Selector *const ms)
{
  for (Kernel_Vector::const_iterator
	 it(m_cks.begin());it!=m_cks.end();++it)
    (*it)->LF()->SetMS(ms);
}

void Shower::AddWeight(const Amplitude &a,const double &t)
{
  double cw(1.0);
  std::vector<double> cv;
  for (size_t i(0);i<a.size();++i) {
    cw*=a[i]->GetWeight(Max(t,m_tmin[a[i]->Beam()?1:0]),cv);
    a[i]->ClearWeights();
  }
  m_weight*=cw;
  if (cv.size()) p_vars->UpdateOrInitialiseWeights
		   (&Shower::GetWeight,*this,cv,Variations_Type::sudakov);
  msg_Debugging()<<a<<" t = "<<t<<" -> w = "<<cw
		 <<" ("<<m_weight<<"), v = "<<cv<<"\n";
}

double Shower::GetWeight(Variation_Parameters *params,
			 Variation_Weights *weights,
			 std::vector<double> &v)
{
  return v[weights->CurrentParametersIndex()];
}

double Shower::VetoWeight(Variation_Parameters *params,
			  Variation_Weights *weights,
			  JetVeto_Args &args)
{
  msg_Debugging()<<METHOD<<"("<<weights<<"){\n";
  int stat(args.m_jcv>=0.0);
  Cluster_Amplitude *ampl(args.p_ampl);
  Jet_Finder *jf(ampl->JF<Jet_Finder>());
  if (stat && jf) {
    double fac(weights?params->m_Qcutfac:1.0);
    if (jf) {
      stat=args.m_jcv<sqr(jf->Qcut()*fac);
      msg_Debugging()<<"  jcv = "<<sqrt(args.m_jcv)<<" vs "
		     <<jf->Qcut()<<" * "<<fac
		     <<" = "<<jf->Qcut()*fac<<"\n";
    }
  }
  if (stat==1) {
    msg_Debugging()<<"} no jet veto\n";
    args.m_acc=1;
    return 1.0;
  }
  msg_Debugging()<<"} jet veto\n";
  return 0.0;
}

int Shower::Evolve(Amplitude &a,unsigned int &nem)
{
  DEBUG_FUNC(this);
  m_weight=1.0;
  msg_Debugging()<<a<<"\n";
  if (nem>=m_maxem) return 1;
  double t(a.T());
  Cluster_Amplitude *ampl(a.ClusterAmplitude());
  for (m_s=Splitting(GeneratePoint(a,t,nem));
       m_s.m_t>Max(a.T0(),m_tmin[m_s.m_type&1]);
       m_s=GeneratePoint(a,m_s.m_t,nem)) {
    for (size_t i(0);i<a.size();++i) a[i]->Store();
    int stat(m_s.p_sk->Construct(m_s,0));
    msg_IODebugging()<<"t = "<<m_s.m_t<<", w = "<<m_s.m_w.MC()
		     <<" / "<<m_s.m_w.Accept()<<" -> "
		     <<(stat==1?"accept\n":"reject\n");
    msg_Debugging()<<"stat = "<<stat<<"\n";
    if (p_gamma) {
      int veto(p_gamma->Reject());
      m_weight*=p_gamma->Weight();
      if (p_vars) {
	std::vector<double> gwa
	  (p_vars->NumberOfParameters(),p_gamma->Weight());
	p_vars->UpdateOrInitialiseWeights
	  (&Shower::GetWeight,*this,gwa,Variations_Type::sudakov);
      }
      if (veto) {
	a.Remove(m_s.p_n);
	m_s.p_c->SetFlav(m_s.p_sk->LF()->Flav(0));
	for (size_t i(0);i<a.size();++i) a[i]->Restore();
	continue;
      }
    }
    JetVeto_Args vwa(ampl,stat?0.0:-1.0);
    Jet_Finder *jf=m_s.p_c->Ampl()->JF<Jet_Finder>();
    if (stat && jf) {
      Cluster_Amplitude *ampl(a.GetAmplitude());
      vwa.m_jcv=jf->JC()->Value(ampl);
      ampl->Delete();
    }
    m_weight*=VetoWeight(NULL,NULL,vwa);
    if (p_vars) p_vars->UpdateOrInitialiseWeights
		  (&Shower::VetoWeight,*this,vwa,Variations_Type::sudakov);
    if (vwa.m_acc==0) return 0;
    AddWeight(a,m_s.m_t);
    a.SetJF(NULL);
    if (++nem>=m_maxem) break;
  }
  AddWeight(a,a.T0());
  return 1;
}

Splitting Shower::GeneratePoint
(const Amplitude &a,const double &t,const unsigned int &nem)
{
  Splitting win;
  double tmin[2]={m_tmin[0],m_tmin[1]};
  for (Amplitude::const_reverse_iterator
	 it(a.rbegin());it!=a.rend();++it) {
    Splitting cur(GeneratePoint(**it,t,nem));
    if (cur.p_c==NULL || cur.p_s==NULL) continue;
    if (cur.m_t<m_tmin[cur.m_type&1]) continue;
    m_tmin[0]=m_tmin[1]=(win=cur).m_t;
  }
  m_tmin[0]=tmin[0];
  m_tmin[1]=tmin[1];
  if (win.p_sk && win.m_t>m_tmin[win.m_type&1])
    msg_Debugging()<<"Emission at "<<win<<"\n";
  return win;
}

Splitting Shower::GeneratePoint
(Parton &p,const double &t,const unsigned int &nem)
{
  Splitting win(&p,NULL,t);
  double sum=0.0, ct=m_rcf*t;
  SKernel_Map::const_iterator kit(m_sks.find(p.Flav()));
  if (kit==m_sks.end()) return Splitting();
  std::vector<std::vector<double> > psum(kit->second.size());
  std::vector<std::vector<size_t> > splits(psum.size());
  while (true) {
    if (win.m_t*m_rcf<=ct) {
      sum=0.0;
      ct=win.m_t;
      for (size_t j(0);j<kit->second.size();++j) {
	psum[j].clear();
	splits[j].clear();
	double csum=0.0;
	for (int i(p.Ampl()->size()-1);i>=0;--i) {
	  if ((*p.Ampl())[i]==&p) continue;
	  Splitting cur(&p,(*p.Ampl())[i]);
	  cur.SetType();
      cur.SetKinSpect(p);
      if (cur.p_kinspec!=NULL) cur.SetKinVar();
	  cur.m_kfac=m_kfac;
	  cur.m_cpl=m_cpl;
	  cur.m_t1=ct;
	  if (kit->second[j]->Allowed(cur)) {
	    double I=kit->second[j]->Integral(cur);
	    psum[j].push_back(csum+=dabs(I));
	    splits[j].push_back(i);
	  }
	}
	if (psum[j].size()) sum+=psum[j].back();
      }
      if (sum==0.0) return Splitting();
      win=Splitting(&p,NULL,ct);
    }
    win.m_t*=exp(log(ran->Get())*Max(2.0*M_PI/sum,1.0e-3));
    if (win.m_t<m_tmin[p.Beam()?1:0]) return win;
    double disc(sum*ran->Get()), csum(0.0);
    for (size_t j(0);j<splits.size();++j)
      if (splits[j].size() &&
	  (csum+=psum[j].back())>=disc) {
	double disc(psum[j].back()*ran->Get());
	for (size_t i(0);i<splits[j].size();++i)
	  if (psum[j][i]>=disc) {
	    win.p_s=(*p.Ampl())[splits[j][i]];
        win.SetKinSpect(p);
        if (win.p_kinspec==NULL) THROW(fatal_error, "Invalid spectator assignment!");
        win.SetKinVar();
	    win.SetType();
	    win.m_kfac=m_kfac;
	    win.m_cpl=m_cpl;
	    win.m_t1=ct;
	    if (!kit->second[j]->GeneratePoint(win)) {
	      msg_Error()<<METHOD<<"(): Error generating point!\n";
	      msg_Debugging()<<win<<"\nQ2 = "<<win.m_Q2
			     <<", eta = "<<win.m_eta
			     <<", t0 = "<<win.m_t0<<"\n";
	      break;
	    }
	    if (!kit->second[j]->LF()->Compute(win)) break;
	    win.m_w=kit->second[j]->GetWeight(win,m_oef);
	    win.m_vars=std::vector<double>(p_vars->NumberOfParameters(),1.0);
	    if (win.m_w.MC()<ran->Get()) {
	      if (p_vars && nem<m_maxrewem && win.m_t>m_rewtmin) {
		const Reweight_Args args(&win,0);
		p_vars->UpdateOrInitialiseWeights
		  (&Shower::Reweight,*this,args,Variations_Type::sudakov);
	      }
	      win.p_c->AddWeight(win,0);
	      msg_IODebugging()<<"t = "<<win.m_t<<", w = "<<win.m_w.MC()
			       <<" / "<<win.m_w.Reject()<<" -> reject ["
			       <<win.p_c->Id()<<"<->"<<win.p_s->Id()<<"]\n";
	      break;
	    }
	    if (p_vars && nem<m_maxrewem && win.m_t>m_rewtmin) {
	      const Reweight_Args args(&win,1);
	      p_vars->UpdateOrInitialiseWeights
		(&Shower::Reweight,*this,args,Variations_Type::sudakov);
	    }
	    win.p_c->AddWeight(win,1);
	    msg_IODebugging()<<"t = "<<win.m_t<<", w = "<<win.m_w.MC()
			     <<" / "<<win.m_w.Accept()<<" -> select ["
			     <<win.p_c->Id()<<"<->"<<win.p_s->Id()<<"]\n";
	    return win;
	  }
	break;
      }
  }
  return win;
}

double Shower::Reweight(Variation_Parameters *params,
			Variation_Weights *weights,
			const Reweight_Args &a)
{
  double rsf(m_rsf), fsf(m_fsf);
  m_rsf*=params->m_muR2fac;
  m_fsf*=params->m_muF2fac;
  MODEL::Running_AlphaS *as(p_as);
  p_as=params->p_alphas;
  PDF::PDF_Base *pdf[2]={p_pdf[0],p_pdf[1]};
  p_pdf[0]=params->p_pdf1;
  p_pdf[1]=params->p_pdf2;
  size_t cur(weights->CurrentParametersIndex());
  if (rsf==m_rsf && fsf==m_fsf && as==p_as &&
      pdf[0]==p_pdf[0] && pdf[1]==p_pdf[1]) {
    a.m_s->m_vars[cur]=1.0;
    return 1.0;
  }
  msg_IODebugging()<<METHOD<<"("<<cur<<") {\n  "
		   <<"\\mu_R -> "<<sqrt(m_rsf)
		   <<", \\mu_F -> "<<sqrt(m_fsf)<<"\n  PDF "
		   <<(p_pdf[0]?p_pdf[0]->LHEFNumber():-1)<<" x "
		   <<(p_pdf[1]?p_pdf[1]->LHEFNumber():-1)<<"\n";
  MC_Weight w(a.m_s->p_sk->GetWeight(*a.m_s,m_oef,&a.m_s->m_w));
  msg_IODebugging()<<"  w_ref = "<<a.m_s->m_w
		   <<"\n  w_new = "<<w<<"\n";
  if (a.m_acc) a.m_s->m_vars[cur]=w.MC()/a.m_s->m_w.MC();
  else a.m_s->m_vars[cur]=w.Reject()/a.m_s->m_w.Reject()
    	 *(1.0-w.MC())/(1.0-a.m_s->m_w.MC());
  msg_IODebugging()<<"} -> w = "<<a.m_s->m_vars[cur]<<"\n";
  p_pdf[0]=pdf[0];
  p_pdf[1]=pdf[1];
  p_as=as;
  m_rsf=rsf;
  m_fsf=fsf;
  return 1.0;
}

double Shower::GetXPDF
(const double &x,const double &Q2,
 const ATOOLS::Flavour &fl,const int b) const
{
  if (p_pdf[b]==NULL) return 1.0;
  if (!p_pdf[b]->Contains(fl.Bar())) {
    if (fl.Strong() || fl.Mass()<10.0) return 0.0;
    return 1.0;
  }
  if (Q2<sqr(2.0*fl.Mass(true))) return 0.0;
  if (x<p_pdf[b]->XMin() || x>p_pdf[b]->XMax() ||
      Q2<p_pdf[b]->Q2Min() || Q2>p_pdf[b]->Q2Max())
    return 0.0;
  p_pdf[b]->Calculate(x,m_fsf*Q2);
  return p_pdf[b]->GetXPDF(fl.Bar());
}

Kernel *Shower::GetKernel(const Splitting &s,const int mode) const
{
  Kernel_Map::const_iterator seit(m_kmap.find(s.m_type|(mode?4:0)));
  if (seit==m_kmap.end()) return NULL;
  SEKernel_Map::const_iterator eit(seit->second.find(s.p_c->Flav()));
  if (eit==seit->second.end()) return NULL;
  EKernel_Map::const_iterator it(eit->second.find(s.p_n->Flav()));
  if (it==eit->second.end()) return NULL;
  if (s.p_s==NULL || it->second->GF()->Allowed(s)) return it->second;
  return NULL;
}
