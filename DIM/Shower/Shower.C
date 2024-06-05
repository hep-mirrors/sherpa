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
  auto pss = s["SHOWER"], nlopss = s["MC@NLO"];
  p_model=model;
  p_as=(MODEL::Running_AlphaS*)p_model->GetScalarFunction("alpha_S");
  for (int i=0;i<2;++i) p_pdf[i]=isr->PDF(i);
  m_tmin[0] = pss["FS_PT2MIN"].Get<double>();
  m_tmin[1] = pss["IS_PT2MIN"].Get<double>();
  m_cplfac[0] = pss["FS_AS_FAC"].Get<double>();
  m_cplfac[1] = pss["IS_AS_FAC"].Get<double>();
  m_rsf=ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
  m_fsf=ToType<double>(rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR"));
  m_rcf = nlopss["RECALC_FACTOR"].Get<double>();
  m_kfac = pss["KFACTOR_SCHEME"].Get<int>();
  m_cpl = pss["COUPLING_SCHEME"].Get<int>();
  m_pdfmin[0]=pss["PDF_MIN"].Get<double>();
  m_pdfmin[1]=pss["PDF_MIN_X"].Get<double>();
  m_maxem=nlopss["MAXEM"].Get<unsigned int>();
  m_maxrewem=nlopss["REWEIGHT_EM"].Get<unsigned int>();
  m_rewtmin=pss["REWEIGHT_SCALE_CUTOFF"].Get<double>();
  m_oef=pss["OEF"].Get<double>();
  if (msg_LevelIsDebugging()) {
    msg_Out()<<METHOD<<"(): {\n\n"
	     <<"   // available gauge calculators\n\n";
    Gauge_Getter::PrintGetterInfo(msg->Out(),25);
    msg_Out()<<"\n   // available lorentz calculators\n\n";
    Lorentz_Getter::PrintGetterInfo(msg->Out(),25);
    msg_Out()<<"\n}"<<std::endl;
  }
  int types(pss["KERNEL_TYPE"].Get<int>());
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
  m_weightsmap["MC@NLO_PS"].Nominal() *= cw;
  if (cv.size()) {
    ATOOLS::Reweight(m_weightsmap["MC@NLO_PS"],
                     [&cv](double varweight,
                           size_t varindex,
                           QCD_Variation_Params& varparams) -> double {
                       return varweight * cv[varindex];
                     });
  }
  msg_Debugging()<<a<<" t = "<<t<<" -> w = "<<cw
		 <<" ("<<m_weightsmap["MC@NLO_PS"].Nominal()<<"), v = "<<cv<<"\n";
}

int Shower::Evolve(Amplitude &a,unsigned int &nem)
{
  DEBUG_FUNC(this);
  m_weightsmap.Clear();
  m_weightsmap["MC@NLO_PS"] = Weights {Variations_Type::qcd};
  m_weightsmap["MC@NLO_QCUT"] = Weights {Variations_Type::qcut};
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
      m_weightsmap["MC@NLO_PS"] *= p_gamma->Weight();
      if (veto) {
	a.Remove(m_s.p_n);
	m_s.p_c->SetFlav(m_s.p_sk->LF()->Flav(0));
	for (size_t i(0);i<a.size();++i) a[i]->Restore();
	continue;
      }
    }
    double jcv {stat ? 0.0 : -1.0};
    Jet_Finder *jf=m_s.p_c->Ampl()->JF<Jet_Finder>();
    if (stat && jf) {
      Cluster_Amplitude *ampl(a.GetAmplitude());
      jcv = jf->JC()->Value(ampl);
      ampl->Delete();
    }
    const bool is_jcv_positive {jcv >= 0.0};
    bool all_vetoed {true};
    ATOOLS::ReweightAll(
        m_weightsmap["MC@NLO_QCUT"],
        [this, jcv, is_jcv_positive, ampl, &all_vetoed](
            double varweight,
            size_t varindex,
            Qcut_Variation_Params* qcutparams) -> double {
          msg_Debugging() << "Applying veto weight to " << varweight << " {\n";
          bool stat {is_jcv_positive};
          Jet_Finder* jf(ampl->JF<Jet_Finder>());
          if (stat && jf) {
            const double fac {
                qcutparams == nullptr ? 1.0 : qcutparams->m_scale_factor};
            stat = jcv < sqr(jf->Qcut() * fac);
            msg_Debugging()
                << "  jcv = " << sqrt(jcv) << " vs " << jf->Qcut() << " * "
                << fac << " = " << jf->Qcut() * fac << "\n";
          }
          if (stat) {
            msg_Debugging() << "} no jet veto\n";
            all_vetoed = false;
            return varweight;
          } else {
            msg_Debugging() << "} jet veto\n";
            return 0.0;
          }
        });
    if (all_vetoed)
      return 0;
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
            win.m_vars = std::vector<double>(s_variations->Size(), 1.0);
            if (win.m_w.MC() < ran->Get()) {
              if (nem < m_maxrewem && win.m_t > m_rewtmin) {
                const Reweight_Args args(&win, 0);
                s_variations->ForEach(
                    [this, &args](size_t varindex,
                                  QCD_Variation_Params& varparams) -> void {
                      Reweight(&varparams, varindex, args);
                    });
              }
              win.p_c->AddWeight(win, 0);
              msg_IODebugging()<<"t = "<<win.m_t<<", w = "<<win.m_w.MC()
			       <<" / "<<win.m_w.Reject()<<" -> reject ["
			       <<win.p_c->Id()<<"<->"<<win.p_s->Id()<<"]\n";
	      break;
            }
            if (nem < m_maxrewem && win.m_t > m_rewtmin) {
              const Reweight_Args args(&win, 1);
              s_variations->ForEach(
                  [this, &args](size_t varindex,
                                QCD_Variation_Params& varparams) -> void {
                    Reweight(&varparams, varindex, args);
                  });
            }
            win.p_c->AddWeight(win, 1);
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

void Shower::Reweight(QCD_Variation_Params* params,
                      size_t varindex,
                      const Reweight_Args& a)
{
  double rsf(m_rsf), fsf(m_fsf);
  m_rsf*=params->m_muR2fac;
  m_fsf*=params->m_muF2fac;
  MODEL::Running_AlphaS *as(p_as);
  p_as=params->p_alphas;
  PDF::PDF_Base *pdf[2]={p_pdf[0],p_pdf[1]};
  p_pdf[0]=params->p_pdf1;
  p_pdf[1]=params->p_pdf2;
  if (rsf==m_rsf && fsf==m_fsf && as==p_as &&
      pdf[0]==p_pdf[0] && pdf[1]==p_pdf[1]) {
    a.m_s->m_vars[varindex]=1.0;
    return;
  }
  msg_IODebugging()<<METHOD<<"("<<varindex<<") {\n  "
		   <<"\\mu_R -> "<<sqrt(m_rsf)
		   <<", \\mu_F -> "<<sqrt(m_fsf)<<"\n  PDF "
		   <<(p_pdf[0]?p_pdf[0]->LHEFNumber():-1)<<" x "
		   <<(p_pdf[1]?p_pdf[1]->LHEFNumber():-1)<<"\n";
  MC_Weight w(a.m_s->p_sk->GetWeight(*a.m_s,m_oef,&a.m_s->m_w));
  msg_IODebugging()<<"  w_ref = "<<a.m_s->m_w
		   <<"\n  w_new = "<<w<<"\n";
  if (a.m_acc)
    a.m_s->m_vars[varindex] = w.MC() / a.m_s->m_w.MC();
  else
    a.m_s->m_vars[varindex] = w.Reject() / a.m_s->m_w.Reject() *
                              (1.0 - w.MC()) / (1.0 - a.m_s->m_w.MC());
  msg_IODebugging()<<"} -> w = "<<a.m_s->m_vars[varindex]<<"\n";
  p_pdf[0]=pdf[0];
  p_pdf[1]=pdf[1];
  p_as=as;
  m_rsf=rsf;
  m_fsf=fsf;
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
