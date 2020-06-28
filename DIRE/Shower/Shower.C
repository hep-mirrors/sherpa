#include "DIRE/Shower/Shower.H"

#include "DIRE/Tools/Amplitude.H"
#include "DIRE/Shower/Cluster_Definitions.H"
#include "DIRE/Tools/Weight.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PDF/Main/Jet_Criterion.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Variations.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace DIRE;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

Shower::Shower():
  p_model(NULL),
  p_cluster(new Cluster_Definitions(this)),
  p_as(NULL)
{
  p_pdf[1]=p_pdf[0]=NULL;
}

Shower::~Shower()
{
  for (Kernel_Vector::const_iterator it(m_cks.begin());
       it!=m_cks.end();++it) delete *it;
  delete p_cluster;
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

bool Shower::Init(MODEL::Model_Base *const model,
		  PDF::ISR_Handler *const isr)
{
  Settings& s = Settings::GetMainSettings();
  DEBUG_FUNC(this);
  p_model=model;
  p_as=(MODEL::Running_AlphaS*)p_model->GetScalarFunction("alpha_S");
  for (int i=0;i<2;++i) p_pdf[i]=isr->PDF(i);
  m_tmin[0] = s["CSS_FS_PT2MIN"].Get<double>();
  m_tmin[1] = s["CSS_IS_PT2MIN"].Get<double>();
  m_cplfac[0] = s["CSS_FS_AS_FAC"].Get<double>();
  m_cplfac[1] = s["CSS_IS_AS_FAC"].Get<double>();
  m_rsf=ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
  m_fsf=ToType<double>(rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR"));
  m_rcf=s["CSS_RECALC_FACTOR"].Get<double>();
  m_tcef=s["CSS_TC_ENHANCE"].Get<double>();
  m_kin=s["CSS_KIN_SCHEME"].Get<int>();
  m_kfac=s["CSS_KFACTOR_SCHEME"].Get<int>();
  m_cpl=s["CSS_COUPLING_SCHEME"].Get<int>();
  m_mec=s["CSS_ME_CORRECTION"].Get<int>();
  m_pdfmin[0]=s["CSS_PDF_MIN"].Get<double>();
  m_pdfmin[1]=s["CSS_PDF_MIN_X"].Get<double>();
  m_maxem=s["CSS_MAXEM"].Get<size_t>();
  m_maxpart=s["CSS_MAXPART"].Get<int>();
  m_reweight=s["CSS_REWEIGHT"].Get<bool>();
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
      if (v->in[0].Kfcode()==6 || v->in[1].Kfcode()==6 || v->in[2].Kfcode()==6) {
	//msg_Out()<<"Do not include "<<v->in[0]<<" --> "<<v->in[1]<<"+"<<v->in[2]
	//	 <<" into shower.\n";
	continue;
      }
      if (sfs.find(FTrip(v->in[0],v->in[1],v->in[2]))
	  !=sfs.end()) continue;
      msg_Indent();
      sfs.insert(FTrip(v->in[0],v->in[1],v->in[2]));
      sfs.insert(FTrip(v->in[0],v->in[2],v->in[1]));
      msg_IODebugging()<<"Add "<<v->in[0].Bar()<<" -> "
		       <<v->in[1]<<" "<<v->in[2]<<" {\n";
      if (!(m_kfac&256)) {
	msg_Indent();
	for (int type(0);type<4;++type)
	  if (types&(1<<type))
	    for (int mode(0);mode<2;++mode)
	      for (int swap(0);swap<2;++swap)
		AddKernel(new Kernel(this,Kernel_Key(v,mode,swap,type)));
      }
      msg_IODebugging()<<"}\n";
    }
  }
  if (!((m_kfac&2) || (m_kfac&1024)) || (m_kfac&512)) return true;
  ATOOLS::Flavour_Vector fls(4);
  for (long int i(-5);i<=5;++i) {
    if (i==0) continue;
    fls[0]=(fls[2]=Flavour(i)).Bar();
    for (long int j(-5);j<=5;++j) {
      if (j==0 || j==i) continue;
      fls[3]=(fls[1]=Flavour(j)).Bar();
      for (int type(0);type<4;++type)
	if (types&(1<<type))
	  AddKernel(new Kernel(this,Kernel_Key(fls,1,type,"FFFF")));
    }
  }
  return true;
}

void Shower::AddKernel(Kernel *const k)
{
  if (k->On()<0) {
    delete k;
    return;
  }
  k->GF()->SetLimits();
  if (k->On()) m_sks[k->LF()->Flav(0)].push_back(k);
  if (k->LF()->Flavs().size()>3) k->SetEF(m_tcef);
  m_cks.push_back(k);
  m_kmap[k->Type()|(k->Type()&1?(k->Mode()?4:0):0)]
    [k->LF()->Flav(1)][k->LF()->Flav(2)]=k;
}

void Shower::SetMS(const ATOOLS::Mass_Selector *const ms)
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
  m_weightsmap["PS"].Nominal() *= cw;
  if (cv.size()) {
    ATOOLS::Reweight(m_weightsmap["PS"],
                     [&cv](double varweight,
                           size_t varindex,
                           QCD_Variation_Params& varparams) -> double {
                       return varweight * cv[varindex];
                     });
  }
  msg_Debugging()<<a<<" t = "<<t<<" -> w = "<<cw
		 <<" ("<<m_weightsmap["PS"].Nominal()<<"), v = "<<cv<<"\n";
}

int Shower::Evolve(Amplitude& a, unsigned int& nem)
{
  DEBUG_FUNC(this);
  m_weightsmap.Clear();
  m_weightsmap["PS"] = MakeWeights(Variations_Type::qcd);
  m_weightsmap["PS_QCUT"] = MakeWeights(Variations_Type::qcut);
  msg_Debugging()<<a<<"\n";
  Cluster_Amplitude *ampl(a.ClusterAmplitude());
  if (ampl->NLO()&128) {
    msg_Debugging()<<"UNLOPS veto (t<t_0)\n";
    a.Reduce();
    return 1;
  }
  if (ampl->NLO()&(4|8)) {
    msg_Debugging()<<"NLO "<<ampl->NLO()<<" path, skip shower\n";
    return 1;
  }
  if (nem>=m_maxem) {
    if (ampl->NLO()&32) {
      msg_Debugging()<<"UNLOPS sign flip\n";
      a.Reduce();
    }
    return 1;
  }
  for (Splitting s(GeneratePoint(a,a.T(),nem));
       s.m_t>Max(a.T0(),m_tmin[s.m_type&1]);
       s=GeneratePoint(a,s.m_t,nem)) {
    for (size_t i(0);i<a.size();++i) a[i]->Store();
    int stat(s.p_sk->Construct(s,1));
    msg_IODebugging()<<"t = "<<s.m_t<<", w = "<<s.m_w.MC()
		     <<" / "<<s.m_w.Accept()<<" -> "
		     <<(stat==1?"accept\n":"reject\n");
    msg_Debugging()<<"stat = "<<stat<<"\n";
    double jcv {stat ? 0.0 : -1.0};
    Jet_Finder *jf(ampl->JF<Jet_Finder>());
    if (stat && jf) {
      Cluster_Amplitude *ampl(a.GetAmplitude());
      jcv = jf->JC()->Value(ampl);
      ampl->Delete();
    }
    if (ampl->Flag()&2) {
      msg_Debugging()<<"Skip UNLOPS veto\n";
      if (s.p_l) a.Remove(s.p_l);
      a.Remove(s.p_n);
      s.p_c->SetFlav(s.p_sk->LF()->Flav(0));
      for (size_t i(0);i<a.size();++i) a[i]->Restore();
      ampl->SetFlag(ampl->Flag()&~2);
      continue;
    }
    if (ampl->NLO()&64) {
      msg_Debugging()<<"UNLOPS projection veto\n";
      if (s.p_l) a.Remove(s.p_l);
      a.Remove(s.p_n);
      s.p_c->SetFlav(s.p_sk->LF()->Flav(0));
      for (size_t i(0);i<a.size();++i) a[i]->Restore();
      a.Reduce();
      s.m_t=a.T0();
      return 1;
    }
    if (ampl->NLO()&16) {
      msg_Debugging()<<"UNLOPS veto\n";
      return 0;
    }
    const bool is_jcv_positive {jcv >= 0.0};
    bool all_vetoed {true};
    const int nqcuts = s_variations->Size(Variations_Type::qcut);
    std::vector<bool> skips (nqcuts + 1, false);
    int nskips {0};
    ATOOLS::ReweightAll(
        m_weightsmap["PS_QCUT"],
        [this, jcv, is_jcv_positive, ampl, &all_vetoed, &skips, &nskips](
            double varweight,
            size_t varindex,
            Qcut_Variation_Params* qcutparams) -> double {
          msg_Debugging() << "Applying veto weight to qcut var #" << varindex
                          << " {\n";
          bool stat {is_jcv_positive};
          Jet_Finder* jf(ampl->JF<Jet_Finder>());
          if (stat && jf) {
            const double fac {
                qcutparams == nullptr ? 1.0 : qcutparams->m_scale_factor};
            stat = jcv < sqr(jf->Qcut() * fac);
            msg_Debugging() << "  jcv = " << sqrt(jcv) << " vs "
                            << jf->Qcut() << " * " << fac << " = "
                            << jf->Qcut() * fac << "\n";
          }
          if (stat) {
            msg_Debugging() << "} no jet veto\n";
            all_vetoed = false;
            return varweight;
          } else if (ampl->NLO() & 2) {
            msg_Debugging() << "  skip emission\n";
            skips[varindex] = true;
            ++nskips;
            all_vetoed = false;
            msg_Debugging() << "} no jet veto\n";
            return varweight;
          } else {
            msg_Debugging() << "} jet veto\n";
            return 0.0;
          }
        });
    if (ampl->NLO()&2) {
      const int nqcdvars = s_variations->Size(Variations_Type::qcd);
      if (skips[0])
        nskips += nqcdvars;
      const double wskip {nskips / double(nqcuts + nqcdvars + 1)};
      if (ran->Get()<=wskip) {
	if (s.p_l) a.Remove(s.p_l);
	a.Remove(s.p_n);
	s.p_c->SetFlav(s.p_sk->LF()->Flav(0));
	for (size_t i(0);i<a.size();++i) a[i]->Restore();
	double lkf(ampl->LKF());
	for (Cluster_Amplitude *campl(ampl);
	     campl;campl=campl->Prev()) {
	  campl->SetLKF(1.0);
	  ampl->SetNLO(ampl->NLO()&~2);
	}
        const double fac {1.0 / lkf / wskip};
        m_weightsmap["PS"] *= fac * skips[0];
        m_weightsmap["PS_QCUT"] *= skips;
	continue;
      }
      else {
        const double fac {1.0 / (1.0 - wskip)};
        skips.flip();
        m_weightsmap["PS"] *= fac * skips[0];
        m_weightsmap["PS_QCUT"] *= skips;
      }
    }
    if (all_vetoed)
      return 0;
    AddWeight(a,s.m_t);
    a.SetJF(NULL);
    if (++nem>=m_maxem) break;
    if (a.size()-a.ClusterAmplitude()->NIn()>m_maxpart) {
      if (s.p_l) a.Remove(s.p_l);
      a.Remove(s.p_n);
      s.p_c->SetFlav(s.p_sk->LF()->Flav(0));
      for (size_t i(0);i<a.size();++i) a[i]->Restore();
    }
    if (a.size()-a.ClusterAmplitude()->NIn()>=m_maxpart) break;
  }
  AddWeight(a,a.T0());
  if (ampl->NLO()&32) {
    msg_Debugging()<<"UNLOPS sign flip\n";
    a.Reduce();
  }
  return 1;
}

Splitting Shower::GeneratePoint
(const Amplitude &a,const double &t,const unsigned int &nem)
{
  Splitting win;
  double tmin[2]={m_tmin[0],m_tmin[1]};
  for (Amplitude::const_reverse_iterator
	 it(a.rbegin());it!=a.rend();++it) {
    for (int cm(0);cm<2;++cm) {
      double ct((*it)->T(cm)>=0.0?(*it)->T(cm):t);
      Splitting cur(GeneratePoint(**it,ct,cm,nem));
      (*it)->SetT(cm,-1.0);
      if (cur.p_c==NULL || cur.p_s==NULL) continue;
      if (cur.m_t<m_tmin[cur.m_type&1]) continue;
      m_tmin[0]=m_tmin[1]=(win=cur).m_t;
    }
  }
  m_tmin[0]=tmin[0];
  m_tmin[1]=tmin[1];
  if (win.p_sk && win.m_t>m_tmin[win.m_type&1])
    msg_Debugging()<<"Emission at "<<win<<"\n";
  return win;
}

Splitting Shower::GeneratePoint
(Parton &p,const double &t,const int &cm,const unsigned int &nem)
{
  Splitting win(&p,NULL,t);
  double sum=0.0, ct=m_rcf*t;
  SKernel_Map::const_iterator kit(m_sks.find(p.Flav()));
  if (kit==m_sks.end()) return Splitting();
  const int nkernels = kit->second.size();
  const int nspects = static_cast<int>(p.Ampl()->size());
  if (nkernels > m_sums.nrows || nspects > m_sums.ncols) {
    m_sums = CumulativeIntegralTable {nkernels, nspects};
  }
  else {
    m_sums.nrows = nkernels;
    m_sums.ncols = nspects;
  }
  while (true) {
    if (win.m_t*m_rcf<=ct) {
      sum=0.0;
      ct=win.m_t;
      for (size_t j(0);j<nkernels;++j) {
        m_sums.Clear(j);
	double csum=0.0;
	for (int i(nspects-1);i>=0;--i) {
	  if ((*p.Ampl())[i]==&p) continue;
	  Splitting cur(&p,(*p.Ampl())[i]);
	  cur.SetType();
	  cur.m_kfac=m_kfac;
	  cur.m_cpl=m_cpl;
	  cur.m_t1=ct;
	  cur.m_cm=cm;
	  if (kit->second[j]->On() &&
	      kit->second[j]->Allowed(cur)) {
            csum += dabs(kit->second[j]->Integral(cur));
            m_sums.AppendSumAndSpect(j, csum, i);
	  }
	}
        if (m_sums.Size(j)) sum += m_sums.LastSum(j);
      }
      if (sum==0.0) return Splitting();
      win=Splitting(&p,NULL,ct);
    }
    win.m_t*=exp(log(ran->Get())*2.0*M_PI/sum);
    if (win.m_t<m_tmin[p.Beam()?1:0]) return win;
    double disc(sum*ran->Get()), csum(0.0);
    for (size_t j(0);j<nkernels;++j)
      if (m_sums.Size(j) && (csum+=m_sums.LastSum(j)) >= disc) {
	double disc(m_sums.LastSum(j)*ran->Get());
	for (size_t i(0);i<m_sums.Size(j);++i)
	  if (m_sums.Sum(j, i) >= disc) {
	    win.p_s=(*p.Ampl())[m_sums.Spect(j, i)];
	    win.SetType();
	    win.m_kin=m_kin;
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
	    if (kit->second[j]->LF()->Construct(win,0)!=1) break;
	    win.m_w=kit->second[j]->GetWeight(win,m_oef);
            win.m_vars = std::vector<double>(s_variations->Size(), 1.0);
            if (win.m_w.MC() < ran->Get()) {
              if (m_reweight && win.m_t > m_rewtmin) {
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
            if (m_reweight && win.m_t > m_rewtmin) {
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
  m_rsf*=params->m_showermuR2fac;
  m_fsf*=params->m_showermuF2fac;
  MODEL::Running_AlphaS *as(p_as);
  p_as=params->p_alphas;
  PDF::PDF_Base *pdf[2]={p_pdf[0],p_pdf[1]};
  p_pdf[0]=params->p_pdf1;
  p_pdf[1]=params->p_pdf2;
  if (rsf==m_rsf && fsf==m_fsf && as==p_as &&
      pdf[0]==p_pdf[0] && pdf[1]==p_pdf[1]) {
    a.m_s->m_vars[varindex] = 1.0;
    return;
  }
  msg_IODebugging()<<METHOD<<"("<<varindex<<") {\n  "
		   <<"\\mu_R -> "<<sqrt(m_rsf)
		   <<", \\mu_F -> "<<sqrt(m_fsf)<<"\n  PDF "
		   <<(p_pdf[0]?p_pdf[0]->LHEFNumber():-1)<<" x "
		   <<(p_pdf[1]?p_pdf[1]->LHEFNumber():-1)<<"\n";
  Weight w(a.m_s->p_sk->GetWeight(*a.m_s,m_oef,&a.m_s->m_w));
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
  if (x<p_pdf[b]->XMin() ||
      x>p_pdf[b]->XMax()*p_pdf[b]->RescaleFactor() ||
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

int Shower::RemnantTest(Parton *const c,const Vec4D &p)
{
  Vec4D pb(rpa->gen.PBeam(c->Beam()-1));
  if (p[0]<0.0 || p.Nan()) return -1;
  if (p[0]>pb[0] && !IsEqual(p[0],pb[0],1.0e-6)) return -1;
  return 1;
}
