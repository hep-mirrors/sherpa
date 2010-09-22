#include "PHASIC++/Selectors/Jet_Finder.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;

Jet_Finder::Jet_Finder
(const int nin,const int nout,Flavour *fl,
 const std::string &ycut):
  Selector_Base("Jetfinder"), m_dparam(0.3), m_cuttag(ycut),
  m_on(true), p_yccalc(NULL)
{
  m_ycut=2.0;
  m_fl=fl;
  m_nin=nin;
  m_nout=nout;
  m_n=m_nin+m_nout;
  m_smax=m_s=sqr(rpa.gen.Ecms());
  if (ycut.find("|")!=std::string::npos) {
    m_dparam=ToType<double>(ycut.substr(ycut.find("|")+1));
    m_cuttag=ycut.substr(0, ycut.find("|"));
  }
  m_pt2min=1.0;
  m_sel_log = new Selector_Log(m_name);
  static bool mets(false);
  if (!mets) {
    mets=true;
    rpa.gen.AddCitation(1,"Matrix element merging with truncated showers is "+
			std::string("published under \\cite{Hoeche:2009rj}."));
  }
}

Jet_Finder::~Jet_Finder() 
{
  if (p_yccalc) delete p_yccalc;
}

size_t Jet_Finder::FillCombinations(const Subprocess_Info &subinfo,
				    size_t &cp,const int fl)
{
  // sum is the sum of the binary indices of all particles which are at the
  // current level of decay-building, e.g. for a simple 2->2 it is going to be
  // 1 + 2 + 4 + 8 or in a decay step it is the sum of all decay products.
  size_t sum(0);
  // pos contains all binary indices involved in the current level of
  // decay-building
  std::vector<int> pos;

  if (subinfo.m_ps.size()==0) return 1<<cp++;

  for (size_t i(0); i<subinfo.m_ps.size(); ++i) {
    pos.push_back(FillCombinations(subinfo.m_ps[i],cp,fl-1));
    if (pos.back()&((1<<m_nin)-1)) 
      m_flavs[pos.back()]=subinfo.m_ps[i].m_fl.Bar();
    else m_flavs[pos.back()]=subinfo.m_ps[i].m_fl;
    sum=sum|pos.back();
  }

  // now that we know about all outer and intermediate particles involved in the
  // process, loop over them and identify all combinations on which we want to
  // cut
  for (size_t i(0);i<pos.size();++i) {
    // a safety pT cut of 1GeV will be added to FS gluons/photons
    // this is purely for technical reasons, no physics and has to do with the
    // choice of reference vectors for their pol vectors in comix
    if ((pos[i]&3)==0 && (m_flavs[pos[i]].Kfcode()==kf_photon ||
                          m_flavs[pos[i]].Kfcode()==kf_gluon
                          || m_flavs[pos[i]].Strong()
                          )) {
      m_pcs.insert(pos[i]);
    }
    
    for (size_t j(i+1);j<pos.size();++j) {
      if (pos[i]>2 || pos[j]>2) { // we can ignore the initial-initial case
        // for all others we check whether they are an allowed combination
        // and if yes, add them to the appropriate decay level "fl"
        
        // check whether i and j are actually combinable in the ME, otherwise
        // they shouldn't be cut on, even if they could be produced by the PS
        // example: e+ e- -> G G d db   here d db are *not* combinable!
        if (p_sproc && p_sproc->Combinable(pos[i], pos[j])) {
          for (size_t k(0);k<pos.size();++k) {
            if (i==k || j==k) continue;
            // un-bar IS flavs for looking for kernel, since that is taken care
            // of by the splitting function type
            Flavour fli((pos[i]&3)?m_flavs[pos[i]].Bar():m_flavs[pos[i]]);
            Flavour flj((pos[j]&3)?m_flavs[pos[j]].Bar():m_flavs[pos[j]]);
            Flavour flk((pos[k]&3)?m_flavs[pos[k]].Bar():m_flavs[pos[k]]);
            // splitting function type:
            // 1 <-> emitter im IS, 2 <-> spectator im IS => 0=FF, 1=IF, 2=FI, 3=II
            int sftype((pos[i]+pos[j])&3?(pos[k]&3)?3:1:(pos[k]&3)?2:0);

            // ask whether the combination i j k is allowed by the shower and if
            // so, whether it's coupling is strong (1), ew (2) or both (3)
            int cpl(p_proc->Process()->Shower()->HasKernel(fli, flj, flk,sftype));
            if (cpl>0) {
              m_fills[fl].push_back(Comb_Key(cpl, pos[i],pos[j],pos[k]));
            }
          }
        }
      }
    }
  }
  m_mcomb.push_back(pos);
  m_mcomb.back().push_back(sum);
  return sum;
}

void Jet_Finder::FillCombinations()
{
  if (p_yccalc==NULL) {
    if (p_proc==NULL) THROW(fatal_error,"Process not set.");
    p_sproc=p_proc->Process()->Get<Single_Process>();
    m_moms.clear();
    m_flavs.clear();
    m_fills.resize(m_nin+m_nout+1);
    Subprocess_Info start(p_proc->Process()->Info().m_ii);
    start.Add(p_proc->Process()->Info().m_fi);
    size_t idx(0);
    FillCombinations(start, idx,m_nin+m_nout);
    p_yccalc = new Algebra_Interpreter();
    p_yccalc->SetTagReplacer(this);
    PrepareMomList(p_proc->Momenta());
    for (size_t i=0;i<p_proc->NIn()+p_proc->NOut();++i) 
      p_yccalc->AddTag("p["+ToString(i)+"]",ToString(p_proc->Momenta()[i]));
    p_yccalc->Interprete(m_cuttag);
    if (msg_LevelIsDebugging()) {
      msg_Out()<<METHOD<<"(): Q_{cut} for '"
	       <<p_proc->Process()->Name()<<"' {\n";
      {
	msg_Indent();
	p_yccalc->PrintEquation();
      }
      msg_Out()<<"}\n";
      msg_Out()<<METHOD<<"(): Identified clusterings {\n";
      for (size_t i(0);i<m_fills.size();++i)
        for (size_t j(0);j<m_fills[i].size();++j)
          msg_Out()<<"  "<<i<<": ["<<ID(m_fills[i][j].first)<<","
                   <<ID(m_fills[i][j].second)<<"] <-> "
                   <<ID(m_fills[i][j].partner)<<" => "
                   <<m_flavs[m_fills[i][j].first]<<" & "
                   <<m_flavs[m_fills[i][j].second]<<" <-> "
                   <<m_flavs[m_fills[i][j].partner]<<"\n";
      msg_Out()<<"}\n";
      msg_Out()<<METHOD<<"(): Momentum combination {\n";
      for (size_t i(0);i<m_mcomb.size();++i) {
	msg_Out()<<"  "<<ID(m_mcomb[i].back())<<" -> {";
	for (size_t j(0);j<m_mcomb[i].size()-1;++j) 
	  msg_Out()<<" "<<ID(m_mcomb[i][j]);
	msg_Out()<<" }\n";
      }
      msg_Out()<<"}\n";
    }
  }
}

void Jet_Finder::PrepareMomList(const Vec4D_Vector &vec)
{
  for (int i(m_nin+m_nout-1);i>=0;--i) {
    m_moms[1<<i]=i<m_nin?-vec[i]:vec[i];
#ifdef DEBUG__Prepare_Moms
    msg_Debugging()<<"p["<<i<<"] = "<<m_moms[1<<i]
		   <<" ("<<m_flavs[1<<i]<<")\n";
#endif
  }
  for (size_t n(0);n<m_mcomb.size()-1;++n) {
    m_moms[m_mcomb[n].back()]=m_moms[m_mcomb[n].front()];
    for (size_t i(1);i<m_mcomb[n].size()-1;++i)
      m_moms[m_mcomb[n].back()]+=m_moms[m_mcomb[n][i]];
#ifdef BOOST_Decays
    Poincare cms(m_moms[m_mcomb[n].back()]);
    for (size_t i(0);i<m_mcomb[n].size()-1;++i) {
      cms.Boost(m_moms[m_mcomb[n][i]]);
      cc+=m_moms[m_mcomb[n][i]];
    }
    static double accu(sqrt(Accu()));
    Vec4D::SetAccu(accu);
    if (!(Vec3D(cc)==Vec3D()) || 
	!IsEqual(cc.Abs2(),m_moms[m_mcomb[n].back()].Abs2())) 
      msg_Error()<<METHOD<<"(): CMS boost failure. sum = "
		 <<cc<<" "<<cc.Abs2()<<" vs. "
		 <<m_moms[m_mcomb[n].back()].Abs2()<<"\n";
    Vec4D::ResetAccu();
#endif
#ifdef DEBUG__Prepare_Moms
    msg_Debugging()<<"p["<<ID(m_mcomb[n].back())<<"] = "
  		   <<m_moms[m_mcomb[n].back()]
  		   <<" ["<<m_flavs[m_mcomb[n].back()]<<"]\n";
#endif
  }
}

bool Jet_Finder::PrepareColList(const std::vector<int> &ci,
				const std::vector<int> &cj)
{
  for (int i(m_nin+m_nout-1);i>=0;--i) {
    m_cols[1<<i]=ColorID(ci[i],cj[i]);
#ifdef DEBUG__Prepare_Cols
    msg_Debugging()<<"c["<<i<<"] = "<<m_cols[1<<i]
		   <<" ("<<m_flavs[1<<i]<<")\n";
#endif
  }
  for (size_t n(0);n<m_mcomb.size()-1;++n) {
    int c[4]={0,0,0,0};
    for (size_t i(0);i<m_mcomb[n].size()-1;++i) {
      ++c[m_cols[m_mcomb[n][i]].m_i];
      --c[m_cols[m_mcomb[n][i]].m_j];
    }
    ColorID cc(0,0);
    for (int i(1);i<4;++i) {
      if (c[i]==1 && cc.m_i==0) cc.m_i=i;
      else if (c[i]==-1 && cc.m_j==0) cc.m_j=i;
      else if (c[i]!=0) {
	msg_Debugging()<<METHOD<<"(): Rejecting invalid intermediate color.\n";
	return false;
      }
    }
    m_cols[m_mcomb[n].back()]=cc;
#ifdef DEBUG__Prepare_Cols
    msg_Debugging()<<"c["<<ID(m_mcomb[n].back())<<"] = "
  		   <<m_cols[m_mcomb[n].back()]
  		   <<" ["<<m_flavs[m_mcomb[n].back()]<<"] "<<cc<<"\n";
#endif
  }
  return true;
}

bool Jet_Finder::Trigger(const Vec4D_Vector &p)
{
  FillCombinations();
  PrepareMomList(p);
  m_ycut=p_yccalc->Calculate()->Get<double>();
  if (!m_on) return true;
  msg_Debugging()<<METHOD<<"(): '"<<p_proc->Process()->Name()
		 <<"' Q_cut = "<<sqrt(m_ycut*m_s)<<(m_on?" {":", off")<<"\n";
  bool uc(false);
  SP(Color_Integrator) ci(p_proc->ColorIntegrator());
  if (ci!=NULL && ci->On()) {
    uc=true;
    std::vector<int> ic(ci->I()), jc(ci->J());
    if (!PrepareColList(ic,jc)) return 1-m_sel_log->Hit(true);
  }
  for (std::set<size_t>::const_iterator it(m_pcs.begin());it!=m_pcs.end();++it){
    const size_t i=*it;
    msg_Debugging()<<"  "<<ID(i)<<"["<<m_flavs[i]
		   <<"] -> "<<m_moms[i].PPerp()<<"\n";
    if (m_moms[i].PPerp2()<m_pt2min) return 1-m_sel_log->Hit(true);
  }
  for (size_t cl(1);cl<m_fills.size();++cl) {
    if (m_fills[cl].empty()) continue;
    msg_Indent();
    msg_Debugging()<<"level = "<<m_fills.size()-cl<<" {\n";
    for (size_t ps(0);ps<m_fills[cl].size();++ps) {
      size_t i(m_fills[cl][ps].first),
	j(m_fills[cl][ps].second), k(m_fills[cl][ps].partner);
      int cpl = m_fills[cl][ps].cpl;
      if (uc && (cpl&1) && !ColorConnected(i,j,k)) continue;
      if (m_flavs[i].IsQuark() && m_flavs[j].IsQuark() &&
	  m_flavs[i]!=m_flavs[j].Bar()) continue;
      msg_Debugging()<<"  "<<ID(i)<<"["<<m_flavs[i]<<"] & "
		     <<ID(j)<<"["<<m_flavs[j]<<"] <-> "
		     <<ID(k)<<"["<<m_flavs[k]<<"], qcut = "
		     <<sqrt(m_ycut*m_s);
      double pt2ij=Qij2(m_moms[i],m_moms[j],m_moms[k],m_flavs[i],m_flavs[j],
                        m_dparam);
      msg_Debugging()<<", ptjk = "<<sqrt(pt2ij)<<" ("
		     <<(pt2ij>=m_ycut*m_s)<<(pt2ij<m_ycut*m_s?")\n":")");
      if (pt2ij<m_ycut*m_s) return 1-m_sel_log->Hit(true);
      msg_Debugging()<<"\n";
    }
    msg_Debugging()<<"}\n";
  }
  msg_Debugging()<<"}\n";
  return 1-m_sel_log->Hit(false);
}

struct Qij2_Key {
public:

  double m_qij2;
  size_t m_i, m_j, m_k;

public:

  inline Qij2_Key(const double &qij2=std::numeric_limits<double>::max(),
		  const size_t &i=0,const size_t &j=0,const size_t &k=0):
    m_qij2(qij2), m_i(i), m_j(j), m_k(k) {}

};// end of struct Qij2_Key

bool Jet_Finder::JetTrigger(const ATOOLS::Vec4D_Vector &p,
                            NLO_subevtlist *const subs)
{
  FillCombinations();
  PrepareMomList(p);
  m_ycut=p_yccalc->Calculate()->Get<double>();
  if (!m_on) return true;
  msg_Debugging()<<METHOD<<"("<<p_proc->Process()->Name()<<"): {\n";
  bool uc(false);
  SP(Color_Integrator) ci(p_proc->ColorIntegrator());
  if (ci!=NULL && ci->On()) {
    uc=true;
    std::vector<int> ic(ci->I()), jc(ci->J());
    if (!PrepareColList(ic,jc)) return 1-m_sel_log->Hit(true);
  }
  Qij2_Key qij2;
  for (size_t cl(1);cl<m_fills.size();++cl) {
    if (m_fills[cl].empty()) continue;
    msg_Indent();
    msg_Debugging()<<"level = "<<m_fills.size()-cl<<" {\n";
    for (size_t ps(0);ps<m_fills[cl].size();++ps) {
      size_t i(m_fills[cl][ps].first),
	j(m_fills[cl][ps].second), k(m_fills[cl][ps].partner);
      int cpl = m_fills[cl][ps].cpl;
      if (uc && (cpl&1) && !ColorConnected(i,j,k)) continue;
      if (m_flavs[i].IsQuark() && m_flavs[j].IsQuark() &&
	  m_flavs[i]!=m_flavs[j].Bar()) continue;
      double pt2ij=Qij2(m_moms[i],m_moms[j],m_moms[k],
			kf_gluon,kf_gluon);
      msg_Debugging()<<"  "<<ID(i)<<"["<<m_flavs[i]<<"] & "
		     <<ID(j)<<"["<<m_flavs[j]<<"] <-> "
		     <<ID(k)<<"["<<m_flavs[k]
		     <<"], ptjk = "<<sqrt(pt2ij)<<"\n";
      if (pt2ij<qij2.m_qij2) qij2=Qij2_Key(pt2ij,i,j,k);
    }
    msg_Debugging()<<"}\n";
    break;
  }
  if (qij2.m_i==qij2.m_j)
    THROW(fatal_error,"No valid clustering");
  size_t i(ID(qij2.m_i).front()), j(ID(qij2.m_j).front());
  size_t k(ID(qij2.m_k).front());
  NLO_subevt *sub(NULL);
  for (size_t l(0);l<subs->size();++l) {
    sub=(*subs)[l];
    if ((sub->m_i==i && sub->m_j==j && sub->m_k==k) ||
	(sub->m_i==j && sub->m_j==i && sub->m_k==k)) break;
    else sub=NULL;
  }
  if (sub==NULL) THROW(fatal_error,"Internal error");
  Process_Base *pi(sub->Proc<Process_Integrator>()->Process());
  Jet_Finder *jf((Jet_Finder*)pi->Selector()->GetSelector("Jetfinder"));
  if (jf==NULL) THROW(critical_error,"No jet finder for "+pi->Name());
  msg_Debugging()<<"} -> jf = "<<jf<<" => '"<<pi->Name()<<"'\n"<<*sub<<"\n";
  Vec4D_Vector pp(CombineMoms(p,sub));
  Vec4D sum;
  for (size_t l(0);l<pp.size();++l)
    sum+=l<(size_t)m_nin?-pp[l]:pp[l];
  if (!IsEqual(sum,Vec4D(),1.0e-6)) {
    msg_Error()<<METHOD<<"(): Momentum not conserved ("
	       <<i<<","<<j<<") <-> "<<k<<" {\n";
    for (size_t l(0);l<pp.size();++l)
      msg_Error()<<"  p["<<l<<"] = "<<pp[l]<<"\n";
    msg_Error()<<"}"<<std::endl;
  }
  bool trg=jf->Trigger(pp);
  return trg;
}

double Jet_Finder::Qij2Min(const ATOOLS::Vec4D_Vector &p,
			   NLO_subevtlist *const subs)
{
  DEBUG_FUNC(subs->back()->m_pname);
  Qij2_Key qij2;
  const ATOOLS::Flavour *f(subs->back()->p_fl);
  for (size_t n(0);n<subs->size()-1;++n) {
    size_t i((*subs)[n]->m_i), j((*subs)[n]->m_j), k((*subs)[n]->m_k);
    Vec4D pi(i<2?-p[i]:p[i]), pj(j<2?-p[j]:p[j]), pk(k<2?-p[k]:p[k]);
    Flavour fi(i<2?f[i].Bar():f[i]), fj(j<2?f[j].Bar():f[j]);
    double pt2ij=Qij2(pi,pj,pk,fi,fj);
    msg_Debugging()<<"("<<i<<")["<<fi<<"] & ("<<j<<")["<<fj
		   <<"] <-> ("<<k<<"), ptjk = "<<sqrt(pt2ij)<<"\n";
    if (pt2ij<qij2.m_qij2) qij2=Qij2_Key(pt2ij,i,j,k);
  }
  msg_Debugging()<<"q_ij = "<<sqrt(qij2.m_qij2)<<"\n";
  return qij2.m_qij2;
}

double Jet_Finder::Lam
(const double &s,const double &sb,const double &sc) const
{
  return sqr(s-sb-sc)-4.0*sb*sc;
}

Vec4D_Vector Jet_Finder::CombineMoms
(const Vec4D_Vector &p,NLO_subevt *const sub)
{
  size_t i(sub->m_i), j(sub->m_j), k(sub->m_k);
  if (i>j) std::swap<size_t>(i,j);
  Flavour flij(kf_none);
  for (size_t m(0);m<sub->m_n;++m)
    if (ID(sub->p_id[m]).size()==2) {
      flij=sub->p_fl[m];
      break;
    }
  double mij2=sqr(flij.Mass()), mk2=sqr(m_flavs[1<<k].Mass());
  if (i>1 && j>1 && k>1) {
    Vec4D pi(p[i]), pj(p[j]), pk(p[k]), Q(pi+pj+pk);
    double Q2=Q.Abs2(), lrat=Lam(Q2,mij2,mk2)/Lam(Q2,(pi+pj).Abs2(),mk2);
    Vec4D pkt(sqrt(lrat)*(pk-(Q*pk/Q2)*Q)+(Q2+mk2-mij2)/(2.*Q2)*Q);
    Vec4D pijt(Q-pkt);
    if (lrat<0.0 || pkt[0]<0.0 || pijt[0]<0.0) return Vec4D_Vector();
    return GetMomMap(p,pijt,pkt,sub,i);
  }
  if (i>1 && j>1 && k<2) {
    Vec4D pi(p[i]), pj(p[j]), pa(-p[k]), Q(pa+pi+pj);
    double Q2=Q.Abs2(), lrat=Lam(Q2,mij2,mk2)/Lam(Q2,(pi+pj).Abs2(),mk2);
    Vec4D pat(sqrt(lrat)*(pa-(Q*pa/Q2)*Q)+(Q2+mk2-mij2)/(2.*Q2)*Q);
    Vec4D pijt(Q-pat), pb(-p[1-k]);
    if (lrat<0.0 || pat[0]>0.0 || pijt[0]<0.0) return Vec4D_Vector();
    Vec4D_Vector pp(GetMomMap(p,pijt,pat,sub,i));
    ZAlign(pp,pat,pb,mk2,sqr(m_flavs[1<<(1-k)].Mass()));
    return pp;
  }
  if (i<2 && j>1 && k>1) {
    Vec4D pa(-p[i]), pj(p[j]), pk(p[k]), Q(pa+pj+pk);
    double Q2=Q.Abs2(), lrat=Lam(Q2,mij2,mk2)/Lam(Q2,(pa+pj).Abs2(),mk2);
    Vec4D pkt(sqrt(lrat)*(pk-(Q*pk/Q2)*Q)+(Q2+mk2-mij2)/(2.*Q2)*Q);
    Vec4D pajt(Q-pkt), pb(-p[1-i]);
    if (lrat<0.0 || pkt[0]<0.0 || pajt[0]>0.0) return Vec4D_Vector();
    Vec4D_Vector pp(GetMomMap(p,pajt,pkt,sub,i));
    ZAlign(pp,pajt,pb,mij2,sqr(m_flavs[1<<(1-i)].Mass()));
    return pp;
  }
  if (i<2 && j>1 && k<2) {
    Vec4D pa(-p[i]), paj(pa+p[j]), pb(-p[k]), Q(paj+pb);
    double Q2=Q.Abs2(), saj=paj.Abs2();
    Vec4D_Vector pp(p);
    ZAlign(pp,paj,pb,saj,mk2,1);
    Q=paj+pb;
    double lrat=Lam(Q2,mij2,mk2)/Lam(Q2,saj,mk2);
    Vec4D pbt(sqrt(lrat)*(pb-(Q*pb/Q2)*Q)+(Q2+mk2-mij2)/(2.*Q2)*Q);
    if (lrat<0.0 || pbt[0]>0.0) return Vec4D_Vector();
    return GetMomMap(pp,Q-pbt,pbt,sub,i);
  }
  return Vec4D_Vector();
}

Vec4D_Vector Jet_Finder::GetMomMap
(const Vec4D_Vector &p,const Vec4D &pijt,const Vec4D &pkt,
 NLO_subevt *const sub,const size_t &i) const
{
  Vec4D_Vector pp(p.size()-1);
  for (size_t m(0);m<sub->m_n;++m) {
    size_t l(ID(sub->p_id[m]).front());
    if (l==i) pp[m]=pijt;
    else if (l==sub->m_k) pp[m]=pkt;
    else pp[m]=p[l];
  }
  return pp;
}

void Jet_Finder::ZAlign
(ATOOLS::Vec4D_Vector &p,ATOOLS::Vec4D &pa,
 const ATOOLS::Vec4D &pb,const double &ma2,const double &mb2,
 const int mode) const
{
  Vec4D Q(pa+pb);
  double Q2=Q.Abs2(), papb=pa*pb, sb=Sign(pb[3]), ea=0.0;
  if (IsZero(mb2)) ea=0.5*(papb+ma2*sqr(pb[3])/papb)/pb[0];
  else ea=(pb[0]*papb+dabs(pb[3])*sqrt(papb*papb-ma2*mb2))/mb2;
  Vec4D pan(ea,0.0,0.0,-sb*sqrt(ea*ea-ma2)), pam(ea,0.0,0.0,-pan[3]);
  if (dabs((pam+pb).Abs2()-Q2)<dabs((pan+pb).Abs2()-Q2)) pan=pam;
  Poincare cmso(-Q), cmsn(-pan-pb);
  cmso.Boost(pa);
  Poincare zrot(pa,-sb*Vec4D::ZVEC);
  for (size_t m(0);m<p.size();++m) {
    cmso.Boost(p[m]);
    zrot.Rotate(p[m]);
    cmsn.BoostBack(p[m]);
  }
}

bool Jet_Finder::NoJetTrigger(const ATOOLS::Vec4D_Vector &p)
{
  return true;
}

double Jet_Finder::Qij2(const Vec4D &pi,const Vec4D &pj,const Vec4D &pk,
			const Flavour &fi,const Flavour &fj,
                        const double &dparam,
                        const int mode)
{
  Vec4D npi(pi), npj(pj);
  Flavour nfi(fi), nfj(fj);
  if (npi[0]<0.0) {
    npi=-pi-pj;
    if (mode==0) nfi=fi==fj.Bar()?Flavour(kf_gluon):(fi.IsVector()?fj.Bar():fi);
  }
  else if (npj[0]<0.0) {
    npj=-pj-pi;
    if (mode==0) nfj=fj==fi.Bar()?Flavour(kf_gluon):(fj.IsVector()?fi.Bar():fj);
  }
  if (nfi.IsQuark() && nfj.IsQuark() && nfi!=nfj.Bar()) return -1.0;
  if ((fi.IsPhoton() && fj.IntCharge()==0) ||
      (fj.IsPhoton() && fi.IntCharge()==0)) return -1.0;
  if (fi.IsPhoton() || fj.IsPhoton()) {
    if (pi[0]<0.0) return pj.PPerp2();
    if (pj[0]<0.0) return pi.PPerp2();
    return Min(pi.PPerp2(), pj.PPerp2())*sqr(pi.DR(pj)/dparam);
  }
  double pipj(dabs(npi*npj)), pipk(dabs(npi*pk)), pjpk(dabs(npj*pk));
  double mti(sqr(Flavour(nfi).Mass())), mtj(sqr(Flavour(nfj).Mass()));
  if (pipj==0.0) {
    if (mti!=0.0 || mtj!=0.0) THROW(fatal_error,"Ill-defined mass term");
  }
  else {
    mti/=2.0*pipj;
    mtj/=2.0*pipj;
  }
  double Cij(nfj.IsVector()?Max(0.0,pipk/(pipj+pjpk)-mti):0.5*pipk/(pipk+pjpk));
  double Cji(nfi.IsVector()?Max(0.0,pjpk/(pipj+pipk)-mtj):0.5*pjpk/(pipk+pjpk));
  return 2.0*dabs(pi*pj)/(Cij+Cji);
}

bool Jet_Finder::ColorConnected
(const size_t &i,const size_t &j,const size_t &k) const
{
  const ColorID &ci(m_cols.find(i)->second), 
    &cj(m_cols.find(j)->second);
  int si(m_flavs.find(i)->second.StrongCharge()), 
    sj(m_flavs.find(j)->second.StrongCharge());
  if (!ColorConnected(ci,cj,si,sj)) return false;
  const ColorID &ck(m_cols.find(k)->second);
  int sk(m_flavs.find(k)->second.StrongCharge());
  const Flavour_Vector &cf(p_sproc->CombinedFlavour(i+j));
  for (size_t f(0);f<cf.size();++f) {
    int sij(cf[f].StrongCharge());
    if (sij==0) {
      if (abs(sk)==3) return true;
      continue;
    }
    ColorID cij(0,0);
    if (sij==3) cij.m_i=si==3?cj.m_i:ci.m_i;
    else if (sij==-3) cij.m_j=si==-3?cj.m_j:ci.m_j;
    else {
      if (ci.m_i==cj.m_j) {
	if (ci.m_j==cj.m_i) {
	  cij=ColorID(ci.m_i,cj.m_j);
	  if (ColorConnected(cij,ck,sij,sk)) return true;
	}
	cij=ColorID(cj.m_i,ci.m_j);
      }
      else {
	cij=ColorID(ci.m_i,cj.m_j);
      }
    }
    if (ColorConnected(cij,ck,sij,sk)) return true;
  }
  return false;
}

bool Jet_Finder::ColorConnected
(const ColorID &ci,const ColorID &cj,const int si,const int sj) const
{
  if (si==3) {
    if (sj==8) 
      if (cj.m_i!=cj.m_j && ci.m_i!=cj.m_j) return false;
  }
  else if (si==-3) {
    if (sj==8) 
      if (cj.m_j!=cj.m_i && ci.m_j!=cj.m_i) return false;
  }
  else {
    if (sj==3) {
      if (ci.m_j!=ci.m_i && ci.m_j!=cj.m_i) return false;
    }
    else if (sj==-3) {
      if (ci.m_i!=ci.m_j && ci.m_i!=cj.m_j) return false;
    }
    else {
      if ((ci.m_i!=cj.m_j && ci.m_j!=cj.m_i) ||
	  (ci.m_i==cj.m_j && ci.m_j==cj.m_i && ci.m_i==ci.m_j)) return false;
    }
  }
  return true;
}

void Jet_Finder::BuildCuts(Cut_Data *cuts) 
{
  FillCombinations();
  if (!m_on) return;
  msg_Debugging()<<METHOD<<"(): {\n";
  for (std::set<size_t>::const_iterator
	 it(m_pcs.begin());it!=m_pcs.end();++it) {
    const size_t id=*it, i=ID(id).front();
    msg_Debugging()<<"  p_{T,"<<i<<"} > "<<sqrt(m_pt2min)<<"\n";
    cuts->etmin[i]=Max
      (cuts->etmin[i],sqrt(m_pt2min+sqr(m_flavs[id].SelMass())));
  }
  for (int i(m_nin); i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = m_fl[i].Mass();
    if (m_fl[i].Resummed()) {
      for (int j(i+1); j<m_nin+m_nout; ++j) {
	if (m_fl[j].Resummed()) {
	  double scut=Max(1.0,sqr(m_fl[i].Mass())+sqr(m_fl[j].Mass()));
	  msg_Debugging()<<"  m_{"<<i<<","<<j<<"} > "<<sqrt(scut)<<"\n";
	  cuts->scut[i][j]=cuts->scut[j][i]=Max(cuts->scut[i][j],scut);
	}
      }
    }
  }
  msg_Debugging()<<"}\n";
}

std::string Jet_Finder::ReplaceTags(std::string &expr) const
{
  return p_yccalc->ReplaceTags(expr);
}

Term *Jet_Finder::ReplaceTags(Term *term) const
{
  term->Set(m_moms.find(1<<term->Id())->second);
  return term;
}

void Jet_Finder::AssignId(Term *term)
{
  term->SetId(ToType<int>
	      (term->Tag().substr
	       (2,term->Tag().length()-3)));
}

namespace PHASIC{

DECLARE_ND_GETTER(Jet_Finder_Getter,"METS",Selector_Base,Selector_Key,false);

Selector_Base *Jet_Finder_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<1) THROW(critical_error,"Invalid syntax");
  Jet_Finder *jf(new Jet_Finder(key.p_proc->NIn(),key.p_proc->NOut(),
				(Flavour*)&key.p_proc->Process()->
				Flavours().front(),key[0][0]));
  jf->SetProcess(key.p_proc);
  static bool menlots(false);
  if (!menlots && key.p_proc->Process()->Info().Has(nlo_type::real)) {
    menlots=true;
    rpa.gen.AddCitation(1,"NLO matrix element merging with truncated showers is "+
			std::string("published under \\cite{Hoeche:2010kg}."));
  }
  if (key.front().size()>1 && key[0][1]=="LO") jf->SetOn(false);
  return jf;
}

void Jet_Finder_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"METS jet finder"; 
}

}
