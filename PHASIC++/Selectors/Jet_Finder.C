#include "PHASIC++/Selectors/Jet_Finder.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PHASIC;
using namespace ATOOLS;

Jet_Finder::Jet_Finder
(const int nin,const int nout,Flavour *fl,
 const std::string &ycut,const std::string &gycut):
  Selector_Base("Jetfinder"), m_value(0.0), 
  m_cuttag(ycut), m_gcuttag(gycut), m_on(true)
{
  m_gycut=m_ycut=2.0;
  /*
  // something better needs to be done, only useful if single ycut, gycut
  if (m_cuttag[0]!='[') {
    Algebra_Interpreter interpreter;
    interpreter.AddTag("E_CMS",ToString(rpa.gen.Ecms()));
    m_ycut=ToType<double>(interpreter.Interprete(ycut));
    m_gycut=ToType<double>(interpreter.Interprete(gycut));
  }
  */
  m_fl=fl;
  m_nin=nin;
  m_nout=nout;
  m_n=m_nin+m_nout;
  /*
  for (size_t i(0);i<m_n;++i)
    if(m_fl[i].Strong()) {
      m_strongflavs.push_back(m_fl[i]);
      m_stronglocs.push_back(i);
    }
  */
  m_smax=m_s=sqr(rpa.gen.Ecms());
  // m_smin=m_ycut*m_s;
  m_pt2min=1.0;
  m_sel_log = new Selector_Log(m_name);
  rpa.gen.AddCitation(1,"Matrix element merging with truncated showers\
 is published under \\cite{Hoeche:2009rj}.");
}

Jet_Finder::~Jet_Finder() 
{
}

Flavour Jet_Finder::GetFlavour(std::string fl)
{
  bool bar(false);
  if (fl=="j") return Flavour(kf_jet);
  if (fl=="Q") return Flavour(kf_quark);
  if (fl=="G") return Flavour(kf_gluon);
  if (fl=="P") return Flavour(kf_photon);
  if (fl.length()>1) {
    if (fl[fl.length()-1]=='b') {
      fl.erase(fl.length()-1,1);
      bar=true;
    }
    else if ((fl[0]=='W' || fl[0]=='H')) {
      if (fl[fl.length()-1]=='-') {
	fl[fl.length()-1]='+';
	bar=true;
      }
    }
    else if (fl[fl.length()-1]=='+') {
      fl[fl.length()-1]='-';
      bar=true;
    }
  }
  if (fl=="Q") return Flavour(kf_quark); // why again?
  Flavour flav(s_kftable.KFFromIDName(fl));
  if (flav.Kfcode()==kf_none) 
    THROW(critical_error,"No flavour for '"+fl+"'.");
  if (bar) flav=flav.Bar();
  return flav;
}

size_t Jet_Finder::FillCombinations(const std::string &name,
				    const std::string &ycut,
				    const std::string &gycut,
				    size_t &cp,const int fl)
{
  bool ex(false);
  size_t sum(0), sp(0);
  std::vector<int> pos;
  std::string cut(ycut), ccut(cut), ncut(cut);
  std::string gcut(gycut), cgcut(gcut), ngcut(gcut);
  for (size_t i(0);i<name.length();++i) {
    if (name[i]=='[') {
      int open(1);
      for (size_t j(i+1);j<name.length();++j) {
	if (name[j]=='[') ++open;
	if (name[j]==']') --open;
	if (open==0) {
	  for (size_t ci(0);ci<cut.length();++ci) {
	    if (cut[ci]=='[') {
	      int copen(1);
	      for (size_t cj(ci+1);cj<cut.length();++cj) {
		if (cut[cj]=='[') ++copen;
		if (cut[cj]==']') --copen;
		if (copen==0) {
		  if (ccut==ycut) ccut=cut.substr(0,ci);
		  ncut=cut.substr(ci+1,cj-ci-1);
		  cut=cut.substr(cj+1);
		}
	      }
	    }
	  }
	  for (size_t ci(0);ci<gcut.length();++ci) {
	    if (gcut[ci]=='[') {
	      int copen(1);
	      for (size_t cj(ci+1);cj<gcut.length();++cj) {
		if (gcut[cj]=='[') ++copen;
		if (gcut[cj]==']') --copen;
		if (copen==0) {
		  if (cgcut==gycut) cgcut=gcut.substr(0,ci);
		  ngcut=gcut.substr(ci+1,cj-ci-1);
		  gcut=gcut.substr(cj+1);
		}
	      }
	    }
	  }
	  pos.push_back(FillCombinations
			(name.substr(i+1,j-i-1),ncut,ngcut,cp,fl-1));
	  m_flavs[pos.back()]=GetFlavour(name.substr(sp,i-sp));
	  if (pos.back()&((1<<m_nin)-1)) 
	    m_flavs[pos.back()]=m_flavs[pos.back()].Bar();
	  sum=sum|pos.back();
	  sp=i=j+3;
	  break;
	}
      }
    }
    else if (name[i]=='_') {
      if (name[i-1]!='_' && name[i+1]=='_') {
	pos.push_back(1<<cp++);
	m_flavs[pos.back()]=GetFlavour(name.substr(sp,i-sp));
	if (pos.back()&((1<<m_nin)-1)) 
	  m_flavs[pos.back()]=m_flavs[pos.back()].Bar();
	sum=sum|pos.back();
	ex=true;
      }
      if (ex && name[i+1]!='_') {
	sp=i+1;
	ex=false;
      }
    }
  }
  if (name[name.length()-1]!=']') {
    pos.push_back(1<<cp++);
    m_flavs[pos.back()]=GetFlavour(name.substr(sp,name.length()-sp));
    if (pos.back()&((1<<m_nin)-1)) 
      m_flavs[pos.back()]=m_flavs[pos.back()].Bar();
    sum=sum|pos.back();
  }
  Algebra_Interpreter interpreter;
  interpreter.AddTag("E_CMS",ToString(rpa.gen.Ecms()));
  m_cycut=ToType<double>(interpreter.Interprete(ccut));
  m_gcycut=ToType<double>(interpreter.Interprete(cgcut));
  int sc[3];
  for (size_t i(0);i<pos.size();++i) {
    m_ycuts[pos[i]][pos[i]]=m_cycut;
    m_gycuts[pos[i]][pos[i]]=m_gcycut;
    sc[0]=m_flavs[pos[i]].StrongCharge();
    if ((pos[i]&3)==0 && sc[0]!=0 &&
	(m_fl[0].Strong() || m_fl[1].Strong())) {
      bool found(false);
      for (size_t l(0);l<m_pcs.size();++l)
	if (m_pcs[l]==(size_t)pos[i]) {
	  found=true;
	  break;
	}
      if (!found) m_pcs.push_back(pos[i]);
    }
    for (size_t j(i+1);j<pos.size();++j) {
      if (pos[i]>2 || pos[j]>2) {
	m_ycuts[pos[i]][pos[j]]=m_cycut;
	m_gycuts[pos[i]][pos[j]]=m_gcycut;
	m_ycut=Min(m_ycut,m_cycut);
	m_gycut=Min(m_gycut,m_gcycut);
	sc[1]=m_flavs[pos[j]].StrongCharge();
	for (size_t k(0);k<pos.size();++k)
	  if (i!=k && j!=k) {
	    sc[2]=m_flavs[pos[k]].StrongCharge();
	    if (sc[0] && sc[1] && sc[2] &&
		(sc[0]==-sc[1] || abs(sc[0])!=3 || abs(sc[1])!=3) &&
		(sc[2]==-sc[1] || abs(sc[2])!=3 || abs(sc[1])!=3)) {
	      if (p_sproc && p_sproc->Combinable(pos[i],pos[j]))
		m_fills[fl].push_back(Comb_Key(pos[i],pos[j],pos[k]));
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
  if (m_ycuts.empty()) {
    if (p_proc==NULL) THROW(fatal_error,"Process not set.");
    m_procname=p_proc->Process()->Name();
    p_sproc=p_proc->Process()->Get<Single_Process>();
    if (m_procname.find("__QCD")!=std::string::npos)
      m_procname=m_procname.substr(0,m_procname.find("__QCD"));
    if (m_procname.find("__EW")!=std::string::npos)
      m_procname=m_procname.substr(0,m_procname.find("__EW"));
    m_moms.clear();
    m_flavs.clear();
    m_ycuts.clear();
    m_gycuts.clear();
    m_fills.resize(m_nin+m_nout+1);
    std::string name(m_procname.substr(m_procname.find('_')+1));
    name=name.substr(name.find("__")+2);
    size_t i(0);
    FillCombinations(name,m_cuttag,m_gcuttag,i,m_nin+m_nout);
    if (msg_LevelIsDebugging()) {
      msg_Out()<<METHOD<<"(): Combinations for '"<<m_procname<<"' {\n";
      double s(sqr(rpa.gen.Ecms()));
      for (std::map<size_t,std::map<size_t,double> >::const_iterator
	     iit(m_ycuts.begin());iit!=m_ycuts.end();++iit) {
	size_t i(iit->first);
	if (iit->second.size()>1)
	  msg_Out()<<"  "<<ID(i)<<"["<<m_flavs[i]<<","
		   <<m_flavs[i].Strong()<<"] & {";
	for (std::map<size_t,double>::const_iterator
	       jit(iit->second.begin());jit!=iit->second.end();++jit) {
	  size_t j(jit->first);
	  if (i!=j) 
	    msg_Out()<<" "<<ID(j)<<"["<<m_flavs[j]<<","
		     <<m_flavs[j].Strong()<<",("<<sqrt(m_ycuts[i][j]*s)
		     <<","<<sqrt(m_gycuts[i][j]*s)<<")]";
	}
	if (iit->second.size()>1) msg_Out()<<" }\n";
      }
      msg_Out()<<"}\n";
      msg_Out()<<METHOD<<"(): Identified clusterings {\n";
      for (size_t i(0);i<m_fills.size();++i)
	for (size_t j(0);j<m_fills[i].size();++j)
	  msg_Out()<<"  ["<<ID(m_fills[i][j].first)<<","
		   <<ID(m_fills[i][j].second)<<"] <-> "
		   <<ID(m_fills[i][j].partner)<<" ("<<i<<")\n";
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
  m_sok=false;
  for (int i(0);i<m_nin+m_nout;++i)
    if (m_flavs[1<<i].Strong() && 
	m_flavs[1<<i].StrongCharge()!=8) {
      m_sok=true;
      break;
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
  if (!m_on) return true;
  if (p_proc->Process()->IsGroup()) {
    bool trigger(false);
    for (size_t i(0);i<p_proc->Process()->Size();++i)
      if (!(*p_proc->Process())[i]->IsMapped())
	if ((*p_proc->Process())[i]->
	    Trigger(p)) trigger=true;
    return trigger;
  }
  return SingleTrigger(p);
}

bool Jet_Finder::SingleTrigger(const Vec4D_Vector &p)
{
  FillCombinations();
  PrepareMomList(p);
  bool uc(false);
  SP(Color_Integrator) ci(p_proc->ColorIntegrator());
  if (ci!=NULL && ci->On()) {
    uc=true;
    std::vector<int> ic(ci->I()), jc(ci->J());
    if (!PrepareColList(ic,jc)) return 1-m_sel_log->Hit(true);
    p_sproc=p_proc->Process()->Get<Single_Process>();
  }
  m_value=2.0;
  msg_Debugging()<<METHOD<<"(): '"
		 <<p_proc->Process()->Name()<<"' {\n";
  for (size_t i(0);i<m_pcs.size();++i) {
    msg_Debugging()<<"  "<<ID(m_pcs[i])<<"["<<m_flavs[m_pcs[i]]
		   <<"] -> "<<m_moms[m_pcs[i]].PPerp()<<"\n";
    if (m_moms[m_pcs[i]].PPerp2()<m_pt2min) return 1-m_sel_log->Hit(true);
  }
  for (size_t cl(1);cl<m_fills.size();++cl) {
    if (m_fills[cl].empty()) continue;
    msg_Indent();
    msg_Debugging()<<"level = "<<m_fills.size()-cl<<" {\n";
    for (size_t ps(0);ps<m_fills[cl].size();++ps) {
      size_t i(m_fills[cl][ps].first),
	j(m_fills[cl][ps].second), k(m_fills[cl][ps].partner);
      if (uc && !ColorConnected(i,j,k)) continue;
      if (m_flavs[i].IsQuark() && m_flavs[j].IsQuark() &&
	  m_flavs[i]!=m_flavs[j].Bar()) continue;
      double ycut(m_ycuts[i][j]);
      msg_Debugging()<<"  "<<ID(i)<<"["<<m_flavs[i]<<"] & "
		     <<ID(j)<<"["<<m_flavs[j]<<"] <-> "
		     <<ID(k)<<"["<<m_flavs[k]<<"], qcut = "
		     <<sqrt(ycut*m_s)<<"/"<<sqrt(m_gycuts[i][j]*m_s);
      double pt2ij=Qij2(m_moms[i],m_moms[j],m_moms[k],m_flavs[i],m_flavs[j]);
      msg_Debugging()<<", ptjk = "<<sqrt(pt2ij)<<" ("
		     <<(pt2ij>=ycut*m_s)<<(pt2ij<ycut*m_s?")\n":")");
      if (pt2ij<ycut*m_s) return 1-m_sel_log->Hit(true);
      if (pt2ij<m_value*m_s) m_value=pt2ij/m_s;
      msg_Debugging()<<"\n";
    }
    msg_Debugging()<<"}\n";
  }
  msg_Debugging()<<"} -> q_min = "<<sqrt(m_value*m_s)<<"\n";
  return 1-m_sel_log->Hit(false);
}

bool Jet_Finder::JetTrigger(const ATOOLS::Vec4D_Vector &p,
                            const ATOOLS::Flavour_Vector &fl, int n)
{
  // needs to be suitably defined, only useful for n = #coloured particles
  ATOOLS::Vec4D_Vector   colp;
  ATOOLS::Flavour_Vector colfl;
  for (size_t i(0);i<fl.size();++i)
    if (fl[i].Strong()) {
      colfl.push_back(fl[i]);
      colp.push_back(p[i]);
    }
  for (size_t i(0);i<colfl.size();++i)
    for (size_t j(0);j<i;++j)
      for (size_t k(0);k<colfl.size();++k)
        if (k!=i && k!=j)
          if (Qij2(colp[i],colp[j],colp[k],colfl[i],colfl[j])<m_smin)
            return false;
  return true;
}

bool Jet_Finder::JetTrigger(const ATOOLS::Vec4D_Vector &p)
{
  return Trigger(p);
  // only useful for n = #coloured particles
  for (size_t i(0);i<m_strongflavs.size();++i)
    for (size_t j(0);j<i;++j)
      for (size_t k(0);k<m_strongflavs.size();++k)
        if (k!=i && k!=j)
          if (Qij2(p[m_stronglocs[i]],p[m_stronglocs[j]],p[m_stronglocs[k]],m_strongflavs[i],m_strongflavs[j])>m_smin)
            return false;
  return true;
}

bool Jet_Finder::NoJetTrigger(const ATOOLS::Vec4D_Vector &p)
{
  // copied from NJet_Finder
  double s=(p[0]+p[1]).Abs2();
  return (s>m_smin*4.);
}

double Jet_Finder::Qij2(const Vec4D &pi,const Vec4D &pj,const Vec4D &pk,
			const Flavour &fi,const Flavour &fj,const int mode)
{
  Vec4D npi(pi), npj(pj);
  Flavour nfi(fi), nfj(fj);
  if (npi[0]<0.0) {
    npi=-pi-pj;
    if (mode==0) nfi=fi==fj.Bar()?Flavour(kf_gluon):(fi.IsGluon()?fj.Bar():fi);
  }
  else if (npj[0]<0.0) {
    npj=-pj-pi;
    if (mode==0) nfj=fj==fi.Bar()?Flavour(kf_gluon):(fj.IsGluon()?fi.Bar():fj);
  }
  if (nfi.IsQuark() && nfj.IsQuark() && nfi!=nfj.Bar()) return -1.0;
  double pipj(dabs(npi*npj)), pipk(dabs(npi*pk)), pjpk(dabs(npj*pk));
  double mti(sqr(Flavour(nfi).Mass())), mtj(sqr(Flavour(nfj).Mass()));
  if (pipj==0.0) {
    if (mti!=0.0 || mtj!=0.0) THROW(fatal_error,"Ill-defined mass term");
  }
  else {
    mti/=2.0*pipj;
    mtj/=2.0*pipj;
  }
  double Cij(nfj.IsGluon()?Max(0.0,pipk/(pipj+pjpk)-mti):1.0);
  double Cji(nfi.IsGluon()?Max(0.0,pjpk/(pipj+pipk)-mtj):1.0);
  return 4.0*dabs(pi*pj)/(Cij+Cji);
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

void Jet_Finder::UpdateCuts(double sprime,double y,Cut_Data *cuts) 
{
  if (!m_on) return;
  msg_Debugging()<<METHOD<<"(): {\n";
  for (int i(m_nin); i<m_nin+m_nout; ++i) {
    cuts->energymin[i] = m_fl[i].Mass();
    if (m_fl[i].Strong())
      for (int j(i+1); j<m_nin+m_nout; ++j) {
	if (m_fl[j].Strong())
	  if (GetYcut(1<<i,1<<j)>0.0) {
	    double scut=GetYcut(1<<i,1<<j)*m_s+
	      sqr(m_fl[i].Mass())+sqr(m_fl[j].Mass());
	    if (m_fl[i].StrongCharge()==8 || m_fl[j].StrongCharge()==8)
	      scut=Min(scut,Max(1.0,sqr(m_fl[i].Mass())+sqr(m_fl[j].Mass())));
	    msg_Debugging()<<"  ("<<i<<","<<j<<") -> "<<sqrt(scut)<<"\n";
	    cuts->scut[i][j]=cuts->scut[j][i]=Max(cuts->scut[i][j],scut);
	  }
      }
  }
  msg_Debugging()<<"}\n";
}

void Jet_Finder::BuildCuts(Cut_Data *cuts) 
{
  FillCombinations();
  UpdateCuts(0.0,0.0,cuts);
}

double Jet_Finder::ActualValue() const 
{
  return m_value; 
}

double Jet_Finder::GetYcut(const size_t& i,const size_t& j) const
{
  std::map<size_t,std::map<size_t,double> >::const_iterator it=m_ycuts.find(i);
  if(it==m_ycuts.end()) return -1.0;
  std::map<size_t,double>::const_iterator jt=(it->second).find(j);
  if(jt==(it->second).end()) return -1.0;
  return jt->second;
}

double Jet_Finder::GetGlobalYcut(const size_t &i,const size_t &j) const
{
  std::map<size_t,std::map<size_t,double> >::const_iterator it=m_gycuts.find(i);
  if(it==m_gycuts.end()) return -1.0;
  std::map<size_t,double>::const_iterator jt=(it->second).find(j);
  if(jt==(it->second).end()) return -1.0;
  return jt->second;
}

namespace PHASIC{

DECLARE_ND_GETTER(Jet_Finder_Getter,"METS",Selector_Base,Selector_Key,false);

Selector_Base *Jet_Finder_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<2) THROW(critical_error,"Invalid syntax");
  Jet_Finder *jf(new Jet_Finder(key.p_proc->NIn(),key.p_proc->NOut(),
				(Flavour*)&key.p_proc->Process()->
				Flavours().front(),key[0][0],key[0][1]));
  jf->SetProcess(key.p_proc);
  if (key.front().size()>2 && key[0][2]=="LO") jf->SetOn(false);
  return jf;
}

void Jet_Finder_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"METS jet finder"; 
}

}
