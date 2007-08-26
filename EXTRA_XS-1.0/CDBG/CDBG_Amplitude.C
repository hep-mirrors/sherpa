#include "CDBG_Amplitude.H"

#include "Message.H"
#include "Exception.H"
#include "MyStrStream.H"
#include "STL_Tools.H"
#include "Shell_Tools.H"

using namespace EXTRAXS;
using namespace ATOOLS;

// #define DEBUG__BG

static const double invsqrttwo(1.0/sqrt(2.0));

CDBG_Amplitude::CDBG_Amplitude():
  p_model(NULL), m_n(0), m_nf(6) {}

CDBG_Amplitude::~CDBG_Amplitude()
{
  CleanUp();
}

void CDBG_Amplitude::CleanUp()
{
  for (size_t i(0);i<m_cur.size();++i) 
    for (size_t j(0);j<m_cur[i].size();++j) delete m_cur[i][j]; 
  m_n=m_maxid=0;
  m_fl=Flavour_Vector();
  m_p=Vec4D_Vector();
  m_ch=Int_Vector();
  m_cl=Int_Matrix();
  m_cur=Current_Matrix();
  m_chirs=Int_Matrix();
  m_ress=Complex_Vector();
  m_hmap=Int_Vector();
  m_hamps=Int_Vector();
}

bool CDBG_Amplitude::MatchIndices(const Int_Vector &ids,const size_t &n,
				  const size_t &i,const size_t &j,
				  const size_t &k)
{
  for (size_t l(0);l<n;++l) {
    bool found(false), twice(false);
    for (size_t m(0);m<m_cur[i][j]->Id().size();++m) 
      if (m_cur[i][j]->Id()[m]==ids[l]) {
	if (found) twice=true;
	found=true;
      }
    for (size_t m(0);m<m_cur[n-i][k]->Id().size();++m) 
      if (m_cur[n-i][k]->Id()[m]==ids[l]) {
	if (found) twice=true;
	found=true;
      }
    if (!found || twice) return false;
  }
  return true;
}

void CDBG_Amplitude::AddCurrent(const Int_Vector &ids,const size_t &n,
				const Flavour &fl)
{
  // add new currents
  Current_Vector curs;
  String_Vector models(p_model->IncludedModels());
  for (size_t i(0);i<models.size();++i) {
    Current_Key key(fl,models[i],p_model);
    Current_Base *cur(Current_Getter::GetObject(key.Type(),key));
    if (cur!=NULL) curs.push_back(cur);
  }
  if (curs.empty()) THROW(fatal_error,"Invalid flavour");
  for (Current_Vector::const_iterator cit(curs.begin());
       cit!=curs.end();++cit) {
    std::set<Vertex_Key> v3;
    // compose current from all possible subcurrents
    for (size_t i(1);i<n;++i) {
      for (size_t j(0);j<m_cur[i].size();++j) {
	for (size_t k(0);k<m_cur[n-i].size();++k) {
	  if (!MatchIndices(ids,n,i,j,k)) continue;
	  for (size_t m(0);m<models.size();++m) {
	    Vertex_Key key(m_cur[i][j],m_cur[n-i][k],*cit,models[m],p_model);
	    if (v3.find(key.SwapAB())!=v3.end()) continue;
	    Vertex *v(Vertex_Getter::GetObject(key.Type(),key));
	    if (v!=NULL) {
	      v->SetJA(key.p_a);
	      v->SetJB(key.p_b);
	      v3.insert(key);
	    }
	  }
	}
      }
    }
    if (!v3.empty() || n==1) {
      Int_Vector isfs(ids.size());
      for (size_t i(0);i<ids.size();++i)
	isfs[i]=m_fl[ids[i]].IsFermion();
      (*cit)->SetId(ids);
      (*cit)->SetFId(isfs);
      (*cit)->FindPermutations();
      (*cit)->SetKey(m_cur[n].size());
      m_cur[n].push_back(*cit);
      (*cit)->Print();
    }
    else delete *cit;
  }
}

void CDBG_Amplitude::Construct(Int_Vector ids,const size_t &n)
{
  if (ids.size()==n) {
    if (n==1) {
      AddCurrent(ids,n,m_fl[ids.back()]);
    }
    else {
      for (size_t i(1);i<=m_nf;++i) {
	AddCurrent(ids,n,Flavour((kf::code)i));
	AddCurrent(ids,n,Flavour((kf::code)i).Bar());
      }
      AddCurrent(ids,n,kf::gluon);
      // add pseudoparticles for 4-gluon vertex
      if (n>1 && n<m_n-1) AddCurrent(ids,n,kf::gluonqgc);
    }
    return;
  }
  // currents are unordered -> use ordered indexing
  size_t last(ids.empty()?-1:ids.back());
  ids.push_back(0);
  if (n==m_n-1) {
    // calculate only one final current
    ids.back()=last+1;
    Construct(ids,n);
    return;
  }
  for (size_t m(n>1?m_n-1:m_n), i(last+1);i<m;++i) {
    // fill currents 0..n-1 for external partons
    // currents 0..n-2 each index for internal partons
    ids.back()=i;
    Construct(ids,n);
  }
}

bool CDBG_Amplitude::Construct(const Flavour_Vector &flavs)
{
  if (msg_LevelIsDebugging()) {
    msg_Out()<<METHOD<<"(): {\n\n   Implemented currents:\n\n";
    Current_Getter::PrintGetterInfo(msg_Out(),15);
    msg_Out()<<"\n   Implemented vertices:\n\n";
    Vertex_Getter::PrintGetterInfo(msg_Out(),15);
    msg_Out()<<"\n}\n";
  }
  m_fl=flavs;
  m_n=m_fl.size();
  m_p.resize(m_n);
  m_ch.resize(m_n);
  m_cl.resize(m_n,Int_Vector(2));
  m_cur.resize(m_n);
  Int_Vector ids;
  size_t tsum(0);
  for (size_t i(1);i<m_n;++i) {
    Construct(ids,i);
    tsum+=m_cur[i].size();
    msg_Debugging()<<"P("<<i<<","<<m_n<<"): "<<m_cur[i].size()<<"\n";
  }
  msg_Debugging()<<"total: cur = "<<tsum<<"\n";
  return true;
}

void CDBG_Amplitude::CalcJL()
{
  for (size_t i(0);i<m_cur[1].size();++i) 
    m_cur[1][i]->ConstructJ(m_p[i],m_ch[i],m_cl[i][0],m_cl[i][1]);
  for (size_t n(2);n<m_n;++n) 
    for (size_t i(0);i<m_cur[n].size();++i) 
      m_cur[n][i]->Evaluate();
}

void CDBG_Amplitude::SetMomenta(const Vec4D_Vector &moms)
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"():\n";
#endif
  Vec4D sum;
  for (size_t i(0);i<m_n;++i) {
    m_p[i]=moms[i];
    sum+=m_p[i];
#ifdef DEBUG__BG
    msg_Debugging()<<"set p["<<i<<"] = "<<m_p[i]<<"\n";
#endif
  }
  static double accu(sqrt(Accu()));
  CVec4D::SetAccu(accu);
  if (!((CVec4D)sum).IsZero()) 
    msg_Error()<<METHOD<<"(): Four momentum not conserved. sum = "
	       <<sum<<"."<<std::endl;
  CVec4D::SetAccu(Accu());
}

void CDBG_Amplitude::SetColors(const Int_Vector &rc,
			       const Int_Vector &ac)
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"():\n";
#endif
  for (size_t i(0);i<m_n;++i) {
    if (m_fl[i].IsGluon()) {
      m_cl[i][0]=rc[i];
      m_cl[i][1]=ac[i];
    }
    else {
      if (m_fl[i].IsAnti()) {
	m_cl[i][1]=ac[i];
	m_cl[i][0]=0;
      }
      else {
	m_cl[i][0]=rc[i];
	m_cl[i][1]=0;
      }
    }
#ifdef DEBUG__BG
    msg_Debugging()<<"m_cl["<<i<<"][0] = "<<m_cl[i][0]
		   <<", m_cl["<<i<<"][1] = "<<m_cl[i][1]<<"\n";
#endif
  }
}

Complex CDBG_Amplitude::Evaluate(const Int_Vector &chirs)
{
  size_t id(MakeId(chirs));
  if (id>m_maxid || m_hamps[id]==0) {
    msg_Debugging()<<METHOD<<"(): Configuration "
		   <<chirs<<" does not exist."<<id<<std::endl;
    return 0.0;
  }
  for (size_t j(0);j<m_n;++j) m_ch[j]=chirs[j];
  CalcJL();
  size_t ihp(MakeId(m_ch)), ihm((1<<m_n)-1-ihp);
  Complex res(m_cur[1].back()->Contract(*m_cur.back()[0],ihm,ihp));
  msg_Debugging()<<"A"<<chirs<<" = "<<res<<" "
		 <<std::abs(res)<<" {"<<ihm<<","<<ihp<<"}\n";
  return res;
}

bool CDBG_Amplitude::EvaluateAll()
{
  for (size_t j(0);j<m_n;++j) m_ch[j]=0;
  CalcJL();
  for (size_t i(0);i<m_ress.size();++i) {
    if (m_hmap[i]>=0) {
      m_ress[i]=std::conj(m_ress[m_hmap[i]]);
      continue;
    }
    size_t ihp(MakeId(m_chirs[i])), ihm((1<<m_n)-1-ihp);
    m_ress[i]=m_cur[1].back()->Contract(*m_cur.back()[0],ihm,ihp);
#ifdef DEBUG__BG
    msg_Debugging()<<"A"<<m_chirs[i]<<" = "<<m_ress[i]<<" "
		   <<std::abs(m_ress[i])<<" {"<<ihm<<","<<ihp<<"}\n";
#endif
  }
  return true;
}

bool CDBG_Amplitude::CheckChirs(const Int_Vector &chirs)
{
  size_t p(0), m(0);
  Int_Vector q(m_nf+1,0);
  for (size_t i(0);i<chirs.size();++i) {
    if (m_fl[i].IsQuark()) q[m_fl[i].Kfcode()]+=chirs[i];
    if (chirs[i]>0) ++p;
    else if (chirs[i]<0) ++m;
    else THROW(fatal_error,"Invalid helicities");
  }
  for (size_t i(0);i<q.size();++i) 
    if (q[i]!=0) return false;
  return p>1 && m>1;
}

bool CDBG_Amplitude::MapChirs(const Int_Vector &chirs)
{
  Int_Vector rchirs(chirs.size());
  for (size_t i(0);i<chirs.size();++i) rchirs[i]=-chirs[i];
  m_hmap.push_back(-1);
  for (size_t i(0);i<m_fl.size();++i) 
    if (!m_fl[i].IsGluon()) return false;
  for (size_t i(0);i<m_chirs.size();++i) {
    bool hit(true);
    for (size_t j(0);j<rchirs.size();++j)
      if (m_chirs[i][j]!=rchirs[j]) hit=false;
    if (!hit) continue;
    m_hmap.back()=i;
    msg_Debugging()<<"mapped helicities "<<(m_hmap.size()-1)
		   <<" "<<chirs<<" -> "<<m_chirs[i]<<"\n";
    return true;
  }
  return false;
}

bool CDBG_Amplitude::ConstructChirs(Int_Vector chirs,const size_t &i)
{
  if (i==chirs.size()) {
    if (CheckChirs(chirs)) {
#ifdef DEBUG__BG
      msg_Debugging()<<METHOD<<"(): Add configuration "<<chirs<<"\n";
#endif
      m_chirs.push_back(chirs);
      m_maxid=Max(m_maxid,MakeId(chirs));
      MapChirs(chirs);
      m_ress.push_back(0.0);
    }
    return true;
  }
  if (chirs[i]!=0) {
    ConstructChirs(chirs,i+1);
  }
  else {
    for (int ch(1);ch>=-1;ch-=2) {
      chirs[i]=ch;
      ConstructChirs(chirs,i+1);
    }
  }
  return true;
}

bool CDBG_Amplitude::Construct(const Flavour_Vector &flavs,
			       XS_Model_Base *const model)
{
  CleanUp();
  p_model=model;
  if (!Construct(flavs)) return false;
  Int_Vector chirs(flavs.size(),0);
  if (!ConstructChirs(chirs,0)) return false;
  m_hamps.resize(m_maxid+1,0);
  for (size_t i(0);i<m_chirs.size();++i)
    m_hamps[MakeId(m_chirs[i])]=1;
  return true;
}

void CDBG_Amplitude::SetGauge(const size_t &n)
{
  Vec4D k(1.0,0.0,1.0,0.0);
  switch(n) {
  case 1: k=Vec4D(1.0,1.0,0.0,0.0); break;
  case 2: k=Vec4D(1.0,invsqrttwo,invsqrttwo,0.0); break;
  case 3: k=Vec4D(1.0,invsqrttwo,-invsqrttwo,0.0); break;
  case 4: k=Vec4D(1.0,invsqrttwo,0.0,invsqrttwo); break;
  case 5: k=Vec4D(1.0,invsqrttwo,0.0,-invsqrttwo); break;
  case 6: k=Vec4D(1.0,0.0,invsqrttwo,invsqrttwo); break;
  case 7: k=Vec4D(1.0,0.0,invsqrttwo,-invsqrttwo); break;
  }
  for (size_t i(0);i<m_cur[1].size();++i) m_cur[1][i]->SetGauge(k);
}

bool CDBG_Amplitude::GaugeTest(const Vec4D_Vector &moms)
{
  msg_Info()<<METHOD<<"(): Performing gauge test ..."<<std::flush;
  msg_Indent();
  SetGauge(0);
  SetMomenta(moms);
  if (!EvaluateAll()) return false;
  Complex_Vector ress(m_ress);
  SetGauge(1);
  SetMomenta(moms);
  if (!EvaluateAll()) return false;
  double mean(0.0);
  for (size_t i(0);i<m_ress.size();++i) mean+=std::abs(ress[i]);
  mean/=m_ress.size();
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_ress.size();++i) {
    msg_Debugging()<<"A("<<m_chirs[i]
		   <<") = "<<m_ress[i]<<" vs. "<<ress[i]<<" -> dev. "
		   <<m_ress[i].real()/ress[i].real()-1.0<<" "
		   <<m_ress[i].imag()/ress[i].imag()-1.0<<"\n";
    double accu(sqrt(Accu()));
    if (!IsEqual(m_ress[i].real(),ress[i].real(),accu) ||
	!IsEqual(m_ress[i].imag(),ress[i].imag(),accu)) {
      double rat(mean/std::abs(m_ress[i])*Accu());
      if (IsEqual(m_ress[i].real(),ress[i].real(),rat) &&
	  IsEqual(m_ress[i].imag(),ress[i].imag(),rat)) {
	msg_Error().precision(12);
	msg_Error()<<METHOD<<"(): Large deviation: "
		   <<m_ress[i]<<" vs. "<<ress[i]<<"\n  => ("
		   <<(m_ress[i].real()/ress[i].real()-1.0)<<","
		   <<(m_ress[i].imag()/ress[i].imag()-1.0)
		   <<") {"<<rat<<"}."<<std::endl;
	msg_Error().precision(6);
      }
      else {
	msg_Error().precision(12);
	msg_Error()<<METHOD<<"(): Gauge test failed. "
		   <<m_ress[i]<<" vs. "<<ress[i]<<"\n  => ("
		   <<(m_ress[i].real()/ress[i].real()-1.0)<<","
		   <<(m_ress[i].imag()/ress[i].imag()-1.0)
		   <<") {"<<rat<<"}."<<std::endl;
	msg_Error().precision(6);
	return false;
      }
    }
  }
  msg_Debugging()<<"}\n";
  msg_Info()<<"satisfied."<<std::endl;
  return true;
}

void CDBG_Amplitude::WriteOutGraph
(std::ostream &str,Graph_Node *graph,size_t &ng) const
{
  if ((*graph)->empty()) {
    size_t fp(0);
    for (size_t j(1);j<graph->size();++j)
      if ((*graph)[j].find("%%")==std::string::npos) {
	std::string cl((*graph)[j]);
	size_t bpos(cl.find("F="));
	if (bpos!=std::string::npos) {
	  size_t epos(bpos+=2);
	  for (;cl[epos]>=48 && cl[epos]<=57;++epos);
	  fp+=ToType<size_t>(cl.substr(bpos,epos-bpos));
	}
      }
    str<<"  \\parbox{"<<(5*m_n+10)<<"mm}{Graph "<<++ng
       <<", $\\sum \\rm F$="<<fp<<" ("<<(fp%2==0?'+':'-')
       <<")\\begin{center}\n";
    str<<"  \\begin{fmfgraph*}("<<(10*m_n)<<","<<(10*m_n)<<")\n";
    str<<"    \\fmfsurround{"<<graph->front()<<"}\n";
    for (size_t j(0);j<m_n;++j) 
      str<<"    \\fmfv{decor.size=0ex,label=$J_{"
	 <<j<<"}$}{j_"<<(1<<j)<<"}\n";
    for (size_t j(1);j<graph->size();++j)
      if ((*graph)[j].find("%%")==std::string::npos) 
	str<<(*graph)[j]<<"\n";
    str<<"  \\end{fmfgraph*}\\end{center}\\vspace*{5mm}}";
    if (ng>0 && ng%3==0) str<<" \\\\\n\n";
    else str<<" &\n\n";
  }
  else {
    for (size_t i(0);i<(*graph)->size();++i)
      WriteOutGraph(str,(*graph)()[i],ng);
  }
}

void CDBG_Amplitude::WriteOutGraphs(const std::string &file) const
{
  Graph_Node graphs("j_"+ToString(1<<(m_n-1)),true);
  graphs.push_back("    %% "+graphs.back());
  m_cur.back().front()->CollectGraphs(&graphs);
  MakeDir(file,448);
  std::ofstream str((file+"/graphs.tex").c_str());
  str<<"\\documentclass[a4paper]{article}\n\n";
  str<<"\\usepackage{feynmp}\n";
  str<<"\\usepackage{amsmath}\n";
  str<<"\\usepackage{amssymb}\n";
  str<<"\\usepackage{longtable}\n";
  str<<"\\usepackage{pst-text}\n\n";
  str<<"\\setlength{\\headheight}{-1cm}\n";
  str<<"\\setlength{\\headsep}{0mm}\n";
  str<<"\\setlength{\\oddsidemargin}{-1in}\n";
  str<<"\\setlength{\\evensidemargin}{-1in}\n";
  str<<"\\setlength{\\textheight}{28truecm}\n";
  str<<"\\setlength{\\textwidth}{21truecm}\n\n";
  str<<"\\begin{document}\n";
  str<<"\\begin{fmffile}{graphs_fg}\n\n";
  str<<"  \\fmfset{thick}{1.25thin}\n";
  str<<"  \\fmfset{arrow_len}{2mm}\n";
  str<<"  \\fmfset{curly_len}{1.5mm}\n";
  str<<"  \\fmfset{wiggly_len}{1.5mm}\n";
  str<<"  \\unitlength=0.5mm\n\n";
  str<<"  \\pagestyle{empty}\n\n";
  str<<"  \\begin{longtable}{ccc}\n\n";
  size_t ng(0);
  WriteOutGraph(str,&graphs,ng);
  str<<"\n  \\end{longtable}\n\n";
  str<<"\\end{fmffile}\n";
  str<<"\\end{document}\n";
}
