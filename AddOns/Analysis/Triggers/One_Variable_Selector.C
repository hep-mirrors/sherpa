#include "AddOns/Analysis/Triggers/Trigger_Base.H"

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "ATOOLS/Math/Variable.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  typedef std::vector<int>        Int_Vector;
  typedef std::vector<Int_Vector> Int_Matrix;

  typedef std::vector<double>        Double_Vector;
  typedef std::vector<Double_Vector> Double_Matrix;

  typedef std::vector<ATOOLS::Flavour_Vector> Flavour_Matrix;

  typedef std::vector<ATOOLS::Variable_Base<double>*> Variable_Vector;

  typedef std::vector<ATOOLS::Histogram*> Histogram_Vector;

  class OVS_Tag_Replacer: public Tag_Replacer {
  private:
    Primitive_Analysis *p_ana;
  public:
    OVS_Tag_Replacer(Primitive_Analysis *const ana): p_ana(ana) {}
    std::string ReplaceTags(std::string &expr) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;
  };// end of class OVS_Tag_Replacer

  class One_Variable_Selector: public Trigger_Base {
  private:
    String_Matrix m_cndlist;
    Flavour_Matrix   m_flavs;
    Int_Matrix       m_items;
    String_Vector    m_vtags;
    Variable_Vector  m_vars;
    Double_Vector    m_mins, m_maxs;
    Double_Matrix    m_histos;
    Histogram_Vector m_dists;
    ATOOLS::Histogram *p_flow;
    OVS_Tag_Replacer m_repl;
  public:
    One_Variable_Selector
    (const std::string &inlist,const std::string &outlist,
     const String_Matrix &cndlist,const int isobs,
     const Flavour_Matrix &flavs,const Int_Matrix &items,
     const String_Vector &vtags,const Double_Vector &mins,
     const Double_Vector &maxs,const Double_Matrix &m_histos,
     Primitive_Analysis *const ana,const std::string &name="");
    ~One_Variable_Selector();
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  ATOOLS::Particle_List &outlist,
		  double weight, double ncount);
    int Evaluate(ATOOLS::Particle_List &reflist,
		 double weight,double ncount,ATOOLS::Particle_List &moms,
		 const size_t i,size_t j,size_t k,size_t o,size_t &eval);
    Analysis_Object &operator+=(const Analysis_Object &obj);
    void EndEvaluation(double scale=1.0);
    void Restore(double scale=1.0);
    void Output(const std::string & pname);
    Analysis_Object *GetCopy() const;    
  };// end of class One_Variable_Selector

} // namespace ANALYSIS

using namespace ANALYSIS;

DECLARE_GETTER(One_Variable_Selector_Getter,"MomSel",
 	       Analysis_Object,Argument_Matrix);

void One_Variable_Selector_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"CndList inlist outlist\n"
     <<std::setw(width+7)<<" "<<"Tags    flavi1,.. itemi1,.. "
     <<"vari mini maxi [mini maxi binsi typei [itemi]]\n"
     <<std::setw(width+7)<<" "<<"Flavs   flav11,.. .. flavN1,..\n"
     <<std::setw(width+7)<<" "<<"Items   item11,.. .. itemN1,..\n"
     <<std::setw(width+7)<<" "<<"Vars    var1      .. varN\n"
     <<std::setw(width+7)<<" "<<"Mins    min1      .. minN\n"
     <<std::setw(width+7)<<" "<<"Maxs    max1      .. maxN\n"
     <<std::setw(width+7)<<" "<<"HTypes  [type1   [.. typeN]]\n"
     <<std::setw(width+7)<<" "<<"HItems  [item1   [.. itemN]]\n"
     <<std::setw(width+7)<<" "<<"HBins   [bins1   [.. binsN]]\n"
     <<std::setw(width+7)<<" "<<"HMins   [min1    [.. minN]]\n"
     <<std::setw(width+7)<<" "<<"HMaxs   [max1    [.. maxN]]\n"
     <<std::setw(width+4)<<" "<<"}";
}

template <class __T>
std::string MakeString(const std::vector<__T> &v);

template <> std::string MakeString<ATOOLS::Flavour>
(const std::vector<ATOOLS::Flavour> &v)
{
  if (v.empty()) return "";
  std::string s((v.front().IsAnti()?"-":"")
		+ToString(v.front().Kfcode()));
  for (size_t i(1);i<v.size();++i)
    if (v[i].IsAnti()) s+=",-"+ToString(v[i].Kfcode());
    else s+=","+ToString(v[i].Kfcode());
  return s;
}

template <class __T>
std::string MakeString(const std::vector<__T> &v)
{
  if (v.empty()) return "";
  std::string s(ToString(v.front()));
  for (size_t i(1);i<v.size();++i) s+=","+ToString(v[i]);
  return s;
}

Analysis_Object *
One_Variable_Selector_Getter::operator()
  (const Argument_Matrix &parameters) const
{
  std::string inlist("FinalState"), outlist("Selected");
  int isobs(0);
  String_Matrix cndlist;
  Flavour_Matrix flavs;
  Int_Matrix items;
  String_Vector vtags;
  Double_Vector mins, maxs;
  Double_Matrix histos(5);
  Data_Reader reader(",",";","!","=");
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="RefList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
    else if (cur[0]=="CndList" && cur.size()>2) {
      cndlist.push_back(String_Vector(2));
      cndlist.back().front()=cur[1];
      cndlist.back().back()=cur[2];
    }
    else if (cur[0]=="IsObs" && cur.size()>1)
      isobs=ToType<int>(cur[1]);
    else if (cur[0]=="Tags" && cur.size()>5) {
      reader.SetString(cur[1]);
      std::vector<int> cfl;
      if (!reader.VectorFromString(cfl,"")) {
	msg_Debugging()<<METHOD<<"(): Invalid flavour syntax '"
		       <<cur[1]<<"'.\n";
	continue;
      }
      flavs.push_back(Flavour_Vector(cfl.size()));
      for (size_t k(0);k<cfl.size();++k) {
	flavs.back()[k]=Flavour((kf_code)abs(cfl[k]));
	if (cfl[k]<0) flavs.back()[k]=flavs.back()[k].Bar();
      }
      reader.SetString(cur[2]);
      std::vector<int> cit;
      if (!reader.VectorFromString(cit,"")) {
	msg_Debugging()<<METHOD<<"(): Invalid item syntax '"
		       <<cur[2]<<"'.\n";
	continue;
      }
      items.push_back(cit);
      vtags.push_back(cur[3]);
      mins.push_back(ToType<double>(cur[4]));
      maxs.push_back(ToType<double>(cur[5]));
      if (cur.size()<=9) {
	for (size_t i(0);i<histos.size();++i) histos[i].push_back(-1);
      }
      else {
	histos[0].push_back(HistogramType(cur[9]));
	histos[1].push_back(ToType<double>(cur[8]));
	histos[2].push_back(ToType<double>(cur[6]));
	histos[3].push_back(ToType<double>(cur[7]));
	if (cur.size()>10) histos[4].push_back(ToType<double>(cur[10]));
	else histos[4].push_back(-1);
      }
    }
    else if (cur[0]=="Flavs" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) {
	reader.SetString(cur[j]);
	std::vector<int> cfl;
	if (!reader.VectorFromString(cfl,"")) {
	  msg_Debugging()<<METHOD<<"(): Invalid flavour syntax '"
			 <<cur[j]<<"'.\n";
	  continue;
	}
	flavs.push_back(Flavour_Vector(cfl.size()));
	for (size_t k(0);k<cfl.size();++k) {
	  flavs.back()[k]=Flavour((kf_code)abs(cfl[k]));
	  if (cfl[k]<0) flavs.back()[k]=flavs.back()[k].Bar();
	}
      }
    }
    else if (cur[0]=="Items" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) {
	reader.SetString(cur[j]);
	std::vector<int> cit;
	if (!reader.VectorFromString(cit,"")) {
	  msg_Debugging()<<METHOD<<"(): Invalid item syntax '"
			 <<cur[j]<<"'.\n";
	  continue;
	}
	items.push_back(cit);
      }
    }
    else if (cur[0]=="Vars" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	vtags.push_back(cur[j]);
    }
    else if (cur[0]=="Mins" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	mins.push_back(ToType<double>(cur[j]));
    }
    else if (cur[0]=="Maxs" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	maxs.push_back(ToType<double>(cur[j]));
    }
    else if (cur[0]=="HTypes" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	histos[0].push_back(HistogramType(cur[j]));
    }
    else if (cur[0]=="HItems" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	histos[4].push_back(ToType<double>(cur[j]));
    }
    else if (cur[0]=="HBins" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	histos[1].push_back(ToType<double>(cur[j]));
    }
    else if (cur[0]=="HMins" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	histos[2].push_back(ToType<double>(cur[j]));
    }
    else if (cur[0]=="HMaxs" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	histos[3].push_back(ToType<double>(cur[j]));
    }
  }
  if (flavs.empty() || items.empty() || vtags.empty() || 
      mins.empty() || maxs.empty()) {
    msg_Debugging()<<METHOD<<"(): Cannot initialize selector.\n";
    return NULL;
  }
  if (histos[0].empty()) histos[0].push_back(-1.0);
  size_t max(Max(vtags.size(),Max(flavs.size(),items.size())));
  max=Max(max,Max(mins.size(),maxs.size()));
  flavs.resize(max,flavs.back());
  items.resize(max,items.back());
  vtags.resize(max,vtags.back());
  mins.resize(max,mins.back());
  maxs.resize(max,maxs.back());
  for (size_t i(0);i<5;++i) histos[i].resize(max,-1);
  for (size_t i(flavs.size());i<max;++i) {
    max=Max(flavs[i].size(),items[i].size());
    for (size_t j(flavs[i].size());j<max;++j) 
      flavs[i].push_back(flavs[i].back());
    for (size_t j(items[i].size());j<max;++j) 
      items[i].push_back(items[i].back());
  }
  msg_Debugging()<<METHOD<<"(): Initialized '"<<isobs<<"' {\n";
  for (size_t i(0);i<max;++i) {
    msg_Debugging()<<"    Tags "<<std::setw(12)<<MakeString(flavs[i])
		   <<" "<<std::setw(9)<<MakeString(items[i])
		   <<" "<<std::setw(9)<<vtags[i]
		   <<" "<<std::setw(9)<<mins[i]
		   <<" "<<std::setw(9)<<maxs[i]
		   <<" "<<std::setw(6)<<histos[2][i]
		   <<" "<<std::setw(6)<<histos[3][i]
		   <<" "<<std::setw(3)<<histos[1][i]
		   <<" "<<std::setw(2)<<histos[0][i]
		   <<" "<<std::setw(2)<<histos[4][i]<<"\n";
  }
  msg_Debugging()<<"}\n";
  return new One_Variable_Selector
    (inlist,outlist,cndlist,isobs,flavs,items,vtags,
     mins,maxs,histos,parameters());
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ATOOLS;

One_Variable_Selector::One_Variable_Selector
(const std::string &inlist,const std::string &outlist,
 const String_Matrix &cndlist,const int isobs,
 const Flavour_Matrix &flavs,const Int_Matrix &items,
 const String_Vector &vtags,const Double_Vector &mins,
 const Double_Vector &maxs,const Double_Matrix &histos,
 Primitive_Analysis *const ana,const std::string &name):
  Trigger_Base(inlist,outlist), m_cndlist(cndlist),
  m_flavs(flavs), m_items(items), m_vtags(vtags), m_mins(mins), m_maxs(maxs), 
  m_histos(histos), m_dists(flavs.size(),NULL), m_repl(ana)
{
  m_isobs=isobs;
  msg_Debugging()<<METHOD<<"(): {\n";
  m_vars.resize(m_vtags.size(),NULL);
  for (size_t i(0);i<m_vtags.size();++i) {
    if (m_vtags[i].find("Calc")==std::string::npos)
      m_vars[i]=ATOOLS::Variable_Getter::GetObject(m_vtags[i],m_vtags[i]);
    else m_vars[i]=ATOOLS::Variable_Getter::GetObject(m_vtags[i],"Calc(1)");
    if (m_vars[i]==NULL) THROW
      (fatal_error,"Variable '"+m_vtags[i]+"' does not exist. Run 'Sherpa"+
       std::string(" SHOW_Analysis_SYNTAX=1' to list variables."));
    ATOOLS::Algebra_Interpreter *inter=m_vars[i]->GetInterpreter();
    if (inter!=NULL) {
      const String_BlobDataBase_Map &data(ana->GetData());
      for (String_BlobDataBase_Map::const_iterator 
	     dit(data.begin());dit!=data.end();++dit) {
	Blob_Data<double> *dat(dynamic_cast<Blob_Data<double>*>(dit->second));
	if (dat!=NULL) inter->AddTag(dit->first,"0");
	else inter->AddTag(dit->first,"(0,0,0,0)");
      }
      m_vars[i]->Init(m_vtags[i]+"{"+ToString((long unsigned int)(&m_repl))+"}");
    }
  }
  if (name!="") m_name=name;
  else {
    size_t n(0);
    while (ana->GetObject(m_inlist+"_"+ToString(n))!=NULL) ++n;
    m_name=m_inlist+"_"+ToString(n);
  }
  p_flow = new ATOOLS::Histogram(1,0.0,(double)m_flavs.size(),m_flavs.size());
  for (size_t i(0);i<m_dists.size();++i)
    if (m_histos[0][i]>-1) {
      msg_Debugging()<<"  init histo "<<i<<" "<<m_vars[i]->Name()<<" for ";
      for (size_t j(0);j<m_flavs[i].size();++j) 
	msg_Debugging()<<m_flavs[i][j].IDName()<<" "<<m_items[i][j]<<" ";
      msg_Debugging()<<"-> type "<<m_histos[0][i]<<", "<<m_histos[1][i]
		     <<" bins, min "<<m_histos[2][i]<<", max "
		     <<m_histos[3][i]<<", item "<<m_histos[4][i]<<"\n";
      m_dists[i] = new ATOOLS::Histogram
	((int)m_histos[0][i],m_histos[2][i],
	 m_histos[3][i],(int)m_histos[1][i]);
      if (m_histos[4].size()<=i) m_histos[4].push_back(-1);
    }
  msg_Debugging()<<"}\n";
}

Analysis_Object &One_Variable_Selector::operator+=
(const Analysis_Object &obj)
{
  const One_Variable_Selector *vob((const One_Variable_Selector*)&obj);
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) *m_dists[i]+=*vob->m_dists[i];
  *p_flow+=*vob->p_flow;
  return *this;
}

void One_Variable_Selector::EndEvaluation(double scale) 
{
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) {
      m_dists[i]->Finalize();
      if (scale!=1.0) m_dists[i]->Scale(scale);
    }
  p_flow->Finalize();
  if (scale!=1.0) p_flow->Scale(scale);
}

void One_Variable_Selector::Restore(double scale) 
{
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) {
      if (scale!=1.0) m_dists[i]->Scale(scale);
      m_dists[i]->Restore();
    }
  if (scale!=1.0) p_flow->Scale(scale);
  p_flow->Restore();
}

void One_Variable_Selector::Output(const std::string & pname) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  std::string bname(pname+"/"+m_name);
  ATOOLS::MakeDir(pname); 
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) {
      std::string name(bname+"_"+m_vars[i]->IDName());
      for (size_t j(0);j<m_flavs[i].size();++j) 
	name+="_"+m_flavs[i][j].IDName()+ToString(m_items[i][j]);
      msg_Debugging()<<"  write '"<<name<<".dat'\n";
      m_dists[i]->Output((name+".dat").c_str());
    }
  if (m_isobs) p_flow->Output((bname+"_eff.dat").c_str());
  msg_Debugging()<<"}\n";
}

One_Variable_Selector::~One_Variable_Selector()
{
  while (m_vars.size()) {
    delete m_vars.back();
    m_vars.pop_back();
  }
  while (m_dists.size()) {
    if (m_dists.back()!=NULL) delete m_dists.back();
    m_dists.pop_back();
  }
  delete p_flow;
}

int One_Variable_Selector::Evaluate
(ATOOLS::Particle_List &reflist,double weight,double ncount,
 ATOOLS::Particle_List &moms,const size_t i,size_t j,size_t k,
 size_t o,size_t &eval) 
{
  bool count(m_vars[i]->IDName()=="Count");
  if (j>=m_flavs[i].size()) {
    ++eval;
    if (count) return 1;
    std::vector<ATOOLS::Vec4D> vmoms(moms.size());
    for (size_t l(0);l<vmoms.size();++l) vmoms[l]=moms[l]->Momentum();
    double val(m_vars[i]->Value(&vmoms.front(),vmoms.size()));
    bool pass(val>=m_mins[i] && val<=m_maxs[i]);
    if (msg_LevelIsDebugging()) {
      msg_Debugging()<<"  "<<m_flavs[i][0].IDName();
      for (size_t k(1);k<m_flavs[i].size();++k) 
	msg_Debugging()<<","<<m_flavs[i][k].IDName();
      msg_Debugging()<<" "<<m_items[i][0];
      for (size_t k(1);k<m_items[i].size();++k) 
	msg_Debugging()<<","<<m_items[i][k];
      msg_Debugging()<<" "<<m_vars[i]->Name()<<"("<<moms.front()->Number();
      for (size_t k(1);k<moms.size();++k) 
	msg_Debugging()<<","<<moms[k]->Number();
      msg_Debugging()<<") = "<<val<<" "<<(pass?"\\in":"\\nin")
		     <<" ["<<m_mins[i]<<","<<m_maxs[i]<<"]\n";
    }
    if (m_dists[i]!=NULL && (eval-1==m_histos[4][i] || m_histos[4][i]<0))
      m_dists[i]->Insert(val,weight,ncount);
    if (pass) return 1;
    if (m_items[i][0]>=0 || m_vars[i]->IDName()=="Count") return 0;
    for (size_t l(0);l<moms.size();++l)
      for (Particle_List::iterator pit(reflist.begin());
	   pit!=reflist.end();++pit) 
	if (*pit==moms[l]) {
	  msg_Debugging()<<"  erase "<<**pit<<"\n";
	  reflist.erase(pit);
	  break;
	}
    moms.clear();
    return -1;
  }
  if (j>0 && m_flavs[i][j]!=m_flavs[i][j-1]) o=k=0;
  for (;k<reflist.size();++k) {
    if (m_flavs[i][j].Includes(reflist[k]->Flav())) {
      if ((m_items[i][j]<0 && -int(o)<=m_items[i][j]+1) || int(o)==m_items[i][j]) {
	moms.push_back(reflist[k]);
	int stat(Evaluate(reflist,weight,ncount,moms,i,j+1,k+1,o+1,eval));
	if (stat<1) return stat;
      	if (int(o)==m_items[i][j]) return 1;
	moms.pop_back();
      }
      ++o;
    }
  }
  return 1;
}

void One_Variable_Selector::Evaluate
(const ATOOLS::Particle_List &inlist,ATOOLS::Particle_List &outlist,
 double weight, double ncount) 
{
  msg_Debugging()<<METHOD<<"(): '"<<m_inlist<<"' -> '"<<m_outlist<<"' {\n";
  std::vector<Particle_List*> cndlist(m_cndlist.size());
  for (size_t i(0);i<m_cndlist.size();++i) {
    cndlist[i] = new Particle_List();
    p_ana->AddParticleList(m_cndlist[i][1],cndlist[i]);
  }
  p_flow->Insert(0.0,weight,ncount);
  Particle_List vreflist(inlist);
  size_t eval(0);
  for (int oldi(0), i(0);i<(int)m_flavs.size();++i) {
    if (i!=oldi) eval=0;
    oldi=i;
    Particle_List moms;
    int stat(Evaluate(vreflist,weight,ncount,moms,i,0,0,0,eval));
    if (m_vars[i]->IDName()=="Count") {
      std::vector<ATOOLS::Vec4D> vmoms(eval);
      double val(m_vars[i]->Value(&vmoms.front(),eval));
      bool pass(val>=m_mins[i] && val<=m_maxs[i]);
      if (msg_LevelIsDebugging()) {
	msg_Debugging()<<"  "<<m_flavs[i][0].IDName();
	for (size_t k(1);k<m_flavs[i].size();++k) 
	  msg_Debugging()<<","<<m_flavs[i][k].IDName();
	msg_Debugging()<<" "<<m_items[i][0];
	for (size_t k(1);k<m_items[i].size();++k) 
	  msg_Debugging()<<","<<m_items[i][k];
	msg_Debugging()<<" "<<m_vars[i]->Name()<<"("<<eval<<") "
		       <<(pass?"\\in":"\\nin")
		       <<" ["<<m_mins[i]<<","<<m_maxs[i]<<"]\n";
      }
      if (m_dists[i]!=NULL) m_dists[i]->Insert(val,weight,ncount);
      if (!pass) {
	msg_Debugging()<<"} failed\n";
	return;
      }
    }
    else if (!eval && msg_LevelIsDebugging()) {
      msg_Debugging()<<"  "<<m_flavs[i][0].IDName();
      for (size_t k(1);k<m_flavs[i].size();++k) 
	msg_Debugging()<<","<<m_flavs[i][k].IDName();
      msg_Debugging()<<" "<<m_items[i][0];
      for (size_t k(1);k<m_items[i].size();++k) 
	msg_Debugging()<<","<<m_items[i][k];
      msg_Debugging()<<" "<<m_vars[i]->Name()<<", range ["
		     <<m_mins[i]<<","<<m_maxs[i]<<"] not checked\n";
    }
    if (stat==0) {
      msg_Debugging()<<"} failed\n";
      return;
    }
    if (stat<0) {
      --i;
      continue;
    }
    p_flow->Insert((double)i+1.5,weight,0);
  }
  msg_Debugging()<<"} passed\n";
  outlist.resize(vreflist.size());
  for (size_t i(0);i<vreflist.size();++i) 
    outlist[i] = new Particle(*vreflist[i]);
  for (size_t i(0);i<m_cndlist.size();++i) {
    Particle_List *list(p_ana->GetParticleList(m_cndlist[i][0]));
    if (list==NULL) {
      msg_Error()<<METHOD<<"(): List '"<<m_cndlist[i][0]
		 <<"' not found."<<std::endl;
      continue;
    }
    cndlist[i]->resize(list->size());
    for (size_t j(0);j<list->size();++j) 
      (*cndlist[i])[j] = new Particle(*(*list)[j]);
  }
}

std::string OVS_Tag_Replacer::ReplaceTags(std::string &expr) const
{
  THROW(fatal_error,"Invalid function call");
}

Term *OVS_Tag_Replacer::ReplaceTags(Term *term) const
{
  size_t bpos(term->Tag().find('['));
  if (bpos!=std::string::npos) {
    Blob_Data<Vec4D> *dat(dynamic_cast<Blob_Data<Vec4D>*>((*p_ana)[term->Tag()]));
    if (dat!=NULL) {
      term->Set(dat->Get());
      return term;
    }
    THROW(critical_error,"Tag '"+term->Tag()+"' not found");
  }
  Blob_Data<double> *dat(dynamic_cast<Blob_Data<double>*>((*p_ana)[term->Tag()]));
  if (dat!=NULL) {
    term->Set(dat->Get());
    return term;
  }
  THROW(critical_error,"Tag '"+term->Tag()+"' not found");
}

Analysis_Object *One_Variable_Selector::GetCopy() const
{
  return new One_Variable_Selector
    (m_inlist,m_outlist,m_cndlist,m_isobs,
     m_flavs,m_items,m_vtags,m_mins,m_maxs,
     m_histos,p_ana,m_name);
}

