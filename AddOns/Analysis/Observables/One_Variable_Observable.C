#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

#include "ATOOLS/Math/Variable.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
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

  class One_Variable_Observable: public Primitive_Observable_Base {
  private:
    Flavour_Matrix   m_flavs;
    Int_Matrix       m_items;
    String_Vector    m_vtags;
    Variable_Vector  m_vars;
    Double_Matrix    m_histos;
    Histogram_Vector m_dists;
  public:
    One_Variable_Observable
    (const std::string &inlist,const Flavour_Matrix &flavs,
     const Int_Matrix &items,const String_Vector &vtags,
     const Double_Matrix &m_histos,Primitive_Analysis *const ana,
     const std::string &name="");
    ~One_Variable_Observable();
    void EvaluateNLOcontrib(double weight,double ncount);
    void Evaluate(const ATOOLS::Particle_List &list, 
		  double weight, double ncount);
    void Evaluate(const ATOOLS::Particle_List &list, 
		  double weight, double ncount,const int mode);
    bool Evaluate(const ATOOLS::Particle_List &list,
		  double weight,double ncount,Particle_List moms,
		  const size_t i,const size_t j,size_t k,
		  size_t o,size_t &eval,const int mode); 
    Analysis_Object &operator+=(const Analysis_Object &obj);
    void EvaluateNLOevt();
    void EndEvaluation(double scale=1.0);
    void Restore(double scale=1.0);
    void Output(const std::string & pname);
    Primitive_Observable_Base *Copy() const;    
  };// end of class One_Variable_Observable

} // namespace ANALYSIS

using namespace ANALYSIS;

DECLARE_GETTER(One_Variable_Observable,"VarObs",
	       Primitive_Observable_Base,Analysis_Key);

void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,One_Variable_Observable>::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList: <list>,\n"
     <<std::setw(width+7)<<" "<<"Flavs: [[flav11, ..], .. [flavN1, ..],\n"
     <<std::setw(width+7)<<" "<<"Items: [[item11, ..], .. [itemN1, ..],\n"
     <<std::setw(width+7)<<" "<<"Vars:  [var1, ..],\n"
     <<std::setw(width+7)<<" "<<"Types: [type1, ..],\n"
     <<std::setw(width+7)<<" "<<"Bins:  [bins1, ..],\n"
     <<std::setw(width+7)<<" "<<"Mins:  [min1,  ..],\n"
     <<std::setw(width+7)<<" "<<"Maxs:  [max1,  ..]\n"
     <<std::setw(width+4)<<" "<<"}";
}

Primitive_Observable_Base *
ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,One_Variable_Observable>::operator()
  (const Analysis_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  s.DeclareVectorSettingsWithEmptyDefault({ "Vars", "Types", "Bins" });
  s.DeclareMatrixSettingsWithEmptyDefault({ "Flavs", "Items" });
  Flavour_Matrix flavs;
  Int_Matrix items;
  String_Vector vtags;
  Double_Matrix histos(4);

  const auto inlist = s["InList"].SetDefault("FinalState").GetVector<std::string>()[0];

  auto cfls = s["Flavs"].GetMatrix<int>();
  for (const auto& cfl : cfls) {
    flavs.push_back(Flavour_Vector(cfl.size()));
    for (size_t k(0);k<cfl.size();++k) {
      flavs.back()[k]=Flavour((kf_code)abs(cfl[k]));
      if (cfl[k]<0) flavs.back()[k]=flavs.back()[k].Bar();
    }
  }

  const auto cits = s["Items"].GetMatrix<int>();
  for (const auto& cit : cits) {
    items.push_back(cit);
  }

  const auto vars = s["Vars"].GetVector<std::string>();
  for (const auto& var : vars) {
    vtags.push_back(var);
  }

  const auto types = s["Types"].GetVector<std::string>();
  for (const auto& type : types) {
    histos[0].push_back(HistogramType(type));
  }

  const auto bins = s["Bins"].GetVector<double>();
  for (const auto& bin : bins) {
    histos[1].push_back(bin);
  }

  const auto mins = s["Mins"].GetVector<double>();
  for (const auto& min : mins) {
    histos[2].push_back(min);
  }

  const auto maxs = s["Maxs"].GetVector<double>();
  for (const auto& max : maxs) {
    histos[3].push_back(max);
  }

  if (flavs.empty() || items.empty() || vtags.empty()) {
    msg_Debugging()<<METHOD<<"(): Cannot initialize selector.\n";
    return NULL;
  }
  if (histos[0].empty()) histos[0].push_back(-1.0);
  size_t max(Max(vtags.size(),Max(flavs.size(),items.size())));
  for (size_t i(flavs.size());i<max;++i) flavs.push_back(flavs.back());
  for (size_t i(items.size());i<max;++i) items.push_back(items.back());
  for (size_t i(vtags.size());i<max;++i) vtags.push_back(vtags.back());
  for (size_t i(histos[0].size());i<max;++i) histos[0].push_back(-1);
  for (size_t i(flavs.size());i<max;++i) {
    max=Max(flavs[i].size(),items[i].size());
    for (size_t j(flavs[i].size());j<max;++j) 
      flavs[i].push_back(flavs[i].back());
    for (size_t j(items[i].size());j<max;++j) 
      items[i].push_back(items[i].back());
  }
  return new One_Variable_Observable(
      inlist,flavs,items,vtags,histos,key.p_analysis);
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ATOOLS;

One_Variable_Observable::One_Variable_Observable
(const std::string &inlist,const Flavour_Matrix &flavs,
 const Int_Matrix &items,const String_Vector &vtags,
 const Double_Matrix &histos,Primitive_Analysis *const ana,
 const std::string &name):
  m_flavs(flavs), m_items(items), m_vtags(vtags), 
  m_histos(histos), m_dists(flavs.size(),NULL)
{
  m_listname=inlist;
  msg_Debugging()<<METHOD<<"(): {\n";
  m_vars.resize(m_vtags.size(),NULL);
  for (size_t i(0);i<m_vtags.size();++i) {
    m_vars[i]=ATOOLS::Variable_Getter::GetObject(m_vtags[i],m_vtags[i]);
    if (m_vars[i]==NULL) THROW
      (fatal_error,"Variable '"+m_vtags[i]+"' does not exist. Run 'Sherpa"+
       std::string(" SHOW_Analysis_SYNTAX=1' to list variables."));
  }
  if (name!="") m_name=name;
  else {
    size_t n(0);
    while (ana->GetObject(m_listname+"_"+ToString(n))!=NULL) ++n;
    m_name=m_listname+"_"+ToString(n);
  }
  for (size_t i(0);i<m_dists.size();++i)
    if (m_histos[0][i]>-1) {
      msg_Debugging()<<"  init histo "<<i<<" for ";
      for (size_t j(0);j<m_flavs[i].size();++j) 
	msg_Debugging()<<m_flavs[i][j].IDName()<<" "<<m_items[i][j]<<" ";
      msg_Debugging()<<"-> type "<<m_histos[0][i]<<", "<<m_histos[1][i]
		     <<" bins, min "<<m_histos[2][i]<<", max "
		     <<m_histos[3][i]<<"\n";
      m_dists[i] = new ATOOLS::Histogram
	((int)m_histos[0][i],m_histos[2][i],
	 m_histos[3][i],(int)m_histos[1][i]);
    }
  msg_Debugging()<<"}\n";
}

Analysis_Object &One_Variable_Observable::operator+=
(const Analysis_Object &obj)
{
  const One_Variable_Observable *vob((const One_Variable_Observable*)&obj);
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) *m_dists[i]+=*vob->m_dists[i];
  return *this;
}

void One_Variable_Observable::EvaluateNLOevt()
{
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) m_dists[i]->FinishMCB();
}

void One_Variable_Observable::EndEvaluation(double scale) 
{
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) {
      m_dists[i]->MPISync();
      m_dists[i]->Finalize();
      if (scale!=1.0) m_dists[i]->Scale(scale);
    }
}

void One_Variable_Observable::Restore(double scale) 
{
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) {
      if (scale!=1.0) m_dists[i]->Scale(scale);
      m_dists[i]->Restore();
    }
}

void One_Variable_Observable::Output(const std::string & pname) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  std::string bname(pname+"/"+m_name);
  std::set<std::string> names;
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) {
      std::string name(bname+"_"+m_vars[i]->IDName());
      for (size_t j(0);j<m_flavs[i].size();++j) 
	name+="_"+m_flavs[i][j].IDName()+ToString(m_items[i][j]);
      int n(0);
      std::string id;
      while (names.find(name+id)!=names.end()) id="_"+ToString(++n);
      msg_Debugging()<<"  write '"<<name<<id<<".dat'\n";
      m_dists[i]->Output((name+id+".dat").c_str());
      names.insert(name+id);
    }
  msg_Debugging()<<"}\n";
}

One_Variable_Observable::~One_Variable_Observable()
{
  while (m_vars.size()) {
    delete m_vars.back();
    m_vars.pop_back();
  }
  while (m_dists.size()) {
    if (m_dists.back()!=NULL) delete m_dists.back();
    m_dists.pop_back();
  }
}

bool One_Variable_Observable::Evaluate
(const ATOOLS::Particle_List &list,double weight,double ncount,
 Particle_List moms,const size_t i,const size_t j,size_t k,
 size_t o,size_t &eval,const int mode) 
{
  if (j>=m_flavs[i].size()) {
    ++eval;
    if (m_vars[i]->IDName()=="Count") return true;
    std::vector<ATOOLS::Vec4D> vmoms(moms.size());
    for (size_t l(0);l<vmoms.size();++l) vmoms[l]=moms[l]->Momentum();
    double val(m_vars[i]->Value(&vmoms.front(),vmoms.size()));
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
      msg_Debugging()<<") = "<<val<<"\n";
    }
    if (m_dists[i]!=NULL) {
      if (mode==1) m_dists[i]->InsertMCB(val,weight,ncount);
      else m_dists[i]->Insert(val,weight,ncount);
    }
    return true;
  }
  if (j>0 && m_flavs[i][j]!=m_flavs[i][j-1]) o=k=0;
  bool pass(false);
  for (;k<list.size();++k) {
    if (m_flavs[i][j].Includes(list[k]->Flav())) {
      if ((m_items[i][j]<0 && -int(o)<=m_items[i][j]+1) || int(o)==m_items[i][j]) {
	moms.push_back(list[k]);
	if (Evaluate(list,weight,ncount,moms,i,j+1,k+1,o+1,eval,mode)) pass=true;
	if (int(o)==m_items[i][j]) return pass;
	moms.pop_back();
      }
      ++o;
    }
  }
  return pass;
}

void One_Variable_Observable::EvaluateNLOcontrib(double weight,double ncount)
{
  Evaluate(*p_ana->GetParticleList(m_listname),weight,ncount,1);
}

void One_Variable_Observable::Evaluate
(const ATOOLS::Particle_List &list,double weight,double ncount)
{
  Evaluate(list,weight,ncount,0);
}

void One_Variable_Observable::Evaluate
(const ATOOLS::Particle_List &list,double weight,double ncount,const int mode)
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_flavs.size();++i) {
    size_t eval(0);
    if (!Evaluate(list,weight,ncount,Particle_List(),i,0,0,0,eval,mode))
      if (m_dists[i]!=NULL) {
	if (mode==1) m_dists[i]->InsertMCB(1.0,0.0,ncount);
	else m_dists[i]->Insert(1.0,0.0,ncount);
      }
    if (m_vars[i]->IDName()=="Count") {
      std::vector<ATOOLS::Vec4D> vmoms(eval);
      double val(m_vars[i]->Value(&vmoms.front(),eval));
      if (msg_LevelIsDebugging()) {
	msg_Debugging()<<"  "<<m_flavs[i][0].IDName();
	for (size_t k(1);k<m_flavs[i].size();++k) 
	  msg_Debugging()<<","<<m_flavs[i][k].IDName();
	msg_Debugging()<<" "<<m_items[i][0];
	for (size_t k(1);k<m_items[i].size();++k) 
	  msg_Debugging()<<","<<m_items[i][k];
	msg_Debugging()<<" "<<m_vars[i]->Name()<<"("<<eval<<")\n";
      }
      if (m_dists[i]!=NULL) {
	if (mode==1) m_dists[i]->InsertMCB(val,weight,ncount);
	else m_dists[i]->Insert(val,weight,ncount);
      }
    }
  }
  msg_Debugging()<<"} done\n";
}

Primitive_Observable_Base *One_Variable_Observable::Copy() const
{
  return new One_Variable_Observable
    (m_listname,m_flavs,m_items,m_vtags,m_histos,p_ana,m_name);
}

