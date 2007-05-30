#include "Primitive_Observable_Base.H"

#include "Variable.H"
#include "Data_Reader.H"
#include "MyStrStream.H"
#include "Message.H"
#include "Histogram.H"
#include "Shell_Tools.H"
#include "Exception.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  typedef std::vector<int>        Int_Vector;
  typedef std::vector<Int_Vector> Int_Matrix;

  typedef std::vector<double>        Double_Vector;
  typedef std::vector<Double_Vector> Double_Matrix;

  typedef std::vector<Flavour>        Flavour_Vector;
  typedef std::vector<Flavour_Vector> Flavour_Matrix;

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
    void Evaluate(const ATOOLS::Particle_List &list, 
		  double weight, int ncount);
    bool Evaluate(const ATOOLS::Particle_List &list,
		  double weight,int ncount,std::vector<ATOOLS::Vec4D> moms,
		  const size_t i,const size_t j,size_t k); 
    Analysis_Object &operator+=(const Analysis_Object &obj);
    void EndEvaluation(double scale);
    void Output(const std::string & pname);
    Primitive_Observable_Base *Copy() const;    
  };// end of class One_Variable_Observable

} // namespae ANALYSIS

using namespace ANALYSIS;

DECLARE_GETTER(One_Variable_Observable_Getter,"VarObs",
 	       Primitive_Observable_Base,Argument_Matrix);

void One_Variable_Observable_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "
     <<"Tags    flavi1,.. itemi1,.. vari mini maxi binsi typei\n"
     <<std::setw(width+7)<<" "<<"Flavs   flav11,.. .. flavN1,..\n"
     <<std::setw(width+7)<<" "<<"Items   item11,.. .. itemN1,..\n"
     <<std::setw(width+7)<<" "<<"Vars    var1      .. varN\n"
     <<std::setw(width+7)<<" "<<"Types   type1     .. typeN\n"
     <<std::setw(width+7)<<" "<<"Bins    bins1     .. binsN\n"
     <<std::setw(width+7)<<" "<<"Mins    min1      .. minN\n"
     <<std::setw(width+7)<<" "<<"Maxs    max1      .. maxN\n"
     <<std::setw(width+4)<<" "<<"}";
}

Primitive_Observable_Base *
One_Variable_Observable_Getter::operator()
  (const Argument_Matrix &parameters) const
{
  std::string inlist("FinalState");
  Flavour_Matrix flavs;
  Int_Matrix items;
  String_Vector vtags;
  Double_Matrix histos(4);
  Data_Reader reader(",",";","!","=");
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
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
	flavs.back()[k]=Flavour((kf::code)abs(cfl[k]));
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
      histos[0].push_back(HistogramType(cur[7]));
      histos[1].push_back(ToType<double>(cur[6]));
      histos[2].push_back(ToType<double>(cur[4]));
      histos[3].push_back(ToType<double>(cur[5]));
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
	  flavs.back()[k]=Flavour((kf::code)abs(cfl[k]));
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
    else if (cur[0]=="Types" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	histos[0].push_back(HistogramType(cur[j]));
    }
    else if (cur[0]=="Bins" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	histos[1].push_back(ToType<double>(cur[j]));
    }
    else if (cur[0]=="Mins" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	histos[2].push_back(ToType<double>(cur[j]));
    }
    else if (cur[0]=="Maxs" && cur.size()>1) {
      for (size_t j(1);j<cur.size();++j) 
	histos[3].push_back(ToType<double>(cur[j]));
    }
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
  return new One_Variable_Observable
    (inlist,flavs,items,vtags,histos,parameters());
}

#include "Primitive_Analysis.H"

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
       std::string(" SHOW_ANALYSIS_SYNTAX=1' to list variables."));
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

void One_Variable_Observable::EndEvaluation(double scale) 
{
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) {
      m_dists[i]->Finalize();
      if (scale!=1.0) m_dists[i]->Scale(scale);
    }
}

void One_Variable_Observable::Output(const std::string & pname) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  std::string bname(pname+"/"+m_name);
  ATOOLS::MakeDir(pname.c_str(),448); 
  for (size_t i(0);i<m_dists.size();++i) 
    if (m_dists[i]!=NULL) {
      std::string name(bname+"_"+m_vars[i]->IDName());
      for (size_t j(0);j<m_flavs[i].size();++j) 
	name+="_"+m_flavs[i][j].IDName()+ToString(m_items[i][j]);
      msg_Debugging()<<"  write '"<<name<<".dat'\n";
      m_dists[i]->Output((name+".dat").c_str());
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
(const ATOOLS::Particle_List &list,double weight,int ncount,
 std::vector<ATOOLS::Vec4D> moms,const size_t i,const size_t j,size_t k) 
{
  if (j>=m_flavs[i].size()) {
    double val(m_vars[i]->Value(&moms.front(),moms.size()));
    msg_Debugging()<<"  "<<m_vars[i]->Name()<<"("<<moms.front();
    for (size_t k(1);k<moms.size();++k) msg_Debugging()<<","<<moms[k];
    msg_Debugging()<<") = "<<val<<"\n";
    if (m_dists[i]!=NULL) m_dists[i]->Insert(val,weight,ncount);
    return true;
  }
  int o(-1);
  bool pass(false);
  for (;k<list.size();++k) {
    if (list[k]->Flav()==m_flavs[i][j]) {
      ++o;
      if (m_items[i][j]<0 || o==m_items[i][j]) {
	moms.push_back(list[k]->Momentum());
	if (Evaluate(list,weight,ncount,moms,i,j+1,k+1)) pass=true;
	if (o==m_items[i][j]) return pass;
      }
    }
  }
  return pass;
}

void One_Variable_Observable::Evaluate
(const ATOOLS::Particle_List &list,double weight,int ncount)
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_flavs.size();++i)
    if (!Evaluate(list,weight,ncount,std::vector<Vec4D>(),i,0,0))
      if (m_dists[i]!=NULL) m_dists[i]->Insert(1.0,0.0,ncount);
  msg_Debugging()<<"} done\n";
}

Primitive_Observable_Base *One_Variable_Observable::Copy() const
{
  return new One_Variable_Observable
    (m_listname,m_flavs,m_items,m_vtags,m_histos,p_ana,m_name);
}

