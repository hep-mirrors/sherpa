#include "Jet_Cone_Distribution.H"
#include "MyStrStream.H"

using namespace ANALYSIS;

Jet_Cone_Distribution::Jet_Cone_Distribution(const int linlog, const double Etcut, 
					     const double etamin, const double etamax, 
					     const double Rmin, const double Rmax, 
					     const int nbins, 
					     Primitive_Calorimeter * const calorimeter) :
  Primitive_Observable_Base(linlog,Rmin,Rmax,nbins,NULL), 
  m_Etcut(Etcut), p_calorimeter(calorimeter)
{
  std::string etname;
  MyStrStream s1;
  s1<<m_Etcut;
  s1>>etname;
  m_name = std::string("ConeNumb_")+etname;
  double dx = (m_xmax-m_xmin)/double(m_nbins);
  for (int i=0;i<nbins;i++) {
    m_cones.push_back(new Calorimeter_Cone(Etcut,m_xmin+i*dx,p_calorimeter));
    m_cones[i]->SetEtaRangeForJets(etamin,etamax,1);
    m_histos.push_back(new ATOOLS::Histogram(0,0.,10.,nbins));
  }
}

Jet_Cone_Distribution::Jet_Cone_Distribution(const int linlog, const double Etcut, 
					     const double Rmin, const double Rmax, 
					     const int nbins, 
					     Primitive_Calorimeter * const calorimeter) :
  Primitive_Observable_Base(linlog,Rmin,Rmax,nbins,NULL), 
  m_Etcut(Etcut), p_calorimeter(calorimeter)
{
  std::string etname;
  MyStrStream s1;
  s1<<m_Etcut;
  s1>>etname;
  m_name = std::string("ConeNumb_")+etname;
  double dx = (m_xmax-m_xmin)/double(m_nbins);
  for (int i=0;i<nbins;i++) {
    m_cones.push_back(new Calorimeter_Cone(Etcut,m_xmin+i*dx,p_calorimeter));
    m_histos.push_back(new ATOOLS::Histogram(0,0.,10.,nbins));
  }
}

Jet_Cone_Distribution::~Jet_Cone_Distribution() 
{
  int size = m_cones.size();
  for (int i=0;i<size;++i) {
    if (m_cones[size-i-1])  { delete m_cones[size-i-1];  m_cones.pop_back();  } 
    if (m_histos[size-i-1]) { delete m_histos[size-i-1]; m_histos.pop_back(); } 
  }
}

Primitive_Observable_Base * Jet_Cone_Distribution::Copy() const 
{
  return new Jet_Cone_Distribution(m_type,m_Etcut,m_xmin,m_xmax,m_nbins,p_calorimeter);
}

void Jet_Cone_Distribution::EndEvaluation(double scale) 
{
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Jet_Cone_Distribution::Reset()
{
  p_histo->Reset();
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Reset();
  }  
}

void Jet_Cone_Distribution::Output(const std::string & pname) {
  int  mode_dir = 448;
  mkdir((pname).c_str(),mode_dir); 
  for (size_t i=0; i<m_histos.size();i++) {
    std::string fname;
    MyStrStream s;
    s<<m_cones[i]->Radius();
    s<<".dat"; 
    s>>fname;
    m_histos[i]->Output((pname+std::string("/")+m_name+std::string("_")+fname).c_str());
  }
  p_histo->Output((pname+std::string("/")+m_name+std::string(".dat")).c_str());
}


void Jet_Cone_Distribution::Fill(double weight, int ncount)
{
  int NofJets;
  for (unsigned int i=0;i<m_cones.size();++i) {
    m_cones[i]->ConstructJets();
    NofJets = m_cones[i]->NumberOfJets();
    m_histos[i]->Insert(NofJets,weight,ncount);
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Jet_Cone_Dependence::Jet_Cone_Dependence(const int linlog, const double Etcut, 
					 const double etamin, const double etamax, 
					 const double Rmin, const double Rmax, 
					 const int njetmin, const int njetmax, 
					 const int nbins, 
					 Primitive_Calorimeter * const calorimeter) :
  Primitive_Observable_Base(linlog,Rmin,Rmax,nbins,NULL), 
  m_Etcut(Etcut), m_njetmin(njetmin), m_njetmax(njetmax), p_calorimeter(calorimeter)
{
  std::string etname;
  MyStrStream s1;
  s1<<m_Etcut;
  s1>>etname;
  m_name = std::string("ConeDep_")+etname;
  double dx = (m_xmax-m_xmin)/double(m_nbins);
  for (int i=0;i<nbins;i++) {
    m_cones.push_back(new Calorimeter_Cone(Etcut,m_xmin+i*dx,p_calorimeter));
    m_cones[i]->SetEtaRangeForJets(etamin,etamax,1);
  }
  for (int i=0;i<m_njetmax-m_njetmin;i++) {
    m_histos.push_back(new ATOOLS::Histogram(0,m_xmin,m_xmax+dx,nbins+1));
  }
}

Jet_Cone_Dependence::Jet_Cone_Dependence(const int linlog, const double Etcut, 
					 const double Rmin, const double Rmax, 
					 const int njetmin, const int njetmax, 
					 const int nbins, 
					 Primitive_Calorimeter * const calorimeter) :
  Primitive_Observable_Base(linlog,Rmin,Rmax,nbins,NULL), 
  m_Etcut(Etcut), m_njetmin(njetmin), m_njetmax(njetmax), p_calorimeter(calorimeter)
{
  std::string etname;
  MyStrStream s1;
  s1<<m_Etcut;
  s1>>etname;
  m_name = std::string("ConeDep_")+etname;
  double dx = (m_xmax-m_xmin)/double(m_nbins);
  for (int i=0;i<nbins;i++) {
    m_cones.push_back(new Calorimeter_Cone(Etcut,m_xmin+i*dx,p_calorimeter));
  }
  for (int i=0;i<m_njetmax-m_njetmin;i++) {
    m_histos.push_back(new ATOOLS::Histogram(0,m_xmin,m_xmax+dx,nbins+1));
  }
}

Jet_Cone_Dependence::~Jet_Cone_Dependence() 
{
  int size = m_cones.size();
  for (int i=0;i<size;++i) {
    if (m_cones[size-i-1])  { delete m_cones[size-i-1];  m_cones.pop_back();  } 
    if (m_histos[size-i-1]) { delete m_histos[size-i-1]; m_histos.pop_back(); } 
  }
}

Primitive_Observable_Base * Jet_Cone_Dependence::Copy() const 
{
  return new Jet_Cone_Dependence(m_type,m_Etcut,m_xmin,m_xmax,
				 m_njetmin,m_njetmax,m_nbins,p_calorimeter);
}

void Jet_Cone_Dependence::EndEvaluation(double scale) 
{
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Jet_Cone_Dependence::Reset()
{
  p_histo->Reset();
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Reset();
  }  
}

void Jet_Cone_Dependence::Output(const std::string & pname) {
  int  mode_dir = 448;
  mkdir((pname).c_str(),mode_dir); 
  for (size_t i=0; i<m_histos.size();++i) {
    std::string fname;
    MyStrStream s;
    s<<m_njetmin+i;
    s<<".dat"; 
    s>>fname;
    m_histos[i]->Output((pname+std::string("/")+m_name+std::string("_")+fname).c_str());
  }
  p_histo->Output((pname+std::string("/")+m_name+std::string(".dat")).c_str());
}


void Jet_Cone_Dependence::Fill(double weight, int ncount)
{
  int NofJets;
  for (unsigned int i=0;i<m_cones.size();++i) {
    m_cones[i]->ConstructJets();
    NofJets = m_cones[i]->NumberOfJets();
    if (NofJets<m_njetmax) m_histos[NofJets-m_njetmin]->Insert(m_cones[i]->Radius(),weight,ncount);
  }
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Jet_Cone_Shape::Jet_Cone_Shape(const int linlog,const double Rmin, const double Rmax,
			       const int jetmin, const int jetmax, const int nbins, 
			       Calorimeter_Cone * cone) :
  Primitive_Observable_Base(linlog,Rmin,Rmax,nbins,NULL), 
  m_jetmin(jetmin), m_jetmax(jetmax), p_cone(cone)
{
  std::string rname,etname;
  MyStrStream s;
  s<<cone->Radius();
  s>>rname;
  MyStrStream s1;
  s1<<cone->Et_cut();
  s1>>etname;
  m_name = std::string("ConeShape_")+etname+std::string("_")+rname;
  for (int i=jetmin;i<jetmax;i++) {
    m_histos.push_back(new ATOOLS::Histogram(linlog,Rmin,Rmax,nbins));
  }
}

Jet_Cone_Shape::~Jet_Cone_Shape()
{
  int size = m_histos.size();
  for (int i=0;i<size;++i) {
    if (m_histos[size-i-1]) { delete m_histos[size-i-1]; m_histos.pop_back(); } 
  }
}

Primitive_Observable_Base * Jet_Cone_Shape::Copy() const 
{
  return new Jet_Cone_Shape(m_type,m_xmin,m_xmax,m_jetmin,m_jetmax,m_nbins,p_cone);
}

void Jet_Cone_Shape::Reset()
{
  p_histo->Reset();
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Reset();
  }  
}


void Jet_Cone_Shape::Output(const std::string & pname)
{
  int  mode_dir = 448;
  mkdir((pname).c_str(),mode_dir); 
  for (size_t i=0; i<m_histos.size();++i) {
    std::string fname;
    MyStrStream s;
    s<<i+m_jetmin;
    s<<".dat"; 
    s>>fname;
    m_histos[i]->Output((pname+std::string("/")+m_name+std::string("_")+fname).c_str());
  }
  p_histo->Output((pname+std::string("/")+m_name+std::string(".dat")).c_str());
}

void Jet_Cone_Shape::EndEvaluation(double scale)
{
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Jet_Cone_Shape::Fill(double weight,int ncount)
{
  p_cone->ConstructJets();
  for (unsigned int i=0; i<m_histos.size();++i) Fill(i,weight,ncount);
}


void Jet_Cone_Shape::Fill(int jetno,double weight,int ncount)
{
  p_cone->FillShape(jetno+m_jetmin,m_histos[jetno],weight,ncount);
}

