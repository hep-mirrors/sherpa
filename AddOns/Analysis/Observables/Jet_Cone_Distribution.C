#include "AddOns/Analysis/Observables/Jet_Cone_Distribution.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "AddOns/Analysis/Triggers/Final_Selector.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Jet_Cone_Distribution_Getter,"JetConeDist",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base * 
Jet_Cone_Distribution_Getter::operator()(const Argument_Matrix &parameters) const
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    return new Jet_Cone_Distribution(HistogramType(parameters[0][6]),
				     ATOOLS::ToType<double>(parameters[0][0]),
				     ATOOLS::ToType<double>(parameters[0][1]),
				     ATOOLS::ToType<double>(parameters[0][2]),
				     ATOOLS::ToType<double>(parameters[0][3]),
				     ATOOLS::ToType<double>(parameters[0][4]),
				     ATOOLS::ToType<int>(parameters[0][5]),
				     parameters());
  }
  else if (parameters.size()<7) return NULL;
  double etcut=0.0, etamin=-10., etamax=10., rmin=0., rmax=10.;
  size_t bins=100;
  std::string scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="ETCUT") etcut=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="ETAMIN") etamin=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="ETAMAX") etamax=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="RMIN")   rmin=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="RMAX")   rmax=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE")  scale=parameters[i][1];
    else if (parameters[i][0]=="NBINS")  bins=ATOOLS::ToType<int>(parameters[i][1]);
  }
  return new Jet_Cone_Distribution(HistogramType(scale),etcut,etamin,etamax,rmin,rmax,bins,parameters());
}									

void Jet_Cone_Distribution_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"etcut etamin etamax rmin rmax bins Lin|LinErr|Log|LogErr";
}

Jet_Cone_Distribution::Jet_Cone_Distribution(const int linlog, const double Etcut, 
					     const double etamin, const double etamax, 
					     const double Rmin, const double Rmax, 
					     const int nbins, 
					     Primitive_Analysis * const ana) :
  Primitive_Observable_Base(linlog,Rmin,Rmax,nbins), 
  m_Etcut(Etcut),m_etamin(etamin),m_etamax(etamax)
{
  p_ana = ana;

  std::string etname;
  MyStrStream s1;
  s1<<m_Etcut;
  s1>>etname;
  m_name = std::string("ConeNumb_")+etname;
  double dx = (m_xmax-m_xmin)/double(m_nbins);
  for (int i=0;i<nbins;i++) {
    m_cones.push_back(new Calorimeter_Cone(Etcut,m_etamin,m_etamax,m_xmin+i*dx));
    m_cones.back()->SetAnalysis(p_ana);
    m_cones[i]->SetEtaRangeForJets(m_etamin,m_etamax,1);
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
  return new Jet_Cone_Distribution(m_type,m_Etcut,m_etamin,m_etamax,
				   m_xmin,m_xmax,m_nbins,p_ana);
}

void Jet_Cone_Distribution::EndEvaluation(double scale) 
{
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->MPISync();
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Jet_Cone_Distribution::Restore(double scale) 
{
  for (size_t i=0; i<m_histos.size();++i) {
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Restore();
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
  ATOOLS::MakeDir(pname); 
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


void Jet_Cone_Distribution::Fill(double weight, double ncount)
{
  int NofJets;
  for (unsigned int i=0;i<m_cones.size();++i) {
    m_cones[i]->ConstructJets();
    NofJets = m_cones[i]->NumberOfJets();
    m_histos[i]->Insert(NofJets,weight,ncount);
    //insert zero in others
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DECLARE_GETTER(Jet_Cone_Dependence_Getter,"JetConeDep",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base * 
Jet_Cone_Dependence_Getter::operator()(const Argument_Matrix &parameters) const
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<9) return NULL;
    return new Jet_Cone_Dependence(HistogramType(parameters[0][8]),
				   ATOOLS::ToType<double>(parameters[0][0]),
				   ATOOLS::ToType<double>(parameters[0][1]),
				   ATOOLS::ToType<double>(parameters[0][2]),
				   ATOOLS::ToType<double>(parameters[0][3]),
				   ATOOLS::ToType<double>(parameters[0][4]),
				   ATOOLS::ToType<int>(parameters[0][5]),
				   ATOOLS::ToType<int>(parameters[0][6]),
				   ATOOLS::ToType<int>(parameters[0][7]),
				   parameters());
  }
  else if (parameters.size()<9) return NULL;
  double etcut=0.0, etamin=-10., etamax=10., rmin=0., rmax=10.;
  size_t bins=100,nmin=1,nmax=10;
  std::string scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="ETCUT") etcut=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="ETAMIN") etamin=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="ETAMAX") etamax=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="RMIN")   rmin=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="RMAX")   rmax=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="NMIN")   nmin=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="NMAX")   nmax=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="NBINS")  bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE")  scale=parameters[i][1];
  }
  return new Jet_Cone_Dependence(HistogramType(scale),etcut,etamin,etamax,rmin,rmax,
				 nmin,nmax,bins,parameters());
}									

void Jet_Cone_Dependence_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"etcut etamin etamax rmin rmax nmin nmax bins Lin|LinErr|Log|LogErr";
}


Jet_Cone_Dependence::Jet_Cone_Dependence(const int linlog, const double Etcut, 
					 const double etamin, const double etamax, 
					 const double Rmin, const double Rmax, 
					 const int njetmin, const int njetmax, 
					 const int nbins, 
					 Primitive_Analysis * const ana) :
  Primitive_Observable_Base(linlog,Rmin,Rmax,nbins), 
  m_Etcut(Etcut), m_etamin(etamin), m_etamax(etamax), m_njetmin(njetmin), m_njetmax(njetmax)
{
  p_ana=ana;
  std::string etname;
  MyStrStream s1;
  s1<<m_Etcut;
  s1>>etname;
  m_name = std::string("ConeDep_")+etname;
  double dx = (m_xmax-m_xmin)/double(m_nbins);
  for (int i=0;i<nbins;i++) {
    m_cones.push_back(new Calorimeter_Cone(Etcut,m_etamin,m_etamax,m_xmin+i*dx));
    m_cones.back()->SetAnalysis(p_ana);
    m_cones[i]->SetEtaRangeForJets(m_etamin,m_etamax,1);
  }
  for (int i=0;i<m_njetmax-m_njetmin;i++) {
    m_histos.push_back(new ATOOLS::Histogram(0,m_xmin,m_xmax+dx,nbins+1));
  }
}

Jet_Cone_Dependence::~Jet_Cone_Dependence() 
{
  int sizec = m_cones.size();
  for (int i=0;i<sizec;++i) {
    if (m_cones[sizec-i-1])  { delete m_cones[sizec-i-1];  m_cones.pop_back();  } 
  }
  int sizeh = m_histos.size(); 
  for (int i=0;i<sizeh;++i) {
  if (m_histos[sizeh-i-1]) { delete m_histos[sizeh-i-1]; m_histos.pop_back(); } 
  }
}

Primitive_Observable_Base * Jet_Cone_Dependence::Copy() const 
{
  /*
  return new Jet_Cone_Dependence(m_type,m_Etcut,m_xmin,m_xmax,m_etamin,m_etamax,
				 m_njetmin,m_njetmax,m_nbins,p_ana);
  */
  return new Jet_Cone_Dependence(m_type,m_Etcut,m_etamin,m_etamax,m_xmin,m_xmax,
				 m_njetmin,m_njetmax,m_nbins,p_ana);
}

void Jet_Cone_Dependence::EndEvaluation(double scale) 
{
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->MPISync();
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Jet_Cone_Dependence::Restore(double scale) 
{
  for (size_t i=0; i<m_histos.size();++i) {
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Restore();
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
  ATOOLS::MakeDir(pname); 
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


void Jet_Cone_Dependence::Fill(double weight, double ncount)
{
  int NofJets;
  for (unsigned int i=0;i<m_cones.size();++i) {
    m_cones[i]->ConstructJets();
    NofJets = m_cones[i]->NumberOfJets();
    if (NofJets<m_njetmax) {
      m_histos[NofJets-m_njetmin]->Insert(m_cones[i]->Radius(),weight,ncount);
      for (size_t j=1; j<m_histos.size();++j) {
	if ((int)j != (NofJets-m_njetmin)) m_histos[j]->Insert(0.,0.,ncount);
      }
    }
    else {
      for (size_t j=1; j<m_histos.size();++j) {
	m_histos[j]->Insert(0.,0.,ncount);
      }
    }
  }
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


DECLARE_GETTER(Jet_Cone_Shape_Getter,"JetConeShape",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base * 
Jet_Cone_Shape_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<10) return NULL;
    return new Jet_Cone_Shape(HistogramType(parameters[0][9]),
			      ATOOLS::ToType<double>(parameters[0][0]),
			      ATOOLS::ToType<double>(parameters[0][1]),
			      ATOOLS::ToType<double>(parameters[0][2]),
			      ATOOLS::ToType<double>(parameters[0][3]),
			      ATOOLS::ToType<double>(parameters[0][4]),
			      ATOOLS::ToType<double>(parameters[0][5]),
			      ATOOLS::ToType<int>(parameters[0][6]),
			      ATOOLS::ToType<int>(parameters[0][7]),
			      ATOOLS::ToType<int>(parameters[0][8]),
			      parameters());
  }
  else if (parameters.size()<10) return NULL;
  double etcut=0.0,radius=0.7, etamin=-10., etamax=10., rmin=0., rmax=10.;
  size_t bins=100, nmin=1, nmax=10;
  std::string scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="ETCUT")   rmin=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="RADIUS") radius=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="ETAMIN") etamin=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="ETAMAX") etamax=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="RMIN")   rmin=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="RMAX")   rmax=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE")  scale=parameters[i][1];
    else if (parameters[i][0]=="NMIN")   nmin=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="NMAX")   nmax=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="NBINS")  bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE")  scale=parameters[i][1];
  }
  return new Jet_Cone_Shape(HistogramType(scale),etcut,radius,etamin,etamax,rmin,rmax,nmin,nmax,bins,parameters());
}									

void Jet_Cone_Shape_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"etcut radius etamin etamax rmin rmax nmin nmax bins Lin|LinErr|Log|LogErr";
}


Jet_Cone_Shape::Jet_Cone_Shape(const int linlog,const double Etcut, const double radius, 
			       const double etamin, const double etamax,
			       const double Rmin, const double Rmax,
			       const int jetmin, const int jetmax, const int nbins, 
			       Primitive_Analysis * const ana) :
  Primitive_Observable_Base(linlog,Rmin,Rmax,nbins), 
  m_Etcut(Etcut), m_R(radius), m_etamin(etamin), m_etamax(etamax), m_jetmin(jetmin), m_jetmax(jetmax)
{
  p_ana=ana;
  std::string etname;
  MyStrStream s1;
  s1<<m_Etcut;
  s1>>etname;
  std::string rname;
  MyStrStream s2;
  s2<<m_R;
  s2>>rname;
  m_name = std::string("ConeShape_")+etname+std::string("_")+rname;
  
  p_cone = new Calorimeter_Cone(Etcut,m_etamin,m_etamax,m_R);
  p_cone->SetAnalysis(p_ana);
  p_cone->SetEtaRangeForJets(m_etamin,m_etamax,1);
   for (int i=jetmin;i<jetmax;i++) {
    m_histos.push_back(new ATOOLS::Histogram(linlog,Rmin,Rmax,nbins));
  }
}

Jet_Cone_Shape::~Jet_Cone_Shape()
{
  
  delete p_cone;
  
  int size = m_histos.size();
  for (int i=0;i<size;++i) {
    if (m_histos[size-i-1]) { delete m_histos[size-i-1]; m_histos.pop_back(); } 
  }
}

Primitive_Observable_Base * Jet_Cone_Shape::Copy() const 
{
  return new Jet_Cone_Shape(m_type,m_Etcut,m_R,m_etamin,m_etamax,m_xmin,m_xmax,m_jetmin,m_jetmax,m_nbins,p_ana);
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
  ATOOLS::MakeDir(pname); 
  for (size_t i=0; i<m_histos.size();++i) {
    std::string fname;
    MyStrStream s1,s2,s3;
    s1<<i+m_jetmin;
    s1<<".dat"; 
    s1>>fname;
    m_histos[i]->Output((pname+std::string("/")+m_name+
			 std::string("_")+fname).c_str());
  }
  p_histo->Output((pname+std::string("/")+m_name+std::string(".dat")).c_str());
}

void Jet_Cone_Shape::EndEvaluation(double scale)
{
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->MPISync();
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Jet_Cone_Shape::Restore(double scale) 
{
  for (size_t i=0; i<m_histos.size();++i) {
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Restore();
  }
}

void Jet_Cone_Shape::Fill(double weight,double ncount)
{
  p_cone->ConstructJets();
  for (unsigned int i=0; i<m_histos.size();++i) Fill(i,weight,ncount);
}


void Jet_Cone_Shape::Fill(int jetno,double weight,double ncount)
{
  p_cone->FillShape(jetno+m_jetmin,m_histos[jetno],weight,ncount);
}
