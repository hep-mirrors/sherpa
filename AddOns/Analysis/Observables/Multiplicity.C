#include "AddOns/Analysis/Observables/Multiplicity.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:finalstate_list;
    std::string reflist=parameters[0].size()>5?parameters[0][5]:finalstate_list;
    return new Class(HistogramType(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),list,reflist);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,CLASS>::operator()(const Argument_Matrix &parameters) const \
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"min max bins Lin|LinErr|Log|LogErr [list [reflist]]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS)					\
  DEFINE_PRINT_METHOD(CLASS)

#include "ATOOLS/Math/MathTools.H"
#include "AddOns/Analysis/Tools/Particle_Qualifier.H"

using namespace ATOOLS;

DEFINE_OBSERVABLE_GETTER(Multiplicity,"Multi")
DEFINE_OBSERVABLE_GETTER(InclMultiplicity,"InclMulti")
DEFINE_OBSERVABLE_GETTER(Hadron_Multiplicities,"Hadron_Multis")




Multiplicity::Multiplicity(int type,double xmin,double xmax,int nbins,
			   const std::string & listname,const std::string &reflist) :
  Primitive_Observable_Base(type,xmin,xmax,nbins)
{
  m_reflist=reflist;
  if (listname!="") {
    m_listname = listname;
    m_name = listname+"_"+m_reflist+"_multi.dat";
  }
  else
    m_name = "multi.dat";
}
 
void Multiplicity::Evaluate(const ATOOLS::Particle_List & pl,
			    double weight, double ncount)
{
  ATOOLS::Particle_List* ref=p_ana->GetParticleList(m_reflist);
  if (ref==NULL || ref->empty()) {
    p_histo->Insert(-1.0,0.0,ncount);
    return;
  }
  p_histo->Insert((double)pl.size(),weight,ncount); 
}

void Multiplicity::EvaluateNLOcontrib(double weight, double ncount)
{
  ATOOLS::Particle_List* ref=p_ana->GetParticleList(m_reflist);
  if (ref==NULL || ref->empty()) {
    p_histo->InsertMCB(-1.0,0.0,ncount);
    return;
  }
  Particle_List *pl=p_ana->GetParticleList(m_listname);
  p_histo->InsertMCB((double)pl->size(),weight,ncount); 
}

void Multiplicity::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * Multiplicity::Copy() const {
  return new Multiplicity(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}



InclMultiplicity::InclMultiplicity(int type,double xmin,double xmax,int nbins,
			   const std::string & listname,const std::string &reflist) :
  Primitive_Observable_Base(type,xmin,xmax,nbins)
{
  m_reflist=reflist;
  if (listname!="") {
    m_listname = listname;
    m_name = listname+"_"+m_reflist+"_inclmulti.dat";
  }
  else
    m_name = "inclmulti.dat";
}
 
void InclMultiplicity::Evaluate(const ATOOLS::Particle_List & pl,
			    double weight, double ncount)
{
  ATOOLS::Particle_List* ref=p_ana->GetParticleList(m_reflist);
  if (ref==NULL || ref->empty()) {
    p_histo->Insert(-1.0,0.0,ncount);
    return;
  }
  for (size_t i(1);i<=pl.size();++i)
    p_histo->Insert((double)i,weight,0); 
  p_histo->Insert(0.0,weight,ncount); 
}

void InclMultiplicity::EvaluateNLOcontrib(double weight, double ncount)
{
  ATOOLS::Particle_List* ref=p_ana->GetParticleList(m_reflist);
  if (ref==NULL || ref->empty()) {
    p_histo->InsertMCB(-1.0,0.0,ncount);
    return;
  }
  Particle_List *pl=p_ana->GetParticleList(m_listname);
  for (size_t i(1);i<=pl->size();++i)
    p_histo->InsertMCB((double)i,weight,0); 
  p_histo->InsertMCB(0.0,weight,ncount); 
}

void InclMultiplicity::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * InclMultiplicity::Copy() const {
  return new InclMultiplicity(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Hadron_Multiplicities::Hadron_Multiplicities(int type,double xmin,double xmax,int nbins,
					     const std::string & listname,const std::string &reflist) :
  Normalized_Observable(4,0.,100.,100)
{
  m_listname = listname; 
  m_name     = m_listname+"_Multis.dat";
}

void Hadron_Multiplicities::Evaluate(const ATOOLS::Particle_List & pl, double weight, double ncount)
{
  //Particle_List * pl = p_ana->GetParticleList(m_listname);
  Flavour flav;
  kf_code kfc;
  for (Particle_List::const_iterator pliter=pl.begin();pliter!=pl.end();pliter++) {
    flav = (*pliter)->Flav();
    if (!flav.IsHadron() && !flav.IsPhoton()) continue;
    kfc  = flav.Kfcode();
    if (kfc==kf_photon)              p_obs->Insert(99.5,weight,ncount); 
    if (kfc==kf_pi)                  p_obs->Insert(1.5,weight,ncount); 
    if (kfc==kf_pi_plus)             p_obs->Insert(2.5,weight,ncount); 
    if (kfc==kf_eta)                 p_obs->Insert(3.5,weight,ncount); 
    if (m_listname=="PrimordialHadrons") 
      if (kfc==kf_K)                 p_obs->Insert(4.5,weight,ncount); 
    //else 
      if (kfc==kf_K_L||kfc==kf_K_S||kfc==kf_K) p_obs->Insert(4.5,weight,ncount); 
    
    if (kfc==kf_K_plus)              p_obs->Insert(5.5,weight,ncount); 
    if (kfc==kf_eta_prime_958)       p_obs->Insert(6.5,weight,ncount); 

    if (kfc==kf_rho_770)             p_obs->Insert(11.5,weight,ncount); 
    if (kfc==kf_rho_770_plus)        p_obs->Insert(12.5,weight,ncount); 
    if (kfc==kf_omega_782)           p_obs->Insert(13.5,weight,ncount); 
    if (kfc==kf_K_star_892)          p_obs->Insert(14.5,weight,ncount); 
    if (kfc==kf_K_star_892_plus)     p_obs->Insert(15.5,weight,ncount); 
    if (kfc==kf_phi_1020)            p_obs->Insert(16.5,weight,ncount); 


    if (kfc==kf_D_plus)              p_obs->Insert(21.5,weight,ncount); 
    if (kfc==kf_D)                   p_obs->Insert(22.5,weight,ncount); 
    if (kfc==kf_D_s_plus)            p_obs->Insert(23.5,weight,ncount); 
    if (kfc==kf_eta_c_1S)            p_obs->Insert(24.5,weight,ncount); 
    if (kfc==kf_B_plus)              p_obs->Insert(25.5,weight,ncount); 
    if (kfc==kf_B)                   p_obs->Insert(26.5,weight,ncount); 
    if (kfc==kf_B_s)                 p_obs->Insert(27.5,weight,ncount); 
    if (kfc==kf_B_c)                 p_obs->Insert(28.5,weight,ncount); 
    if (kfc==kf_eta_b)               p_obs->Insert(29.5,weight,ncount); 

    if (kfc==kf_D_star_2010_plus)    p_obs->Insert(31.5,weight,ncount); 
    if (kfc==kf_D_star_2007)         p_obs->Insert(32.5,weight,ncount); 
    if (kfc==kf_D_s_star_plus)       p_obs->Insert(33.5,weight,ncount); 
    if (kfc==kf_J_psi_1S)            p_obs->Insert(34.5,weight,ncount); 
    if (kfc==kf_B_star_plus)         p_obs->Insert(35.5,weight,ncount); 
    if (kfc==kf_B_star)              p_obs->Insert(36.5,weight,ncount); 
    if (kfc==kf_B_s_star)            p_obs->Insert(37.5,weight,ncount); 
    if (kfc==kf_B_c_star)            p_obs->Insert(38.5,weight,ncount); 
    if (kfc==kf_Upsilon_1S)          p_obs->Insert(39.5,weight,ncount); 

    if (kfc==kf_f_2_1270)            p_obs->Insert(41.5,weight,ncount); 
    if (kfc==kf_f_2_prime_1525)      p_obs->Insert(42.5,weight,ncount); 
    if (kfc==kf_psi_2S)              p_obs->Insert(43.5,weight,ncount); 

    if (kfc==kf_p_plus)              p_obs->Insert(51.5,weight,ncount); 
    if (kfc==kf_n)                   p_obs->Insert(52.5,weight,ncount); 
    if (kfc==kf_Sigma_plus)          p_obs->Insert(53.5,weight,ncount); 
    if (kfc==kf_Sigma)               p_obs->Insert(54.5,weight,ncount); 
    if (kfc==kf_Sigma_minus)         p_obs->Insert(55.5,weight,ncount); 
    if (kfc==kf_Lambda)              p_obs->Insert(56.5,weight,ncount); 
    if (kfc==kf_Xi)                  p_obs->Insert(57.5,weight,ncount); 
    if (kfc==kf_Xi_minus)            p_obs->Insert(58.5,weight,ncount); 

    if (kfc==kf_Delta_1232_plus_plus) p_obs->Insert(61.5,weight,ncount); 
    if (kfc==kf_Delta_1232_plus)     p_obs->Insert(62.5,weight,ncount); 
    if (kfc==kf_Delta_1232)          p_obs->Insert(63.5,weight,ncount); 
    if (kfc==kf_Delta_1232_minus)    p_obs->Insert(64.5,weight,ncount); 
    if (kfc==kf_Sigma_1385_plus)     p_obs->Insert(65.5,weight,ncount); 
    if (kfc==kf_Sigma_1385)          p_obs->Insert(66.5,weight,ncount); 
    if (kfc==kf_Sigma_1385_minus)    p_obs->Insert(67.5,weight,ncount); 
    if (kfc==kf_Xi_1530)             p_obs->Insert(68.5,weight,ncount); 
    if (kfc==kf_Xi_1530_minus)       p_obs->Insert(69.5,weight,ncount); 
    if (kfc==kf_Omega_minus)         p_obs->Insert(70.5,weight,ncount); 

    if (kfc==kf_Lambda_c_plus)       p_obs->Insert(81.5,weight,ncount); 
    if (kfc==kf_Sigma_c_2455_plus)   p_obs->Insert(82.5,weight,ncount); 
    if (kfc==kf_Sigma_c_2455_plus_plus) p_obs->Insert(83.5,weight,ncount); 
    if (kfc==kf_Xi_c_2466)           p_obs->Insert(84.5,weight,ncount);     
    if (kfc==kf_Xi_c_2466_plus)      p_obs->Insert(85.5,weight,ncount);     

    if (kfc==kf_Lambda_b)            p_obs->Insert(91.5,weight,ncount); 
    if (kfc==kf_Sigma_b_5820)        p_obs->Insert(92.5,weight,ncount); 
    if (kfc==kf_Sigma_b_5820_plus)   p_obs->Insert(93.5,weight,ncount); 
    if (kfc==kf_Xi_b_5840)           p_obs->Insert(94.5,weight,ncount);     
    if (kfc==kf_Xi_b_5840_minus)     p_obs->Insert(95.5,weight,ncount);     
  }
  for (int i=0;i<100;++i) p_norm->Insert(i+.5,weight,ncount);
}

Primitive_Observable_Base * Hadron_Multiplicities::Copy() const
{
  return new Hadron_Multiplicities(1,0.,100.,100,m_listname);
}

