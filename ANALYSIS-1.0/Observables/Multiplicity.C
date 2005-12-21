#include "Multiplicity.H"

using namespace ANALYSIS;

#include "MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const String_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:"Analysed";
    return new Class(HistogramType(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),list);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list="Analysed", scale="Lin";
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

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

#include "MathTools.H"
#include "Particle_Qualifier.H"

using namespace ATOOLS;

DEFINE_OBSERVABLE_GETTER(Multiplicity,Multiplicity_Getter,"Multi")
DEFINE_OBSERVABLE_GETTER(Hadron_Multiplicities,Hadron_Multiplicities_Getter,"Hadron_Multis")




Multiplicity::Multiplicity(int type,double xmin,double xmax,int nbins,
			   const std::string & listname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins,NULL)
{
  if (listname!="") {
    m_listname = listname;
    m_name = listname+"_multi.dat";
  }
  else
    m_name = "multi.dat";
}
 
void Multiplicity::Evaluate(const ATOOLS::Particle_List & pl,
			    double weight, int ncount)
{
  p_histo->Insert(pl.size(),weight,ncount); 
}


Primitive_Observable_Base * Multiplicity::Copy() const {
  return new Multiplicity(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Hadron_Multiplicities::Hadron_Multiplicities(int type,double xmin,double xmax,int nbins,
					     const std::string & listname) :
  Primitive_Observable_Base(1,0.,100.,100,NULL)
{
  m_listname = listname; 
  m_name     = m_listname+"_Multis.dat";
  std::cout<<"Init Hadron_Multiplicities: listname = |"<<m_listname<<"|"<<std::endl;
}

void Hadron_Multiplicities::Evaluate(const ATOOLS::Particle_List & pl, double weight, int ncount)
{
  //Particle_List * pl = p_ana->GetParticleList(m_listname);
  Flavour flav;
  kf::code kfc;

  for (Particle_List::const_iterator pliter=pl.begin();pliter!=pl.end();pliter++) {
    flav = (*pliter)->Flav();
    if (!flav.IsHadron() && !flav.IsPhoton()) continue;
    kfc  = flav.Kfcode();
    if (kfc==kf::photon)              p_histo->Insert(0,weight,ncount);
    if (kfc==kf::pi)                  p_histo->Insert(1,weight,ncount);
    if (kfc==kf::pi_plus)             p_histo->Insert(2,weight,ncount);
    if (kfc==kf::eta)                 p_histo->Insert(3,weight,ncount);
    if (m_listname=="PrimordialHadrons") {
      if (kfc==kf::K)                 p_histo->Insert(4,weight,ncount);
    }
    else {
      if (kfc==kf::K_L||kfc==kf::K_S) p_histo->Insert(4,weight,ncount);
    }
    if (kfc==kf::K_plus)              p_histo->Insert(5,weight,ncount);
    if (kfc==kf::eta_prime_958)       p_histo->Insert(6,weight,ncount);

    if (kfc==kf::D_plus)              p_histo->Insert(11,weight,ncount);
    if (kfc==kf::D)                   p_histo->Insert(12,weight,ncount);
    if (kfc==kf::D_s_plus)            p_histo->Insert(13,weight,ncount);

    if (kfc==kf::rho_770)             p_histo->Insert(21,weight,ncount);
    if (kfc==kf::rho_770_plus)        p_histo->Insert(22,weight,ncount);
    if (kfc==kf::omega_782)           p_histo->Insert(23,weight,ncount);
    if (kfc==kf::K_star_892)          p_histo->Insert(24,weight,ncount);
    if (kfc==kf::K_star_892_plus)     p_histo->Insert(25,weight,ncount);
    if (kfc==kf::phi_1020)            p_histo->Insert(26,weight,ncount);

    if (kfc==kf::D_star_2010_plus)    p_histo->Insert(31,weight,ncount);
    if (kfc==kf::D_star_2007)         p_histo->Insert(32,weight,ncount);
    if (kfc==kf::D_s_star_plus)       p_histo->Insert(33,weight,ncount);

    if (kfc==kf::p_plus)                        p_histo->Insert(51,weight,ncount);
    if (kfc==kf::n)                             p_histo->Insert(52,weight,ncount);
    if (kfc==kf::Sigma_plus)                    p_histo->Insert(53,weight,ncount);
    if (kfc==kf::Sigma)                         p_histo->Insert(54,weight,ncount);
    if (kfc==kf::Sigma_minus)                   p_histo->Insert(55,weight,ncount);
    if (kfc==kf::Lambda)                        p_histo->Insert(56,weight,ncount);
    if (kfc==kf::Xi)                            p_histo->Insert(57,weight,ncount);
    if (kfc==kf::Xi_minus)                      p_histo->Insert(58,weight,ncount);

    if (kfc==kf::Delta_1232_plus_plus)          p_histo->Insert(61,weight,ncount);
    if (kfc==kf::Delta_1232_plus)               p_histo->Insert(62,weight,ncount);
    if (kfc==kf::Delta_1232)                    p_histo->Insert(63,weight,ncount);
    if (kfc==kf::Delta_1232_minus)              p_histo->Insert(64,weight,ncount);
    if (kfc==kf::Sigma_1385_plus)               p_histo->Insert(65,weight,ncount);
    if (kfc==kf::Sigma_1385)                    p_histo->Insert(66,weight,ncount);
    if (kfc==kf::Sigma_1385_minus)              p_histo->Insert(67,weight,ncount);
    if (kfc==kf::Xi_1530)                       p_histo->Insert(68,weight,ncount);
    if (kfc==kf::Xi_1530_minus)                 p_histo->Insert(69,weight,ncount);
    if (kfc==kf::Omega_minus)                   p_histo->Insert(70,weight,ncount);
  }
}

Primitive_Observable_Base * Hadron_Multiplicities::Copy() const
{
  return new Hadron_Multiplicities(1,0.,100.,100,m_listname);
}

