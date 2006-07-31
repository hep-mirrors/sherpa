#include "Hadron_Multiplicity.H"
#include "Message.H"
#include "MathTools.H"
#include "Particle_Qualifier.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

#define DEFINE_GETTER_METHOD(CLASS,NAME)                                \
  Primitive_Observable_Base *                                   \
  NAME::operator()(const String_Matrix &parameters) const               \
{ return new CLASS(parameters); }

#define DEFINE_PRINT_METHOD(NAME)                                       \
  void NAME::PrintInfo(ostream &str,const size_t width) const   \
{ str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)                        \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);     \
  DEFINE_GETTER_METHOD(CLASS,NAME)                                      \
  DEFINE_PRINT_METHOD(NAME)

DEFINE_OBSERVABLE_GETTER(Hadron_Multiplicity,Hadron_Multiplicity_Getter,"Hadron_Multi")

Hadron_Multiplicity::Hadron_Multiplicity(const String_Matrix & parameters):
    Primitive_Observable_Base(parameters)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable Hadron_Multi in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"IntermediateHadrons";
    m_name="";
    if (m_listname!="IntermediateHadrons") m_name=m_listname+string("_");
    m_name+="Hadron_Multi.dat";
  }
  else {
    if (m_listname=="") m_listname="IntermediateHadrons";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="IntermediateHadrons") m_name=m_listname+string("_");
      m_name+="Hadron_Multi";
    }
    if (p_histo->Title()=="SherpaDefault") {
      string title = "Hadron_Multiplicity";
      p_histo->SetTitle(title);
    }
  }
}

Hadron_Multiplicity::Hadron_Multiplicity(const Hadron_Multiplicity * old) :
  Primitive_Observable_Base(*old)
{
  m_listname = old->m_listname; // ?
  m_name     = old->m_name;     // ?
}

Primitive_Observable_Base * Hadron_Multiplicity::Copy() const 
{ return new Hadron_Multiplicity(this); }

void Hadron_Multiplicity::Evaluate(const Particle_List & pl,
				   double weight, int ncount)
{
  Flavour flav;
  kf::code kfc;

  for (Particle_List::const_iterator pliter=pl.begin();pliter!=pl.end();pliter++) {
    flav = (*pliter)->Flav();
    if (!flav.IsHadron() && !flav.IsPhoton()) continue;
    kfc  = flav.Kfcode();
    if (kfc==kf::photon)              p_histo->Insert("#gamma",weight,ncount);

    if (kfc==kf::pi)                  p_histo->Insert("#pi^{0}",weight,ncount);
    if (kfc==kf::pi_plus)             p_histo->Insert("#pi^{#pm}",weight,ncount);
    if (kfc==kf::eta)                 p_histo->Insert("#eta",weight,ncount);
    if (m_listname=="PrimordialHadrons") {
      if (kfc==kf::K)                 p_histo->Insert("K^{0}",weight,ncount);
    }
    else {
      if (kfc==kf::K_L||kfc==kf::K_S) p_histo->Insert("K^{0}",weight,ncount);
    }
    if (kfc==kf::K_plus)              p_histo->Insert("K^{#pm}",weight,ncount);
    if (kfc==kf::eta_prime_958)       p_histo->Insert("#eta'",weight,ncount);
    if (kfc==kf::D_plus)              p_histo->Insert("D^{#pm}",weight,ncount);
    if (kfc==kf::D)                   p_histo->Insert("D^{0}",weight,ncount);
    if (kfc==kf::D_s_plus)            p_histo->Insert("D_{s}",weight,ncount);
    if (kfc==kf::B_plus)              p_histo->Insert("B^{#pm}",weight,ncount);
    if (kfc==kf::B)                   p_histo->Insert("B^{0}",weight,ncount);
    if (kfc==kf::B_s)                 p_histo->Insert("B_{s}",weight,ncount);

    if (kfc==kf::rho_770)             p_histo->Insert("#rho^{0}",weight,ncount);
    if (kfc==kf::rho_770_plus)        p_histo->Insert("#rho^{#pm}",weight,ncount);
    if (kfc==kf::omega_782)           p_histo->Insert("#omega",weight,ncount);
    if (kfc==kf::K_star_892)          p_histo->Insert("K^{*0}",weight,ncount);
    if (kfc==kf::K_star_892_plus)     p_histo->Insert("K^{*#pm}",weight,ncount);
    if (kfc==kf::phi_1020)            p_histo->Insert("#phi",weight,ncount);
    if (kfc==kf::D_star_2010_plus)    p_histo->Insert("D^{*#pm}",weight,ncount);
    if (kfc==kf::D_star_2007)         p_histo->Insert("D^{*0}",weight,ncount);
    if (kfc==kf::D_s_star_plus)       p_histo->Insert("D_{s}^{*}",weight,ncount);
    if (kfc==kf::B_star_plus)         p_histo->Insert("D^{*#pm}",weight,ncount);
    if (kfc==kf::B_star)              p_histo->Insert("D^{*0}",weight,ncount);
    if (kfc==kf::B_s_star)            p_histo->Insert("D_{s}^{*}",weight,ncount);

    if (kfc==kf::J_psi_1S)            p_histo->Insert("J/#Psi",weight,ncount);
    if (kfc==kf::psi_2S)              p_histo->Insert("#Psi'(3685)",weight,ncount);

    if (kfc==kf::p_plus)                        p_histo->Insert("p",weight,ncount);
    if (kfc==kf::n)                             p_histo->Insert("n",weight,ncount);
    if (kfc==kf::Sigma_plus)                    p_histo->Insert("#Sigma^{+}",weight,ncount);
    if (kfc==kf::Sigma)                         p_histo->Insert("#Sigma^{0}",weight,ncount);
    if (kfc==kf::Sigma_minus)                   p_histo->Insert("#Sigma^{-}",weight,ncount);
    if (kfc==kf::Lambda)                        p_histo->Insert("#Lambda",weight,ncount);
    if (kfc==kf::Xi)                            p_histo->Insert("#Xi^{0}",weight,ncount);
    if (kfc==kf::Xi_minus)                      p_histo->Insert("#Xi^{-}",weight,ncount);

    if (kfc==kf::Delta_1232_plus_plus)          p_histo->Insert("#Delta^{++}",weight,ncount);
    if (kfc==kf::Delta_1232_plus)               p_histo->Insert("#Delta^{+}",weight,ncount);
    if (kfc==kf::Delta_1232)                    p_histo->Insert("#Delta^{0}",weight,ncount);
    if (kfc==kf::Delta_1232_minus)              p_histo->Insert("#Delta^{-}",weight,ncount);
    if (kfc==kf::Sigma_1385_plus)               p_histo->Insert("#Sigma(1385)^{#pm}",weight,ncount);
    if (kfc==kf::Sigma_1385)                    p_histo->Insert("#Sigma(1385)^{0}",weight,ncount);
    if (kfc==kf::Sigma_1385_minus)              p_histo->Insert("#Sigma(1385)^{#pm}",weight,ncount);
    if (kfc==kf::Xi_1530)                       p_histo->Insert("#Xi(1530)^{0}",weight,ncount); 
    if (kfc==kf::Xi_1530_minus)                 p_histo->Insert("#Xi(1530)^{-}",weight,ncount);
    if (kfc==kf::Omega_minus)                   p_histo->Insert("#Omega^{-}",weight,ncount);
  }
}





