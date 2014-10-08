#ifndef PDF_NNPDF_PDF_NNPDF_H
#define PDF_NNPDF_PDF_NNPDF_H

#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "NNPDFDriver.h"

// This is all copied from the MSTW code
namespace PDF {

  class PDF_NNPDF : public PDF_Base {
  private:

    // Use the driver supplied by NNPDF  
    NNPDFDriver *p_pdf;

    std::string m_path, m_file;

    int    m_anti;
    double m_x, m_Q2;

  public:

    PDF_NNPDF(const ATOOLS::Flavour &bunch,const std::string &file,int set); 

    ~PDF_NNPDF(); 

    PDF_Base * GetCopy();

    void   CalculateSpec(double,double);
    double GetXPDF(const ATOOLS::Flavour); 

  };// end of class PDF_NNPDF

}  

#endif

#include "NNPDFDriver.h"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

using namespace PDF;
using namespace ATOOLS;

PDF_NNPDF::PDF_NNPDF
(const ATOOLS::Flavour &bunch,
 const std::string &bfile,int set):
  m_path(rpa->gen.Variable("SHERPA_SHARE_PATH")),
  m_file(bfile), m_anti(1)
{
  m_member=set;
  std::string file(m_file);
  p_pdf = new NNPDFDriver(m_path+"/"+file, m_member); // Path to the file to load
  
  m_bunch=bunch; // This is the beam
  m_type="NNPDF"; // A wild guess
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  for (int i=1;i<6;i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());
  // Use NNPDFdriver's hasPhoton method to decide whether a set has a photon PDF
  if (p_pdf->hasPhoton()) m_partons.insert(Flavour(kf_photon));


  // Read more stuff from .info
  m_xmin=p_pdf->GetXMin(); 
  m_xmax=p_pdf->GetXMax(); 
  m_q2min=pow(p_pdf->GetQMin(),2);
  m_q2max=pow(p_pdf->GetQMax(),2);
  m_asinfo.m_order=p_pdf->GetOrderAlphaS();
  m_asinfo.m_asmz =p_pdf->GetAlphaSMz();
  
  // Quark masses
  int nf(p_pdf->GetNFL());
  //int nf(p_pdf->info().get_entry_as<int>("NumFlavors"));
  if (nf<0) m_asinfo.m_flavs.resize(5);
  else      m_asinfo.m_flavs.resize(nf);
  // for now assume thresholds are equal to masses, as does LHAPDF-6.0.0
  for (size_t i(0);i<m_asinfo.m_flavs.size();++i) {
    m_asinfo.m_flavs[i]=PDF_Flavour((kf_code)i+1);
    if      (i==0)
      m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
          =p_pdf->GetMDown();
    else if (i==1)
      m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
          =p_pdf->GetMUp();
    else if (i==2)
      m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
          =p_pdf->GetMStrange();
    else if (i==3)
      m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
          =p_pdf->GetMCharm();
    else if (i==4)
      m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
          =p_pdf->GetMBottom();
    else if (i==5)
      m_asinfo.m_flavs[i].m_mass=m_asinfo.m_flavs[i].m_thres
          =p_pdf->GetMTop();
  }

}

PDF_NNPDF::~PDF_NNPDF()
{
  delete p_pdf;
}

// Necessary?
PDF_Base *PDF_NNPDF::GetCopy() 
{
  PDF_Base *copy = new PDF_NNPDF(m_bunch,m_file,m_member);
  m_copies.push_back(copy);
  return copy;
}

// More like set, but ok - TODO: some way to to the m_calculated bit
void PDF_NNPDF::CalculateSpec(double x,double Q2)
{
  // TODO: Marek wants some more efficiency here using m_calculated
  m_x=x;
  m_Q2=Q2;
}


// Return x*f(x) for flavour infl
double PDF_NNPDF::GetXPDF(const ATOOLS::Flavour infl) 
{
  // TODO: isn't this implemented in the driver?
  //if(m_x<m_xmin) m_x=m_xmin;
  //if (m_x/m_rescale>m_xmax || m_rescale<0.0) return 0.0;
 
  // Some bizarre logic about flavours  
  int kfc=m_anti*int(infl);
  if (abs(kfc)==kf_gluon) kfc=0;
  else if (abs(kfc)==kf_photon) kfc=13;
  // Get xfx from NNPDFdriver, apply m_rescale (TODO, where is m_rescale set?)
  // TODO: get rid of sqrt in the following?
  return m_rescale*p_pdf->xfx(m_x/m_rescale, sqrt(m_Q2), kfc);
}

DECLARE_PDF_GETTER(NNPDF_Getter);


PDF_Base *NNPDF_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  //int set=args.p_read->GetValue<int>("PDF_SET_VERSION",1); // 1 is the default value I think this is not needed
  int ibeam=args.m_ibeam;
  int member = args.p_read->GetValue<int>("PDF_SET_MEMBER", 0); // 0 is the default value
  std::string gfile;
  if (m_key == "NNPDF30NLO") {
    gfile = std::string("NNPDF30_nlo_as_0118");
    if (member>100 || member <0) {
      std::cerr << "Error, PDF_SET_MEMBER not in [0,100] --- exiting!" << std::endl;
      exit(1);
    }
  }
  else if (m_key == "NNPDF30NNLO") {
    gfile = std::string("NNPDF30_nnlo_as_0118");
    if (member>100 || member <0) {
      std::cerr << "Error, PDF_SET_MEMBER not in [0,100] --- exiting!" << std::endl;
      exit(1);
    }
  }
  else {
    std::cerr << "Requested PDF_SET " << m_key << " not available --- exiting!" << std::endl;
    std::cerr << "Run Sherpa with SHOW_PDF_SETS=1" << std::endl;
  }
  return new PDF_NNPDF(args.m_bunch, gfile, member);
}

void NNPDF_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"NNPDF fit, see arXiv"; // TODO: change to proper reference
}

NNPDF_Getter *p_get_nnpdf[2]; // TODO: wtf? // TODO: find out number of things


extern "C" void InitPDFLib()
{
  p_get_nnpdf[0] = new NNPDF_Getter("NNPDF30NLO");
  p_get_nnpdf[1] = new NNPDF_Getter("NNPDF30NNLO");
}

extern "C" void ExitPDFLib()
{
  for (int i(0);i<2;++i) delete p_get_nnpdf[i];
}
