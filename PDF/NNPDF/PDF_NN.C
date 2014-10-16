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
    std::map<int, double> m_xfx;
    std::map<int, bool>   m_calculated;

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
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  m_type="NNPDF"; // A wild guess
  // initialise all book-keep arrays etc.
  // This is copied from LHAPDF_CPP_Interface.C
  std::vector<int> kfcs;
  kfcs.push_back(-kf_t);
  kfcs.push_back(-kf_b);
  kfcs.push_back(-kf_c);
  kfcs.push_back(-kf_s);
  kfcs.push_back(-kf_u);
  kfcs.push_back(-kf_d);
  kfcs.push_back(kf_d);
  kfcs.push_back(kf_u);
  kfcs.push_back(kf_s);
  kfcs.push_back(kf_c);
  kfcs.push_back(kf_b);
  kfcs.push_back(kf_t);
  kfcs.push_back(kf_gluon);
  for (int i=0;i<kfcs.size();i++)  {
    m_partons.insert(Flavour(abs(kfcs[i]),kfcs[i]<0));
    m_xfx[kfcs[i]]=0.;
    m_calculated[kfcs[i]]=false;
  }
  // Quark masses
  int nf(p_pdf->GetNFlavors());
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

  // Read more stuff from .info
  m_xmin=p_pdf->GetXMin(); 
  m_xmax=p_pdf->GetXMax(); 
  m_q2min=pow(p_pdf->GetQMin(),2);
  m_q2max=pow(p_pdf->GetQMax(),2);
  m_asinfo.m_order=p_pdf->GetOrderAlphaS();
  m_asinfo.m_asmz =p_pdf->GetAlphaSMz();

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

// This is resets the x and Q^2 infromation and erases all calculated values
void PDF_NNPDF::CalculateSpec(double x, double Q2)
{
  for (std::map<int,bool>::iterator it=m_calculated.begin();
       it!=m_calculated.end();++it) it->second=false;
  m_x=x/m_rescale;
  m_Q2=Q2;
}


// Return x*f(x) for flavour infl
double PDF_NNPDF::GetXPDF(const ATOOLS::Flavour infl) 
{
  int kfc = m_anti*int(infl);

  int kfc_nn(kfc);

  // nn pdf requires 0 for gluons
  if (kfc==21) kfc_nn = 0;
  if (!m_calculated[kfc]) {
    m_xfx[kfc]=p_pdf->xfx(m_x, m_Q2, kfc_nn);
    m_calculated[kfc]=true; // This is the caching bit
  }
  return m_rescale*m_xfx[kfc];
}

DECLARE_PDF_GETTER(NNPDF_Getter);


PDF_Base *NNPDF_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
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

NNPDF_Getter *p_get_nnpdf[2]; // TODO: find out number of things


extern "C" void InitPDFLib()
{
  p_get_nnpdf[0] = new NNPDF_Getter("NNPDF30NLO");
  p_get_nnpdf[1] = new NNPDF_Getter("NNPDF30NNLO");
}

extern "C" void ExitPDFLib()
{
  for (int i(0);i<2;++i) delete p_get_nnpdf[i];
}
