#ifndef PDF_NNPDF_PDF_NNPDF_H
#define PDF_NNPDF_PDF_NNPDF_H

#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "NNPDFDriver.h"

// This is all copied from the MSTW code
namespace PDF {

  class PDF_NNPDF : public PDF_Base {
  private:

    // Use the driver supplied by NNPDF  
    NNPDFDriver *p_pdf;

    std::string m_path, m_file;

    int    m_anti;
    int m_lookup[28];
    int m_prefix;


    double m_x, m_Q2;
    std::map<int, double> m_xfx;
    std::map<int, bool>   m_calculated;

  public:

    PDF_NNPDF(const ATOOLS::Flavour &bunch,const std::string &file,
              const std::string &set,int member, int prefix);

    ~PDF_NNPDF(); 

    PDF_Base * GetCopy();

    void   CalculateSpec(const double&,const double&);
    double GetXPDF(const ATOOLS::Flavour&);
    double GetXPDF(const kf_code&,bool);

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
 const std::string &bfile,
 const std::string &set,int member, int prefix):
  m_path(rpa->gen.Variable("SHERPA_SHARE_PATH")),
  m_file(bfile), m_anti(1)
{
  m_set=set;
  m_member=member;
  m_prefix=prefix;
  m_lhef_number=prefix+member;
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
    //m_calculated[kfcs[i]]=false;
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
  m_asinfo.m_asmz=p_pdf->GetAlphaSMz();
  m_asinfo.m_mz2=sqr(p_pdf->GetMz());

  // Hopefully efficient lookup table
  // 0 entspricht -6, also anti top
  // 1 entspricht -5, also anti bottom
  // ...
  // 12 entspricht 6, also top
  // LHAPDF5 style: 6 entspraeche 0, also gluon, aber sherpa nimmt gluons als 21
  for (int i=0;i<13;++i) {
    m_lookup[i] = i;

  }
  for (int i=13;i<27;++i) {
    m_lookup[i] = -1; // A safety, NNPDFdriver expects a number between 0 and 12
  }
  m_lookup[27] = 6; // The gluon in LHAPDF5/NNPDFDriver convention

}

PDF_NNPDF::~PDF_NNPDF()
{
  delete p_pdf;
}

// Necessary?
PDF_Base *PDF_NNPDF::GetCopy() 
{
  PDF_Base *copy = new PDF_NNPDF(m_bunch,m_file,m_set,m_member,m_prefix);
  m_copies.push_back(copy);
  return copy;
}

// This is resets the x and Q^2 infromation and erases all calculated values
void PDF_NNPDF::CalculateSpec(const double& x, const double& Q2)
{
  //for (std::map<int,bool>::iterator it=m_calculated.begin();
       //it!=m_calculated.end();++it) it->second=false;
  m_x=x/m_rescale;
  m_Q2=Q2;
}


// Return x*f(x) for flavour infl
double PDF_NNPDF::GetXPDF(const ATOOLS::Flavour& infl)
{
  // Parton flavour IDs
  int kfc = (infl == 21 || infl == 22) ? int(infl) : m_anti*int(infl);
  // Hopefully efficient lookup --- relate 21 to 0
  int kfc_nn(m_lookup[kfc+6]); // kfc runs from -6 to 6 and also 21
                               // While the driver wants 
                               // numbers from 0 to 12 to access
                               // array elements
 
  return m_rescale*p_pdf->xfx(m_x, m_Q2, kfc_nn);
}

double PDF_NNPDF::GetXPDF(const kf_code& kf, bool anti)
{
  // Parton flavour IDs
  int kfc = (kf == 21 || kf == 22) ? kf : m_anti*(anti?-kf:kf);
  // Hopefully efficient lookup --- relate 21 to 0
  int kfc_nn(m_lookup[kfc+6]); // kfc runs from -6 to 6 and also 21
                               // While the driver wants
                               // numbers from 0 to 12 to access
                               // array elements

  return m_rescale*p_pdf->xfx(m_x, m_Q2, kfc_nn);
}

DECLARE_PDF_GETTER(NNPDF_Getter);


PDF_Base *NNPDF_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  std::string gfile;
  int pdfsetprefix=-1;
  if (args.m_set == "NNPDF30NLO") {
    gfile = std::string("NNPDF30_nlo_as_0118");
    pdfsetprefix=260000;
    if (args.m_member>100 || args.m_member <0) {
      THROW(fatal_error,"PDF_SET_MEMBER out of range [0,100].");
    }
  }
  else if (args.m_set == "NNPDF30NNLO") {
    gfile = std::string("NNPDF30_nnlo_as_0118");
    pdfsetprefix=261000;
    if (args.m_member>100 || args.m_member <0) {
      THROW(fatal_error,"PDF_SET_MEMBER out of range [0,100].");
    }
  }
  else THROW(not_implemented,"Requested PDF_SET not available.");
  return new PDF_NNPDF(args.m_bunch, gfile, args.m_set, args.m_member, pdfsetprefix);
}

void NNPDF_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"NNPDF fit, see arXiv:1410.8849 [hep-ph]";
}

NNPDF_Getter *p_get_nnpdf[2];


extern "C" void InitPDFLib()
{
  p_get_nnpdf[0] = new NNPDF_Getter("NNPDF30NLO");
  p_get_nnpdf[1] = new NNPDF_Getter("NNPDF30NNLO");
}

extern "C" void ExitPDFLib()
{
  for (int i(0);i<2;++i) delete p_get_nnpdf[i];
}
