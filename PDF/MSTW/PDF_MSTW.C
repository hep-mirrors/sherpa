#ifndef PDF_MSTW_PDF_MSTW_H
#define PDF_MSTW_PDF_MSTW_H

#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "mstwpdf.h"

namespace PDF {

  class PDF_MSTW : public PDF_Base {
  private:

    c_mstwpdf *p_pdf;

    std::string m_path, m_file;

    int    m_anti;
    double m_x, m_Q2;

  public:

    PDF_MSTW(const ATOOLS::Flavour &bunch,const std::string &file,int set); 

    ~PDF_MSTW(); 

    PDF_Base * GetCopy();

    void   CalculateSpec(double,double);
    double GetXPDF(const ATOOLS::Flavour); 

  };// end of class PDF_MSTW

}  

#endif

#include "mstwpdf.h"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

using namespace PDF;
using namespace ATOOLS;

PDF_MSTW::PDF_MSTW
(const ATOOLS::Flavour &bunch,
 const std::string &bfile,int set):
  m_path(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/MSTW08Grid"),
  m_file(bfile), m_anti(1)
{
  m_member=set;
  m_type="MSTW";
  m_bunch=bunch;
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  for (int i=1;i<6;i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());
  m_partons.insert(Flavour(kf_photon));
  std::string file(m_file);
  if (m_member<100) file+=".00.dat";
  else if (m_member<200)
    file+=(m_member<110?".68cl.0":".68cl.")+ToString(m_member-100)+".dat";
  else file+=(m_member<210?".90cl.0":".90cl.")+ToString(m_member-200)+".dat";
  if (m_member>100) msg_Info()<<METHOD<<"(): Init member "<<m_member
			      <<", file '"<<file<<"'."<<std::endl;
  p_pdf = new c_mstwpdf(m_path+"/"+file);
  m_xmin=p_pdf->xmin;
  m_xmax=p_pdf->xmax;
  m_q2min=p_pdf->qsqmin;
  m_q2max=p_pdf->qsqmax;
  m_asinfo.m_order=p_pdf->alphaSorder;
  m_asinfo.m_asmz=p_pdf->alphaSMZ;
}

PDF_MSTW::~PDF_MSTW()
{
  delete p_pdf;
}

PDF_Base *PDF_MSTW::GetCopy() 
{
  PDF_Base *copy = new PDF_MSTW(m_bunch,m_file,m_member);
  m_copies.push_back(copy);
  return copy;
}

void PDF_MSTW::CalculateSpec(double x,double Q2)
{
  m_x=x;
  m_Q2=Q2;
}


double PDF_MSTW::GetXPDF(const ATOOLS::Flavour infl) 
{
  if(m_x<m_xmin) m_x=m_xmin;
  if (m_x/m_rescale>m_xmax || m_rescale<0.0) return 0.0;
  int kfc=m_anti*int(infl);
  if (abs(kfc)==kf_gluon) kfc=0;
  else if (abs(kfc)==kf_photon) kfc=13;
  return m_rescale*p_pdf->parton(kfc,m_x/m_rescale,sqrt(m_Q2));
}

DECLARE_PDF_GETTER(MSTW_Getter);

PDF_Base *MSTW_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  int set=args.p_read->GetValue<int>("PDF_SET_VERSION",1);
  int ibeam=args.m_ibeam;
  set=args.p_read->GetValue<int>("PDF_SET_VERSION_"+ToString(ibeam+1),set);
  return new PDF_MSTW(args.m_bunch,m_key,set);
}

void MSTW_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"MSTW 2008 fit, see arXiv:0901:0002 [hep-ph]";
}

MSTW_Getter *p_get_mstw08[63];

extern "C" void InitPDFLib()
{
  p_get_mstw08[0] = new MSTW_Getter("mstw2008lo");
  p_get_mstw08[1] = new MSTW_Getter("mstw2008nlo");
  p_get_mstw08[2] = new MSTW_Getter("mstw2008nnlo");
  for (int i(10);i<=30;++i) 
    p_get_mstw08[i-7] = new MSTW_Getter("mstw2008nlo_asmz1"+ToString(i));
  for (int i(7);i<=9;++i) 
    p_get_mstw08[i+17] = new MSTW_Getter("mstw2008nnlo_asmz10"+ToString(i));
  for (int i(10);i<=27;++i) 
    p_get_mstw08[i+17] = new MSTW_Getter("mstw2008nnlo_asmz1"+ToString(i));
  p_get_mstw08[45] = new MSTW_Getter("mstw2008lo_nf3");
  p_get_mstw08[46] = new MSTW_Getter("mstw2008lo_nf4");
  p_get_mstw08[47] = new MSTW_Getter("mstw2008nlo_nf3");
  p_get_mstw08[48] = new MSTW_Getter("mstw2008nlo_nf4");
  p_get_mstw08[49] = new MSTW_Getter("mstw2008nnlo_nf3");
  p_get_mstw08[50] = new MSTW_Getter("mstw2008nnlo_nf4");
  p_get_mstw08[51] = new MSTW_Getter("mstw2008nlo_mb4.00GeV.00");
  p_get_mstw08[52] = new MSTW_Getter("mstw2008nlo_mb4.00GeV_nf4.00");
  p_get_mstw08[53] = new MSTW_Getter("mstw2008nlo_mb4.25GeV.00");
  p_get_mstw08[54] = new MSTW_Getter("mstw2008nlo_mb4.25GeV_nf4.00");
  p_get_mstw08[55] = new MSTW_Getter("mstw2008nlo_mb4.50GeV.00");
  p_get_mstw08[56] = new MSTW_Getter("mstw2008nlo_mb4.50GeV_nf4.00");
  p_get_mstw08[57] = new MSTW_Getter("mstw2008nlo_mb5.00GeV.00");
  p_get_mstw08[58] = new MSTW_Getter("mstw2008nlo_mb5.00GeV_nf4.00");
  p_get_mstw08[59] = new MSTW_Getter("mstw2008nlo_mb5.25GeV.00");
  p_get_mstw08[60] = new MSTW_Getter("mstw2008nlo_mb5.25GeV_nf4.00");
  p_get_mstw08[61] = new MSTW_Getter("mstw2008nlo_mb5.50GeV.00");
  p_get_mstw08[62] = new MSTW_Getter("mstw2008nlo_mb5.50GeV_nf4.00");
}

extern "C" void ExitPDFLib()
{
  for (int i(0);i<63;++i) delete p_get_mstw08[i];
}
