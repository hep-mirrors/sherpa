#include "SHERPA/Tools/Scale_Variations.H"

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/PDF_Base.H"
#include "PHASIC++/Process/Process_Base.H"

#if defined USING__LHAPDF && defined USING__LHAPDF6
#include "LHAPDF/LHAPDF.h"
#endif

namespace SHERPA {
  std::ostream& operator<<(std::ostream &s,const Scale_Variation &sv)
  {
    return s<<"Scale_Variation[mu_{R/F}^2-facs=("
            <<sv.MuR2Fac()<<","<<sv.MuF2Fac()
            <<"),PDF="<<sv.PdfId()<<",val="<<sv.Value()<<"]";
  }

  std::ostream& operator<<(std::ostream &s,const NamedScaleVariationMap &nsvm)
  {
    s<<"Named scale variations:"<<std::endl;
    for (NamedScaleVariationMap::const_iterator it=nsvm.begin();
         it!=nsvm.end();++it) s<<it->first<<" : "<<*it->second<<std::endl;
    return s;
  }

  std::ostream& operator<<(std::ostream &s,const NamedScaleVariationMap *nsvm)
  {
    return s<<nsvm->size()<<" variations";
  }

  std::ostream& operator<<(std::ostream &s,const Scale_Variations &svs)
  {
    return s<<svs.GetNamedScalesMap();
  }

}

using namespace ATOOLS;
using namespace SHERPA;

Scale_Variation::Scale_Variation(const double &muR2fac, const double &muF2fac,
                                 PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
                                 bool deletepdfs) :
  m_deletepdfs(deletepdfs),
  m_muR2fac(muR2fac), m_muF2fac(muF2fac), m_val(0.),
  p_pdf1(pdf1), p_pdf2(pdf2),
  m_pdf1id(pdf1->LHEFNumber()), m_pdf2id(pdf2->LHEFNumber()),
  m_pdf1set(pdf1->Set()), m_pdf2set(pdf2->Set()),
  m_pdf1setmember(pdf1->Member()), m_pdf2setmember(pdf2->Member()),
  m_name(GenerateName())
{
}

Scale_Variation::~Scale_Variation()
{
  if (m_deletepdfs) {
    if (p_pdf1) { delete p_pdf1; p_pdf1=NULL; }
    if (p_pdf2) { delete p_pdf2; p_pdf1=NULL; }
  }
}

std::string Scale_Variation::GenerateName()
{
  if (m_pdf1id==m_pdf2id)
    return std::string("MUR")+ATOOLS::ToString(sqrt(m_muR2fac))+std::string("_")
           +std::string("MUF")+ATOOLS::ToString(sqrt(m_muF2fac))+std::string("_")
           +std::string("PDF")+ATOOLS::ToString(m_pdf1id);
  else
    return std::string("MUR")+ATOOLS::ToString(sqrt(m_muR2fac))+std::string("_")
           +std::string("MUF")+ATOOLS::ToString(sqrt(m_muF2fac))+std::string("_")
           +std::string("PDF")+ATOOLS::ToString(m_pdf1id)+std::string("_")
           +std::string("PDF")+ATOOLS::ToString(m_pdf2id);
}

Scale_Variations::Scale_Variations() :
  m_on(false), m_loadlhapdf(true), p_nsvmap(new NamedScaleVariationMap())
{
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.AddWordSeparator("\t");
  std::vector<std::string> vars,pdfs;
  reader.VectorFromFile(vars,"SCALE_VARIATIONS");
  reader.VectorFromFile(pdfs,"PDF_VARIATIONS");
  for (size_t i(0);i<pdfs.size();++i) vars.push_back("1.,1.,"+pdfs[i]);
  if (vars.size()) m_on=true;
  if (!m_on) return;
  PRINT_FUNC(vars.size());
#if defined USING__LHAPDF && defined USING__LHAPDF6
  // assume that LHAPDF is already loaded and interface initialised
  std::string path;
  if (reader.ReadFromFile(path,"LHAPDF_GRID_PATH")) LHAPDF::setPaths(path);
  const std::vector<std::string>& avsets(LHAPDF::availablePDFSets());
#endif
  for (size_t i(0);i<vars.size();++i) {
    std::string cur(vars[i]);
    msg_Debugging()<<cur<<std::endl;
    size_t pos1(cur.find(",",0));
    double muR2fac(ToType<double>(std::string(cur,0,pos1)));
    size_t pos2(cur.find(",",pos1+1));
    double muF2fac(ToType<double>(std::string(cur,pos1+1,pos2)));
    if (pos2==std::string::npos) {
      // only use nominal PDF
      Scale_Variation * scvar(new Scale_Variation(muR2fac,muF2fac,
                                                  rpa->gen.PDF(0),
                                                  rpa->gen.PDF(1),false));
      (*p_nsvmap)[scvar->Name()]=scvar;
    }
    else {
      // switch to other PDFs, only works with LHAPDF6
      // only then can initialise full set of error PDFs at once,
      // otherwise need to be given explicitly
      std::string pdfname(cur,pos2+1,std::string::npos);
      bool fullset(false);
      if (pdfname.find("[all]")!=std::string::npos) fullset=true;
      PDF::PDF_Base *pdf1(NULL),*pdf2(NULL);
      if (!fullset) {
        int member(0);
        if (pdfname.find("/")!=std::string::npos) {
          member=ToType<int>(std::string(pdfname,pdfname.find("/")+1));
          pdfname=pdfname.substr(0,pdfname.find("/"));
        }
        PDF::PDF_Arguments pa1(rpa->gen.Bunch(0),&reader,0,pdfname,member);
        PDF::PDF_Arguments pa2(rpa->gen.Bunch(1),&reader,1,pdfname,member);
        pdf1=PDF::PDF_Base::PDF_Getter_Function::GetObject(pdfname,pa1);
        pdf2=PDF::PDF_Base::PDF_Getter_Function::GetObject(pdfname,pa2);
        if (!pdf1 || !pdf2)
          THROW(fatal_error,"PDF set "+pdfname+" not available.");
        Scale_Variation * scvar(new Scale_Variation(muR2fac,muF2fac,
                                                    pdf1,pdf2,true));
        (*p_nsvmap)[scvar->Name()]=scvar;
      }
      else {
#if defined USING__LHAPDF && defined USING__LHAPDF6
        // check whether interface is loaded
        if (!ATOOLS::s_loader->LibraryIsLoaded("LHAPDFSherpa"))
          THROW(fatal_error,"LHAPDF interface not initialised. "+std::string("")
                            +"Add LHAPDFSherpa to PDF_LIBRARY");
        // assume members are labeled 0..n
        pdfname=pdfname.substr(0,pdfname.find("[all]"));
        bool fnd(false);
        for (size_t j(0);j<avsets.size();++j) if (avsets[j]==pdfname) fnd=true;
        if (!fnd) THROW(fatal_error,"PDF set "+pdfname+" not available.");
        LHAPDF::PDFSet set(pdfname);
        for (size_t j(0);j<set.size();++j) {
          PDF::PDF_Arguments pa1(rpa->gen.Bunch(0),&reader,0,pdfname,j);
          PDF::PDF_Arguments pa2(rpa->gen.Bunch(1),&reader,1,pdfname,j);
          pdf1=PDF::PDF_Base::PDF_Getter_Function::GetObject(pdfname,pa1);
          pdf2=PDF::PDF_Base::PDF_Getter_Function::GetObject(pdfname,pa2);
          Scale_Variation * scvar(new Scale_Variation(muR2fac,muF2fac,
                                                      pdf1,pdf2,true));
          (*p_nsvmap)[scvar->Name()]=scvar;
        }
#else
        THROW(not_implemented,"Full set reweightings only work with LHAPDF6."
                              +std::string(" Otherwise specify separately."));
#endif
      }
    }
  }
  msg_Debugging()<<*p_nsvmap<<std::endl;
}

Scale_Variations::~Scale_Variations()
{
  if (p_nsvmap) {
    for (NamedScaleVariationMap::iterator it=p_nsvmap->begin();
         it!=p_nsvmap->end();++it) {
      delete it->second;
    }
    delete p_nsvmap;
  }
}

void Scale_Variations::ResetValues()
{
  for (NamedScaleVariationMap::iterator it=p_nsvmap->begin();
       it!=p_nsvmap->end();++it) it->second->SetValue(0.);
}

void Scale_Variations::ExtractParameters(const ATOOLS::Weight_Info &winfo,
                                         PHASIC::Process_Base * proc)
{
  const ME_Weight_Info * const mewgt(proc->GetMEwgtinfo());
  m_params.oqcd=proc->OrderQCD();
  m_params.oew=proc->OrderEW();
  m_params.x1=winfo.m_pdf.m_x1;
  m_params.x2=winfo.m_pdf.m_x2;
  m_params.x1p=mewgt->m_y1;
  m_params.x2p=mewgt->m_y2;
}

bool Scale_Variations::Calculate(Scale_Variation * sv)
{
  return true;
}

double Scale_Variations::Calculate
(const double& B, const double& V,
 const double& I, const std::vector<double>& wgts,
 const Flavour& fl1, const Flavour& fl2,
 PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
 const double& x1, const double& x2,
 const double& x1p, const double& x2p,
 const double& muR2,
 const double& muF12, const double& muF22,
 const double& muR2fac, const double& muF2fac,
 const double& oqcd, const double& oew,
 const size_t& mode)
{
  DEBUG_FUNC("mode="<<mode);
  // calculate new event weight
  double muR2new(muR2*muR2fac);
  double muF12new(muF12*muF2fac),muF22new(muF22*muF2fac);
  pdf1->Calculate(x1,muF12new);
  pdf2->Calculate(x2,muF22new);
  double fa=pdf1->GetXPDF(fl1)/x1;
  double fb=pdf2->GetXPDF(fl2)/x2;
  msg_Debugging()
      <<"pdf1 = ("<<fl1<<","<<x1<<","<<sqrt(muF12new)<<":"<<x1*fa<<") , "
      <<"pdf2 = ("<<fl2<<","<<x2<<","<<sqrt(muF22new)<<":"<<x2*fb<<")\n";
  double asf=pow((*MODEL::as)(muR2new)/(*MODEL::as)(muR2),oqcd);
  msg_Debugging()<<"asf="<<asf<<std::endl;
  if (mode==0) { // B,R,S
    msg_Debugging()<<"B,R,S event: new wgt="<<B*asf*fa*fb<<std::endl;
    return B*asf*fa*fb;
  }
  else if (mode==1) { // B,V,I
    if (wgts.size()<2) THROW(fatal_error,"Not enough weight information.");
    std::vector<double> w(9,0.);
    // B term (if only born order already the correct one)
    double asfborn(wgts[0]==0.?asf:pow(asf,(oqcd-1.)/oqcd));
    double bnew(B?B*asfborn*fa*fb:0.);
    msg_Debugging()<<" -> "<<bnew<<std::endl;
    // VI terms
    double lr=log(muR2fac), lf=log(muF2fac);
    w[0]=V+I+wgts[0]*lr+wgts[1]*0.5*sqr(lr);
    msg_Debugging()<<"VI = "<<w[0]*fa*fb;
    // KP terms
    double w0(w[0]*fa*fb);
    if (wgts.size()>=18) {
      for (int i(1);i<9;++i) {
        w[i]=wgts[i+1]+wgts[i+9]*lf;
      }
      double faq(0.0), faqx(0.0), fag(0.0), fagx(0.0);
      double fbq(0.0), fbqx(0.0), fbg(0.0), fbgx(0.0);
      Flavour quark(kf_quark), gluon(kf_gluon);
      if (w[1]!=0. || w[2]!=0. || w[3]!=0. || w[4]!=0.) {
        if (fl1.IsQuark()) {
          faq=fa;
          fag=pdf1->GetXPDF(gluon)/x1;
          pdf1->Calculate(x1/x1p,muF12new);
          faqx=pdf1->GetXPDF(fl1)/x1;
          fagx=pdf1->GetXPDF(gluon)/x1;
        }
        else if (fl1.IsGluon()) {
          fag=fa;
          for (size_t i=0;i<quark.Size();++i)
            faq+=pdf1->GetXPDF(quark[i])/x1;
          pdf1->Calculate(x1/x1p,muF12new);
          fagx=pdf1->GetXPDF(fl1)/x1;
          for (size_t i=0;i<quark.Size();++i)
            faqx+=pdf1->GetXPDF(quark[i])/x1;
        }
        else THROW(not_implemented,
                   std::string("Change of scales not implemented for ")
                   +ToString(fl1));
      }
      if (w[5]!=0. || w[6]!=0. || w[7]!=0. || w[8]!=0.) {
        if (fl2.IsQuark()) {
          fbq=fb;
          fbg=pdf2->GetXPDF(gluon)/x2;
          pdf2->Calculate(x2/x2p,muF22new);
          fbqx=pdf2->GetXPDF(fl2)/x2;
          fbgx=pdf2->GetXPDF(gluon)/x2;
        }
        else if (fl2.IsGluon()) {
          fbg=fb;
          for (size_t i=0;i<quark.Size();++i)
            fbq+=pdf2->GetXPDF(quark[i])/x2;
          pdf2->Calculate(x2/x2p,muF22new);
          fbgx=pdf2->GetXPDF(fl2)/x2;
          for (size_t i=0;i<quark.Size();++i)
            fbqx+=pdf2->GetXPDF(quark[i])/x2;
        }
        else THROW(not_implemented,
                   std::string("Change of scales not implemented for ")
                   +ToString(fl2));
      }
      w0+=(faq*w[1]+faqx*w[2]+fag*w[3]+fagx*w[4])*fb;
      w0+=(fbq*w[5]+fbqx*w[6]+fbg*w[7]+fbgx*w[8])*fa;
      msg_Debugging()<<" -> "<<w0;
    }
    else
    msg_Debugging()<<" -> "<<w0*asf<<std::endl;
    return bnew+w0*asf;
  }
  else {
    THROW(not_implemented,"Invalid option");
  }
  return 0.;
}

bool Scale_Variations::ComputeVariations(const ATOOLS::Weight_Info &winfo,
                                         PHASIC::Process_Base * proc)
{
  DEBUG_FUNC(proc);
  if (!m_on) return true;
  ResetValues();
  ExtractParameters(winfo,proc);
  for (NamedScaleVariationMap::iterator it=p_nsvmap->begin();
       it!=p_nsvmap->end();++it) {
    if (!Calculate(it->second)) return false;
  }
  return true;
}

namespace ATOOLS {
  template <> Blob_Data<NamedScaleVariationMap*>::~Blob_Data() {}
  template class Blob_Data<NamedScaleVariationMap*>;
}

