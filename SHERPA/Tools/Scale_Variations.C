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

#include <limits>

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

typedef void (*PDF_Init_Function)();
typedef void (*PDF_Exit_Function)();

Scale_Variation::Scale_Variation(const double &muR2fac, const double &muF2fac,
                                 PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
                                 MODEL::One_Running_AlphaS * as,
                                 bool deletepdfs, bool deleteas) :
  m_deletepdfs(deletepdfs), m_deleteas(deleteas),
  m_muR2fac(muR2fac), m_muF2fac(muF2fac), m_val(0.), m_RSvals(0,0.),
  p_pdf1(pdf1), p_pdf2(pdf2),
  m_pdf1id(pdf1->LHEFNumber()), m_pdf2id(pdf2->LHEFNumber()),
  m_pdf1set(pdf1->Set()), m_pdf2set(pdf2->Set()),
  m_pdf1setmember(pdf1->Member()), m_pdf2setmember(pdf2->Member()),
  p_as(as), m_name(GenerateName())
{
}

Scale_Variation::~Scale_Variation()
{
  if (m_deletepdfs) {
    if (p_pdf1) { delete p_pdf1; p_pdf1=NULL; }
    if (p_pdf2) { delete p_pdf2; p_pdf1=NULL; }
  }
  if (m_deleteas) {
    if (p_as) { delete p_as; p_as=NULL; }
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
  m_on(false), m_loadlhapdf(true), m_ckkw(false),
  m_quark(Flavour(kf_quark)), m_gluon(Flavour(kf_gluon)),
  p_nsvmap(new NamedScaleVariationMap())
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
  // check whether LHAPDF is already loaded, if not load and init interface
  if (!ATOOLS::s_loader->LibraryIsLoaded("LHAPDFSherpa")) {
    s_loader->AddPath(std::string(LHAPDF_PATH)+"/lib");
    ATOOLS::s_loader->LoadLibrary("LHAPDF");
    void *init(s_loader->GetLibraryFunction("LHAPDFSherpa","InitPDFLib"));
    if (init==NULL) THROW(fatal_error,"Cannot load PDF library LHAPDFSherpa");
    ((PDF_Init_Function)init)();
  }
  std::string path;
  if (reader.ReadFromFile(path,"LHAPDF_GRID_PATH")) LHAPDF::setPaths(path);
  const std::vector<std::string>& avsets(LHAPDF::availablePDFSets());
#endif
  for (size_t i(0);i<vars.size();++i) {
    std::string cur(vars[i]);
    size_t pos1(cur.find(",",0));
    double muR2fac(ToType<double>(std::string(cur,0,pos1)));
    size_t pos2(cur.find(",",pos1+1));
    double muF2fac(ToType<double>(std::string(cur,pos1+1,pos2)));
    if (pos2==std::string::npos) {
      // only use nominal PDF
      Scale_Variation * scvar(new Scale_Variation(muR2fac,muF2fac,
                                                  rpa->gen.PDF(0),
                                                  rpa->gen.PDF(1),
                                                  MODEL::as->
                                                  GetAs(PDF::isr::hard_process),
                                                  false,false));
      (*p_nsvmap)[scvar->Name()]=scvar;
    }
    else {
      // switch to other PDFs, only works with LHAPDF6
      // only then can initialise full set of error PDFs at once,
      // otherwise need to be given explicitly
      // for each PDF, extract PDF_AS_Info and init new One_Running_AlphaS
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
        MODEL::One_Running_AlphaS * as
            = new MODEL::One_Running_AlphaS(pdf1);
        if (!pdf1 || !pdf2)
          THROW(fatal_error,"PDF set "+pdfname+" not available.");
        if (!as)
          THROW(fatal_error,"AlphaS for "+pdfname+" could not be initialised.");
        Scale_Variation * scvar(new Scale_Variation(muR2fac,muF2fac,
                                                    pdf1,pdf2,as,true,true));
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
          MODEL::One_Running_AlphaS * as
              = new MODEL::One_Running_AlphaS(pdf1);
          Scale_Variation * scvar(new Scale_Variation(muR2fac,muF2fac,
                                                      pdf1,pdf2,as,true,true));
          (*p_nsvmap)[scvar->Name()]=scvar;
        }
#else
        THROW(not_implemented,"Full set reweightings only work with LHAPDF6."
                              +std::string(" Otherwise specify separately."));
#endif
      }
    }
  }
  msg_Info()<<*p_nsvmap;
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
       it!=p_nsvmap->end();++it) {
    it->second->SetValue(0.);
    it->second->DeleteRSValues();
  }
  m_params.dads.clear();
  m_params.csi=Cluster_Sequence_Info(1.,0.,0.);
}

void Scale_Variations::ExtractParameters(const ATOOLS::Weight_Info &winfo,
                                         PHASIC::Process_Base * proc)
{
  DEBUG_FUNC(proc->Name());
  const ME_Weight_Info * const mewgt(proc->GetMEwgtinfo());
  const NLO_subevtlist * const sevtlist(proc->GetSubevtList());
  if (!mewgt) THROW(fatal_error,"No ME_Weight_Info found.");
  msg_Debugging()<<*mewgt<<std::endl;
  m_params.type=mewgt->m_type;
  m_params.B=mewgt->m_B;
  m_params.VI=mewgt->m_VI;
  m_params.KP=mewgt->m_KP;
  m_params.oqcd=mewgt->m_oqcd;
  m_params.oew=mewgt->m_oew;
  m_params.fl1=mewgt->m_fl1;
  m_params.fl2=mewgt->m_fl2;
  m_params.x1=mewgt->m_x1;
  m_params.x2=mewgt->m_x2;
  m_params.x1p=mewgt->m_y1;
  m_params.x2p=mewgt->m_y2;
  m_params.renwgts=mewgt->m_wren;
  m_params.kpwgts=mewgt->m_wfac;
  m_params.muR2=mewgt->m_mur2;
  m_params.muF12=mewgt->m_muf2;
  m_params.muF22=mewgt->m_muf2;
  m_params.dads=mewgt->m_dadsinfos;
  m_params.rda=mewgt->m_rdainfos;
  m_params.csi=mewgt->m_clusseqinfo;
  if (sevtlist) {
    msg_Debugging()<<"NLO RS event contains "<<sevtlist->size()
                   <<" subevents"<<std::endl;
    for (NamedScaleVariationMap::iterator it=p_nsvmap->begin();
         it!=p_nsvmap->end();++it) {
      it->second->InitialisRSValues(sevtlist->size());
    }
    m_params.rswgts.resize(sevtlist->size(),0.);
    m_params.rsmuR2s.resize(sevtlist->size(),0.);
    m_params.rsmuF2s.resize(sevtlist->size(),0.);
    for (size_t i(0);i<sevtlist->size();++i) {
      m_params.rswgts[i]=(*sevtlist)[i]->m_mewgt;
      m_params.rsmuR2s[i]=(*sevtlist)[i]->m_mu2[stp::ren];
      m_params.rsmuF2s[i]=(*sevtlist)[i]->m_mu2[stp::fac];
    }
  }
  m_ckkw=m_params.csi.m_txfl.size();
}

bool Scale_Variations::Calculate(Scale_Variation * sv,
                                 PHASIC::Process_Base * proc)
{
  DEBUG_FUNC("event type: "<<m_params.type);
  if (proc->GetSubevtList()) {
    // NLO RS
    std::vector<double> renwgtsdummy,kpwgtsdummy;
    std::vector<DADS_Info> dadsdummy;
    std::vector<RDA_Info> rdadummy;
    Cluster_Sequence_Info csidummy;
    for (size_t i(0);i<proc->GetSubevtList()->size();++i) {
      sv->SetValue(i,Calculate(mewgttype::none,
                               m_params.rswgts[i],0.,
                               renwgtsdummy,kpwgtsdummy,
                               dadsdummy,rdadummy,csidummy,
                               m_params.x1,m_params.x2,
                               0.,0.,
                               m_params.fl1,m_params.fl2,
                               m_params.rsmuR2s[i],
                               m_params.rsmuF2s[i],m_params.rsmuF2s[i],
                               sv->MuR2Fac(),sv->MuF2Fac(),
                               m_params.oqcd,m_params.oew,
                               sv->PDF1(),sv->PDF2(),
                               sv->AlphaS()));
    }
  }
  else {
    // LO, LOPS, MEPS, MC@NLO S, NLO BVI
    sv->SetValue(Calculate(m_params.type,
                           m_params.B,m_params.VI,
                           m_params.renwgts,m_params.kpwgts,
                           m_params.dads,m_params.rda,m_params.csi,
                           m_params.x1,m_params.x2,
                           m_params.x1p,m_params.x2p,
                           m_params.fl1,m_params.fl2,
                           m_params.muR2,
                           m_params.muF12,m_params.muF22,
                           sv->MuR2Fac(),sv->MuF2Fac(),
                           m_params.oqcd,m_params.oew,
                           sv->PDF1(),sv->PDF2(),
                           sv->AlphaS()));
  }
  return true;
}

double Scale_Variations::Calculate(const ATOOLS::mewgttype::code& type,
                                   const double& B, const double& VI,
                                   const std::vector<double>& renwgts,
                                   const std::vector<double>& kpwgts,
                                   const std::vector<ATOOLS::DADS_Info>& dads,
                                   const std::vector<ATOOLS::RDA_Info>& rda,
                                   const Cluster_Sequence_Info& csi,
                                   const double& x1, const double& x2,
                                   const double& x1p, const double& x2p,
                                   const int& fl1, const int& fl2,
                                   const double& muR2,
                                   const double& muF12, const double& muF22,
                                   const double& muR2fac, const double& muF2fac,
                                   const size_t& oqcd, const size_t& oew,
                                   PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
                                   MODEL::One_Running_AlphaS * as)
{
  DEBUG_FUNC("factors (muR/muF)=("<<muR2fac<<","<<muF2fac<<"), "
             <<"pdf1="<<pdf1->LHEFNumber()<<", pdf2="<<pdf2->LHEFNumber());
  size_t precision(msg->Out().precision());
  if (msg_LevelIsDebugging()) msg->SetPrecision(15);
  // calculate new event weight
  double muR2new(muR2*muR2fac);
  double muF12new(muF12*muF2fac),muF22new(muF22*muF2fac);
  msg_Debugging()<<"B = "<<B<<", VI = "<<VI<<std::endl;
  msg_Debugging()<<"muR: "<<muR2<<" -> "<<muR2new<<" , "
                 <<"muF1: "<<muF12<<" -> "<<muF12new<<" , "
                 <<"muF2: "<<muF22<<" -> "<<muF22new<<std::endl;
  msg_Debugging()<<"oqcd: "<<oqcd<<", oew: "<<oew<<std::endl;
  msg_Debugging()<<"cluster sequence: "<<csi.m_txfl.size()<<" steps"<<std::endl;
  // build type minus METS
  mewgttype::code nometstype=((type&mewgttype::METS)?type^mewgttype::METS:type);
  pdf1->Calculate(x1,muF12new);
  pdf2->Calculate(x2,muF22new);
  double fa(pdf1->GetXPDF(fl1)/x1);
  double fb(pdf2->GetXPDF(fl2)/x2);
  msg_Debugging()
      <<"pdf1 = ("<<fl1<<","<<x1<<","<<sqrt(muF12new)<<") = "<<fa<<"\n"
      <<"pdf2 = ("<<fl2<<","<<x2<<","<<sqrt(muF22new)<<") = "<<fb<<std::endl;
  msg_Debugging()<<"fa*fb="<<fa*fb<<std::endl;
  double pdffac(PDFRatioFactor(fa,fb,csi,muF2fac,pdf1,pdf2));
  msg_Debugging()<<"fa*fb*pdffac="<<fa*fb*pdffac<<std::endl;
  // reset MODEL::as to hard process
  MODEL::as->SetActiveAs(PDF::isr::hard_process);
  double asnew((*as)(muR2new)),asold((*MODEL::as)(muR2));
  if (nometstype==mewgttype::none) { // B,R,S
    double asf=pow(asnew/asold,oqcd);
    msg_Debugging()<<"asf = "<<asf<<std::endl;
    msg_Debugging()<<"B,R,S event: new wgt="<<B*asf*fa*fb*pdffac<<std::endl;
    if (msg_LevelIsDebugging()) msg->SetPrecision(precision);
    return B*asf*fa*fb*pdffac;
  }
  else if (nometstype==mewgttype::H) { // H
    double RDAnew(0.);
    for (size_t i(0);i<rda.size();++i) {
      double RDAinew(0.);
      if (rda[i].m_wgt!=0.) {
        double asrdainew((*as)(rda[i].m_mur2*muR2fac)),
               asrdaiold((*MODEL::as)(rda[i].m_mur2));
        double asfrdai=pow(asrdainew/asrdaiold,oqcd);
        msg_Debugging()<<"asf = "<<asfrdai<<std::endl;
        pdf1->Calculate(x1,rda[i].m_muf12*muF2fac);
        pdf2->Calculate(x2,rda[i].m_muf22*muF2fac);
        double farda=pdf1->GetXPDF(fl1)/x1;
        double fbrda=pdf2->GetXPDF(fl2)/x2;
        RDAinew=farda*fbrda*asfrdai*rda[i].m_wgt*pdffac;
      }
      msg_Debugging()<<"  new RDA_"<<i<<" = "<<RDAinew<<std::endl;
      RDAnew+=RDAinew;
    }
    msg_Debugging()<<"H event: new wgt="<<RDAnew<<std::endl;
    if (msg_LevelIsDebugging()) msg->SetPrecision(precision);
    return RDAnew;
  }
  else { // B,VI,KP,DADS
    // B term (if only born order already the correct one)
    double asf=pow(asnew/asold,oqcd);
    double asfborn((renwgts[0]==0.&&renwgts[1]==0.)?asf:
                                                    pow(asnew/asold,oqcd-1));
    msg_Debugging()<<"asf(B) = "<<asfborn<<std::endl;
    msg_Debugging()<<"asf(VI,KP) = "<<asf<<std::endl;
    double Bnew(B*asfborn*fa*fb*pdffac);
    msg_Debugging()<<"new B = "<<Bnew<<std::endl;
    // VI terms
    double lr=log(muR2fac);
    double VInew(VI+renwgts[0]*lr+renwgts[1]*0.5*ATOOLS::sqr(lr));
    VInew*=asf*fa*fb*pdffac;
    msg_Debugging()<<"new VI = "<<VInew<<std::endl;
    // KP terms
    double lf=log(muF2fac);
    std::vector<double> w(8,0.);
    for (int i(0);i<8;++i) w[i]=kpwgts[i]+kpwgts[i+8]*lf;
    double faq(0.0), faqx(0.0), fag(0.0), fagx(0.0);
    double fbq(0.0), fbqx(0.0), fbg(0.0), fbgx(0.0);
    if (w[0]!=0. || w[1]!=0. || w[2]!=0. || w[3]!=0.) {
      Flavour flav1(abs(fl1),fl1<0);
      if (flav1.IsQuark()) {
        faq=fa;
        fag=pdf1->GetXPDF(m_gluon)/x1;
        pdf1->Calculate(x1/x1p,muF12new);
        faqx=pdf1->GetXPDF(fl1)/x1;
        fagx=pdf1->GetXPDF(m_gluon)/x1;
      }
      else if (flav1.IsGluon()) {
        fag=fa;
        for (size_t i=0;i<m_quark.Size();++i)
          faq+=pdf1->GetXPDF(m_quark[i])/x1;
        pdf1->Calculate(x1/x1p,muF12new);
        fagx=pdf1->GetXPDF(fl1)/x1;
        for (size_t i=0;i<m_quark.Size();++i)
          faqx+=pdf1->GetXPDF(m_quark[i])/x1;
      }
      else THROW(not_implemented,
                 std::string("Change of scales not implemented for ")
                 +ToString(fl1));
    }
    if (w[4]!=0. || w[5]!=0. || w[6]!=0. || w[7]!=0.) {
      Flavour flav2(abs(fl2),fl2<0);
      if (flav2.IsQuark()) {
        fbq=fb;
        fbg=pdf2->GetXPDF(m_gluon)/x2;
        pdf2->Calculate(x2/x2p,muF22new);
        fbqx=pdf2->GetXPDF(fl2)/x2;
        fbgx=pdf2->GetXPDF(m_gluon)/x2;
      }
      else if (flav2.IsGluon()) {
        fbg=fb;
        for (size_t i=0;i<m_quark.Size();++i)
          fbq+=pdf2->GetXPDF(m_quark[i])/x2;
        pdf2->Calculate(x2/x2p,muF22new);
        fbgx=pdf2->GetXPDF(fl2)/x2;
        for (size_t i=0;i<m_quark.Size();++i)
          fbqx+=pdf2->GetXPDF(m_quark[i])/x2;
      }
      else THROW(not_implemented,
                 std::string("Change of scales not implemented for ")
                 +ToString(fl2));
    }
    double KPnew(0.);
    KPnew+=(faq*w[0]+faqx*w[1]+fag*w[2]+fagx*w[3])*fb;
    KPnew+=(fbq*w[4]+fbqx*w[5]+fbg*w[6]+fbgx*w[7])*fa;
    KPnew*=asf*pdffac;
    msg_Debugging()<<"new KP = "<<KPnew<<std::endl;
    // DADS terms
    double DADSnew(0.);
    for (size_t i(0);i<dads.size();++i) {
      double DADSinew(0.);
      if (dads[i].m_wgt!=0.) {
        // muR, muF must be same as in B
        pdf1->Calculate(dads[i].m_x1,muF12new);
        pdf2->Calculate(dads[i].m_x2,muF22new);
        double fadadsi=pdf1->GetXPDF(dads[i].m_fl1)/dads[i].m_x1;
        double fbdadsi=pdf2->GetXPDF(dads[i].m_fl2)/dads[i].m_x2;
        DADSinew=fadadsi*fbdadsi*asf*dads[i].m_wgt*pdffac;
      }
      msg_Debugging()<<"  new DADS_"<<i<<" = "<<DADSinew<<std::endl;
      DADSnew+=DADSinew;
    }
    msg_Debugging()<<"new DADS = "<<DADSnew<<std::endl;
    msg_Debugging()<<"B,V,I event: new wgt="<<Bnew+VInew+KPnew+DADSnew
                   <<std::endl;
    if (msg_LevelIsDebugging()) msg->SetPrecision(precision);
    return Bnew+VInew+KPnew+DADSnew;
  }
  if (msg_LevelIsDebugging()) msg->SetPrecision(precision);
  return 0.;
}

double Scale_Variations::PDFRatioFactor
(const double& fa, const double& fb, const ATOOLS::Cluster_Sequence_Info& csi,
 const double& muF2fac, PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2)
{
  DEBUG_FUNC("ckkw="<<m_ckkw<<", #steps="<<csi.m_txfl.size());
  if (!m_ckkw) return 1.;
  // if no cluster history, something is wrong
  if (csi.m_txfl.size()<2) THROW(fatal_error,"Insufficient cluster history.");
  // want to calculate (i=0 -> core, i=N -> ext)
  // wn-ext * [\prod_{i=0}^{N-1} wn_i/wd_i]
  // = [wn-ext * \prod_{i=1}^{N-1} wn_i/wd_i * 1/wd_0] * wn-core
  // = [\prod_{i=1}^N wn_i/wd_{i-1}] * wn-core
  // and vary only wn-core, but have varied fa*fb = wn-ext so far
  // -------------------------------------------------------------------
  // this computes wn-ext * \prod_{i=1}^{N-1} wn_i/wd_i * 1/wd_0
  // multiply wn-ext by 1/fa*fb*flux to cancel factor coming from outside
  double fac(csi.m_txfl[0].m_pdfrationumerator/(fa*fb));
  double t(std::numeric_limits<double>::max());
  msg_Debugging()<<"pdffac="<<fac<<std::endl;
  bool lastoneordered(false);
  for (size_t i(1);i<csi.m_txfl.size();++i) {
    msg_Debugging()<<i<<": "<<csi.m_txfl[i];
    if (i<2 || t<csi.m_txfl[i].m_t) {
      msg_Debugging()<<": "<<csi.m_txfl[i].m_pdfrationumerator
                     <<" / "<<csi.m_txfl[i].m_pdfratiodenominator;
      fac*=(csi.m_txfl[i].m_pdfrationumerator/
            csi.m_txfl[i].m_pdfratiodenominator);
      t=csi.m_txfl[i].m_t;
      if (i==csi.m_txfl.size()-1) lastoneordered=true;
    }
    else msg_Debugging()<<": Skip. Unordered history "
                        <<sqrt(t)<<" > "<<sqrt(csi.m_txfl[i].m_t)<<std::endl;
    else t=std::numeric_limits<double>::max();
    msg_Debugging()<<std::endl;
  }
  msg_Debugging()<<"pdffac="<<fac<<std::endl;
  // multiply with new wn-core if last step was ordered
  if (lastoneordered) {
    msg_Debugging()<<"last step ordered, apply new core PDF:"<<std::endl;
    double x1(csi.m_txfl.back().m_xa),x2(csi.m_txfl.back().m_xb);
    pdf1->Calculate(x1,t*muF2fac);
    pdf2->Calculate(x2,t*muF2fac);
    double fcorea(pdf1->GetXPDF(csi.m_txfl.back().m_fla)/x1);
    double fcoreb(pdf2->GetXPDF(csi.m_txfl.back().m_flb)/x2);
    msg_Debugging()
        <<"  pdf1 = ("<<csi.m_txfl.back().m_fla<<","<<x1<<","
        <<sqrt(t*muF2fac)<<") = "<<fcorea<<"\n"
        <<"  pdf2 = ("<<csi.m_txfl.back().m_fla<<","<<x2<<","
        <<sqrt(t*muF2fac)<<") = "<<fcoreb<<std::endl;
    msg_Debugging()<<"fa*fb="<<fcorea*fcoreb<<std::endl;
    fac*=fcorea*fcoreb;
  }
  else msg_Debugging()<<"last step unordered, no core PDF variation"<<std::endl;
  msg_Debugging()<<"pdffac="<<fac<<std::endl;
  return fac;
}


bool Scale_Variations::ComputeVariations(const ATOOLS::Weight_Info &winfo,
                                         PHASIC::Process_Base * proc)
{
  DEBUG_FUNC(proc->Name());
  if (!m_on) return true;
  ResetValues();
  ExtractParameters(winfo,proc);
  for (NamedScaleVariationMap::iterator it=p_nsvmap->begin();
       it!=p_nsvmap->end();++it) {
    if (!Calculate(it->second,proc)) return false;
  }
  return true;
}

namespace ATOOLS {
  template <> Blob_Data<NamedScaleVariationMap*>::~Blob_Data() {}
  template class Blob_Data<NamedScaleVariationMap*>;
}

