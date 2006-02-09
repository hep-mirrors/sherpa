//bof
//Version: 4 ADICIC++-0.0/2006/02/09

//Implementation of Sudakov_Calculator.H.



#include "Running_AlphaS.H"////////////////////////////////////////////////////
#include "Data_Read.H"
#include "PDF_Handler.H"
#include "Sudakov_Calculator.H"





using namespace std;
using namespace ATOOLS;
using namespace PDF;
using namespace ADICIC;





//#include "Sudakov_Calculator.tpt.cc"





//=============================================================================



//Mimic Ariadne.
const bool Sudakov_Calculator::sf_ariadne=false;//true;//false;
const bool& Sudakov_Calculator::Ariadne=Sudakov_Calculator::sf_ariadne;

//Temporary PDF switch.
const bool Sudakov_Calculator::sf_pdf=true;//false;

//So far there is no static Sudakov_Calculator.
int Sudakov_Calculator::s_count=0;
const int& Sudakov_Calculator::InStore=Sudakov_Calculator::s_count;

//Approximation variables.
double Sudakov_Calculator::s_asapprox=dpa.sud.AlphaSFix();

//
Sudakov_Calculator::Toolbox Sudakov_Calculator::s_box
=Sudakov_Calculator::Toolbox();

//
Sudakov_Calculator::Double_Double_Func
Sudakov_Calculator::GetAlphaSCorr=&Sudakov_Calculator::FixAlphaSCorr;
Sudakov_Calculator::Int_Double_Func
Sudakov_Calculator::GetNf=&Sudakov_Calculator::FixNf;

//
Sudakov_Calculator::PDF_Corr_Func
Sudakov_Calculator::GetPDFCorr=&Sudakov_Calculator::NoPDFCorr;



//=============================================================================



int Sudakov_Base::s_count=0;
const int& Sudakov_Base::InStore=Sudakov_Base::s_count;

//
const double Sudakov_Base::s_average=2.0;



//=============================================================================



Sudakov_Calculator::~Sudakov_Calculator() {    //Virtual.
  assert(p_dip==NULL);
  assert(p_sur==NULL);
  --s_count;
#ifdef TEMP_OUTPUT
  cout<<"~Sudakov_Calculator"<<endl;///////////////////////////////////////////
#endif
}





void Sudakov_Calculator::ShowEnvironment() {    //Static.
  cout<<endl;
  cout<<"====================================================="<<endl;
  cout<<"         Valid Sudakov_Calculator parameters"<<endl;
  cout<<"-----------------------------------------------------"<<endl;
  cout<<"Running AlphaS is wished for........."
      <<dpa.sud.RunAlphaS()<<".\n";
  cout<<"Running AlphaS is initialized........"
      <<bool(s_box.m_ras[0])<<".\n";
  cout<<"Fixed AlphaS is set to..............."
      <<dpa.sud.AlphaSFix()<<".\n";
  cout<<"Fixed Nf is set to..................."
      <<dpa.sud.NfFix()<<".\n";
  cout<<"AlphaS approximation is set to......."
      <<s_asapprox<<".\n";
  cout<<"PDF treatment is wished for.........."
      <<bool(sf_pdf)<<".\n";
  cout<<"PDFs are initialized................."
      <<bool(s_box.m_pdf[0] && s_box.m_pdf[1])<<".\n";
  cout<<"Radiation type is set to............."
      <<dpa.sud.RadiationType()<<".\n";
  cout<<"FF dipole shower cut-off scale is set to...."
      <<dpa.sud.MinK2t()<<" GeV^2.\n";
  cout<<"FF dipole shower maximum scale is set to...."
      <<dpa.sud.MaxK2t()<<" GeV^2.\n";
  cout<<"II dipole shower cut-off scale is set to...."
      <<dpa.sud.MinIIK2t()<<" GeV^2.\n";
  cout<<"II dipole shower maximum scale is set to...."
      <<dpa.sud.MaxIIK2t()<<" GeV^2.\n";
  cout<<"II eff. enhancement exponent is............."
      <<dpa.sud.IIEffExp()<<".\n";
  cout<<"====================================================="<<endl;
}





const Trio Sudakov_Calculator::AdjustEnvironment(const string& path,
						 MODEL::Model_Base* pmod) {

  //Static method.

  //static MODEL::Running_AlphaS locras(0.1188,8315.25,1);
  static MODEL::Running_AlphaS locras(0.118,8315.0,1);

  Trio ret;

#ifdef DIPOLE_PARAMETER_OUTPUT
  cout<<"{ "<<__PRETTY_FUNCTION__<<" ...\n";
#endif

  if(dpa.sud.RunAlphaS()==false) {
    s_box.m_ras[0]=NULL;
    s_asapprox=dpa.sud.AlphaSFix();
    GetAlphaSCorr=&FixAlphaSCorr;
    GetNf=&FixNf;
    ret=Nil;
  } else {
    if(s_box.m_ras[0]) {
      if(pmod) {
#ifdef DIPOLE_PARAMETER_OUTPUT
	cout<<"  Re-initialize alphaS treatment as global.\n";
#endif
	s_box.m_ras[0]=pmod->GetScalarFunction("alpha_S");
	ret=Positive;
      } else {
	string s;
	if(s_box.m_ras[0]==&locras) { ret=Negative; s="local";}
	else { ret=Positive; s="global";}
#ifdef DIPOLE_PARAMETER_OUTPUT
	cout<<"  Keep the alphaS treatment "<<s<<".\n";
#endif
      }
    } else {
      if(pmod) {
	//The Running_AlphaS object physically resides in the initialized
	//Model ==> global treatment.
	//So the following is simply an assignment (no new-operator is needed).
#ifdef DIPOLE_PARAMETER_OUTPUT
	cout<<"  Initialize alphaS treatment as global.\n";
#endif
	s_box.m_ras[0]=pmod->GetScalarFunction("alpha_S");
	ret=Positive;
      } else {
	//Local treatment, see on top of this method.
#ifdef DIPOLE_PARAMETER_OUTPUT
	cout<<"  Initialize alphaS treatment as local.\n";
#endif
	s_box.m_ras[0]=&locras;
	ret=Negative;
      }
    }
    assert(s_box.m_ras[0]);
    double scmin=Min(dpa.sud.MinK2t(),dpa.sud.MinIIK2t());
    double cutq2=static_cast<MODEL::Running_AlphaS*>(s_box.m_ras[0])->CutQ2();
    if(scmin<cutq2) s_asapprox=(*s_box.m_ras[0])(cutq2)+0.0001;
    else            s_asapprox=(*s_box.m_ras[0])(scmin)+0.0001;
    //s_asapprox=(*s_box.m_ras[0])(cutq2)+0.0001;//////////////////////////////
    //cout<<cutq2<<" : "<<scmin<<" :: "<<s_asapprox<<endl;
    //for(int i=0; i<100; ++i)
    //  cout<<(1.0-i/100.0)<<"\t"<<(*s_box.m_ras[0])(1.0-i/100.0)<<endl;
    assert(s_asapprox>(*s_box.m_ras[0])(scmin) &&
	   s_asapprox>(*s_box.m_ras[0])(dpa.sud.MinK2t()) &&
	   s_asapprox>(*s_box.m_ras[0])(dpa.sud.MinIIK2t()));
    GetAlphaSCorr=&RunAlphaSCorr;
    GetNf=&RunNf;
  }

  if(sf_pdf==false) {
    if(s_box.m_pdf[0]) { delete s_box.m_pdf[0]; s_box.m_pdf[0]=NULL;} 
    if(s_box.m_pdf[1]) { delete s_box.m_pdf[1]; s_box.m_pdf[1]=NULL;}
    GetPDFCorr=&NoPDFCorr;
#ifdef DIPOLE_PARAMETER_OUTPUT
    cout<<"  PDFs are not initialized.\n";
#endif
  } else {
    bool h=true;
    if(s_box.m_pdf[0] && s_box.m_pdf[1]) {
#ifdef DIPOLE_PARAMETER_OUTPUT
      cout<<"  PDFs are already initialized.\n";
#endif
    } else {
      assert(!s_box.m_pdf[0] && !s_box.m_pdf[1]);
      Flavour fl0, fl1;
      string thepath;
      if(path=="default") thepath="../TestIt/data/ISR.dat";
      else thepath=path+"ISR.dat";
      Data_Read dataread(thepath,true);
      if(dataread.FileExists()) {
	if(dataread.GetValue<Switch::code>("ISR_1")==Switch::On &&
	   dataread.GetValue<Switch::code>("ISR_2")==Switch::On) {
	  PDF_Handler pdfhandler;
	  s_box.m_pdf[0]=pdfhandler.GetPDFLib(&dataread,fl0,0);
	  s_box.m_pdf[1]=pdfhandler.GetPDFLib(&dataread,fl1,1);
	  s_box.m_pdf[0]->SetRenormalizationScaleFactor(1.0);    //To be sure.
	  s_box.m_pdf[1]->SetRenormalizationScaleFactor(1.0);    //To be sure.
	  assert(fl0==Flavour(kf::p_plus));
	  assert(fl1==Flavour(kf::p_plus,1));
#ifdef DIPOLE_PARAMETER_OUTPUT
	  cout<<"  PDFs have been initialized.\n";
#endif
	} else {
	  h=false;
#ifdef DIPOLE_PARAMETER_OUTPUT
	  cout<<"  PDFs could not be initialized. ISR switched off!\n";
#endif
	}
      } else {
	h=false;
#ifdef DIPOLE_PARAMETER_OUTPUT
	cout<<"  PDFs could not be initialized. File not found!\n";
#endif
      }
    }
    if(h) {
      assert(s_box.m_pdf[0] && s_box.m_pdf[1]);
      GetPDFCorr=&IsPDFCorr;
    } else {
      assert(!s_box.m_pdf[0] && !s_box.m_pdf[1]);
      GetPDFCorr=&NoPDFCorr;
    }
  }

#ifdef DIPOLE_PARAMETER_OUTPUT
  cout<<"}\n";    //assert(0);
#endif

  return ret;

}





void Sudakov_Calculator::Which() const {    //Virtual.
  cout<<"Incomplete Sudakov_Group object!"<<endl;
}



//-----------------------------------------------------------------------------



const double Sudakov_Calculator::IsPDFCorr(bool z, const Multiflavour& mufl,
					   const Multidouble& mudo) {

  //Static method.
  //cout<<sf::plusini+z<<","<<sf::plusfin+z<<endl;/////////////////////////////
  //cout<<sr::xpini+z<<","<<sr::xpfin+z<<endl;/////////////////////////////////

  static bool speedup=true;    //Allow for speedup or not.
  static PDF_Base* pdf[2]={NULL,NULL};
  static Flavour fla[2]={Flavour(),Flavour()};
  static double xin[2]={0.0,0.0};
  static double q2i[2]={0.0,0.0};
  static double wden[2];

  bool same=false;

  s_box.m_pdf[z]->Calculate(mudo[sr::xpfin+z],mudo[sr::shat]);    //Or mperp^2?
  double wnum=s_box.m_pdf[z]->GetXPDF(mufl[sf::plusfin+z])/mudo[sr::xpfin+z];

  if( !speedup || pdf[z]!=s_box.m_pdf[z] ||
      xin[z]!=mudo[sr::xpini+z] || q2i[z]!=sqr(mudo[sr::mdip]) ) {
    s_box.m_pdf[z]->Calculate(mudo[sr::xpini+z],sqr(mudo[sr::mdip]));
    wden[z]=s_box.m_pdf[z]->GetXPDF(mufl[sf::plusini+z])/mudo[sr::xpini+z];
    pdf[z]=s_box.m_pdf[z];
    fla[z]=mufl[sf::plusini+z];
    xin[z]=mudo[sr::xpini+z];
    q2i[z]=sqr(mudo[sr::mdip]);
  } else {
    if(fla[z]!=mufl[sf::plusini+z]) {
      wden[z]=s_box.m_pdf[z]->GetXPDF(mufl[sf::plusini+z])/mudo[sr::xpini+z];
      fla[z]=mufl[sf::plusini+z];
    } else {
      same=true;
      //Then take the kept value of wden[z].
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //cout<<"       "<<same<<"   "
  //    <<mufl[sf::plusfin+z]<<", "<<mufl[sf::plusini+z]<<"    "
  //    <<mudo[sr::xpfin+z]<<", "<<mudo[sr::xpini+z]<<"    "
  //    <<mudo[sr::shat]<<", "<<sqr(mudo[sr::mdip])<<"   |   "
  //    <<wnum<<", "<<wden[z]<<"  =>  "<<(wnum/wden[z])<<"\n";
  /////////////////////////////////////////////////////////////////////////////

  wnum/=wden[z];
  assert(wnum>=0.0);
  //assert(wnum<=value);
  return wnum;

}



//=============================================================================



Sudakov_Calculator::Toolbox::Toolbox() : m_ras(1), m_pdf(2) {
  assert(m_ras[0]==NULL && m_pdf[0]==NULL && m_pdf[1]==NULL);
}



Sudakov_Calculator::Toolbox::~Toolbox() {
#ifdef TEMP_OUTPUT
  cout<<"~Sudakov_Calculator::Toolbox(Funcs:"<<s_box.m_ras[0]
      <<";Pdfs:"<<s_box.m_pdf[0]<<","<<s_box.m_pdf[1]<<")\n";
#endif
  //Responsibility of deleting running alphaS is elsewhere!
  //for(size_t i=0; i<m_ras.size(); ++i) if(m_ras[i]) delete m_ras[i];
  //Responsibility of deleting the PDFs is here!
  for(size_t i=0; i<m_pdf.size(); ++i) if(m_pdf[i]) delete m_pdf[i];
}



//=============================================================================



Sudakov_Base::~Sudakov_Base() {    //Virtual.
  --s_count;
#ifdef TEMP_OUTPUT
  std::cout<<"~Sudakov_Base"<<std::endl;//////////////////////////
#endif
}





void Sudakov_Base::Which() const {    //Virtual.
  cout<<"Incomplete Sudakov object!"<<endl;
}



//=============================================================================





//eof
