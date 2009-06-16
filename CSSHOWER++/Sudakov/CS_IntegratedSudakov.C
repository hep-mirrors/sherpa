#include "CSSHOWER++/Sudakov/CS_IntegratedSudakov.H"

#include "CSSHOWER++/Sudakov/CS_SingleSudakov_QCD.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Org/Exception.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

CS_IntegratedSudakov::CS_IntegratedSudakov(int _mode,MODEL::Running_AlphaS * _runas,
					   PDF::ISR_Handler * _isr) :
  m_mode(_mode), p_runas(_runas)
{
  if (_isr) p_isr=_isr;
  /*
  if (m_mode==0 && p_isr) 
    THROW(fatal_error,"IS shower disabled despite resolved initial state!");
  */
  for (int i=0;i<2; i++) p_pdf[i] = p_isr->PDF(i);

  for (int i=1;i<6;++i) {
    Flavour fl = Flavour(i);
    if (fl.IsOn() && fl.Strong()) {
      //the FF case
      Add(new DeltaFF_q_qg(fl));
      Add(new DeltaFF_q_qg(fl.Bar()));
      if (fl.Mass()<100.) Add(new DeltaFF_g_qq(fl));
      //initial state splittings
      if (m_mode==1) {
	//the FI case
	//Add(new q_qg_FI(fl,p_pdf));
	//Add(new q_qg_FI(fl.Bar(),p_pdf));
	//if (fl.Mass()<100.) Add(new g_qq_FI(fl,p_pdf));
	//the IF case
	//Add(new q_qg_IF(fl,p_pdf));
	//Add(new q_qg_IF(fl.Bar(),p_pdf));
	//Add(new q_gq_IF(fl,p_pdf));
	//Add(new q_gq_IF(fl.Bar(),p_pdf));
	//if (fl.Mass()<100.) Add(new g_qq_IF(fl,p_pdf));
	//if (fl.Mass()<100.) Add(new g_qq_IF(fl.Bar(),p_pdf));
	//the II case
	///Add(new q_qg_II(fl,p_pdf));
	//Add(new q_qg_II(fl.Bar(),p_pdf));
	//Add(new q_gq_II(fl,p_pdf));
	//Add(new q_gq_II(fl.Bar(),p_pdf));
	//if (fl.Mass()<100.) Add(new g_qq_II(fl,p_pdf));
	//if (fl.Mass()<100.) Add(new g_qq_II(fl.Bar(),p_pdf));
      }
    }
  }
  Add(new DeltaFF_g_gg());
  if (m_mode==1) {
    //Add(new g_gg_FI(p_pdf));
    //Add(new g_gg_IF(p_pdf));
    //Add(new g_gg_II(p_pdf));
  }

  //top splitting functions 
  Flavour fltop = Flavour(6);
  if (fltop.IsOn()) {
    //the FF case
      Add(new DeltaFF_q_qg(fltop));
      Add(new DeltaFF_q_qg(fltop.Bar()));
      //initial state splittings
      if (m_mode==1) {
	//the FI case
	//Add(new q_qg_FI(fltop,p_pdf));
	//Add(new q_qg_FI(fltop.Bar(),p_pdf));
      }
  }
  /*
  //susy splitting functions
  if (s_model->Name()==std::string("MSSM")) {
    for (int i=51;i<67;i++) {
      if (i==57) i=61;
      Flavour fl = Flavour(i);
      if (fl.IsOn()) {
	//the FF case
	Add(new sQ_sQg_FF(fl));
	Add(new sQ_sQg_FF(fl.Bar()));
	if (m_mode==1) {
	  //the FI case
	  Add(new sQ_sQg_FI(fl,p_pdf));
	  Add(new sQ_sQg_FI(fl.Bar(),p_pdf));
	}
      }
    }
    Add(new sG_sGg_FF());
    if (m_mode==1) {
      Add(new sG_sGg_FI(p_pdf));
    }
  }
  */
}

CS_IntegratedSudakov::~CS_IntegratedSudakov() {
  m_ssud=m_suds.begin();
  do {
    if (*m_ssud) { delete (*m_ssud); (*m_ssud=NULL); }
    m_ssud = m_suds.erase(m_ssud);
  } while (m_ssud!=m_suds.end());
  m_suds.clear();
}


void CS_IntegratedSudakov::Add(CS_SingleSudakov_Base * sud) {
  m_suds.push_back(sud);
}

double CS_IntegratedSudakov::Delta(const cstp::code dtype,const ATOOLS::Flavour & fl,  
				   const double Q2, const double mk2,const double kt2prod,
				   const double kt2dec) 
{
  double delta=1.;
  std::list<CS_SingleSudakov_Base *>::iterator sit=m_suds.begin();
  for (;sit!=m_suds.end();sit++) {
    if ((*sit)->GetFlavour()==fl && (*sit)->GetType()==dtype) {
      (*sit)->SetQ2(Q2);
      (*sit)->SetSpecMass(mk2);
      delta*=(*sit)->Delta(kt2dec,kt2prod);
    }
  }
  return delta;
}
