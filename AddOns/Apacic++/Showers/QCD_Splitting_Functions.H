#ifndef APACIC_Showers_QCD_Splitting_Functions_H
#define APACIC_Showers_QCD_Splitting_Functions_H
//----------------------------------------------------------------------
//               QCD_Splitting_Functions.H
//----------------------------------------------------------------------
#include "AddOns/Apacic++/Showers/Splitting_Function.H"

namespace APACIC {
  class Sudakov_Tools;

  // quark to quark + gluon  splitting function
  class q_qg : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha, m_kfactor, m_fmed;
  public:
    q_qg(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour quarkflavour,double fmed=0.);
    q_qg(ATOOLS::Mass_Selector *&ms,
	 ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools,double fmed=0.);
    double operator()(double z);
    double GetZ();
    double GetPhi(double z);
    const ATOOLS::Simple_Polarisation_Info GetPolB(double z, double phi);
    const ATOOLS::Simple_Polarisation_Info GetPolC(double z, double phi, double phi_b);
    
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool massterm);
    double CrudeInt(double _zmin, double _zmax);
    double Integral(double zmin, double zmax);
  };

  // gluon to gluon + gluon splitting function (needed twice for initial state shower)
  class g_gg : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha, m_kfactor, m_fmed;
  public:
    g_gg(ATOOLS::Mass_Selector *&ms,double fmed=0.);
    g_gg(ATOOLS::Mass_Selector *&ms,Sudakov_Tools * _tools,double fmed=0.);
    double operator()(double z);
    double GetZ();
    double GetPhi(double z);
    const ATOOLS::Simple_Polarisation_Info GetPolB(double z, double phi);
    const ATOOLS::Simple_Polarisation_Info GetPolC(double z, double phi, double phi_b);
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool masses);
    double CrudeInt(double _zmin, double _zmax);
    double Integral(double zmin, double zmax);
  };

  // gluon to quark + anti-quark splitting function
  class g_qq : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha, m_fmed;
  public:
    g_qq(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour quarkflavour,double fmed=0.);
    g_qq(ATOOLS::Mass_Selector *&ms,
	 ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools,double fmed=0.);
    double operator()(double z);
    double GetZ();
    double GetPhi(double z);
    const ATOOLS::Simple_Polarisation_Info GetPolB(double z, double phi);
    const ATOOLS::Simple_Polarisation_Info GetPolC(double z, double phi, double phi_b);
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool masses);
    double CrudeInt(double _zmin, double _zmax);
    double Integral(double zmin, double zmax);
  };


  // quark to qluon + quark splitting function (only used in Initial State Shower)
  class q_gq : public Splitting_Function {
    Sudakov_Tools * p_tools;
    double          m_alpha, m_kfactor, m_fmed;
  public:
    q_gq(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour quarkflavour,double fmed=0.);
    q_gq(ATOOLS::Mass_Selector *&ms,
	 ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools,double fmed=0.);
    double operator()(double z);
    double GetZ();
    double GetPhi(double z);
    const ATOOLS::Simple_Polarisation_Info GetPolB(double z, double phi);
    const ATOOLS::Simple_Polarisation_Info GetPolC(double z, double phi, double phi_b);
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool massterm);
    double CrudeInt(double _zmin, double _zmax);
    double Integral(double zmin, double zmax);
  };

  /*!
    \file
    \brief contains the classes APACIC::q_qg, APACIC::q_gq, 
           APACIC::g_gg and APACIC::g_qq
  */

  //##############################################################################

  /*!
    \class q_qg
    \brief represents a quark to quark + gluon  splitting function

    This class represents the splitting function for a \f$q \to q \; g \f$
    splitting:
    \f[
       P_{q \to q g} (z) = C_F \frac{1 + z^2}{1-z}
    \f]

  */
  /*!

    \var Sudakov_Tools * p_tools

  */
  /*!

    \var double m_alpha

  */
  /*!

    \var double m_kfactor

  */
  /*!

    \var double m_fmed

  */
  /*!

    \fn q_qg::q_qg(ATOOLS::Flavour quarkflavour)
    \brief Constructor, taking the quark flavour as argument (obsolete)

  */
  /*!

    \fn q_qg::q_qg(ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools)
    \brief Constructor, takng the quark flavour and a class to
           deliver running alpha_s 

  */
  /*!

    double operator()(double z)

  */
  /*!

    double GetZ()

  */
  /*!

    double GetPhi(double z)

  */
  /*!

    const ATOOLS::Simple_Polarisation_Info GetPolB(double z, double phi)

  */
  /*!

    const ATOOLS::Simple_Polarisation_Info GetPolC(double z, double phi, double phi_b)

  */
  /*!
    
    double GetCoupling()

  */
  /*!

    double GetCoupling(double t)

  */
  /*!

    double GetWeight(double z,double pt2,bool massterm)

  */
  /*!

    double CrudeInt(double _zmin, double _zmax)

  */
  /*!

    double Integral(double zmin, double zmax)

  */



  //##############################################################################

  /*!
    \class g_gg
    \brief represents a gluon to gluon + gluon splitting function 
           (needed twice for initial state shower)

    This class represents the splitting function for a \f$g \to g \ g \f$
    splitting:
    \f[
       P_{g \to g g} (z) = C_A \frac{[1 - z (1 - z)]^2}{z(1-z)}
    \f]
  */
  /*!
    \fn g_gg::g_gg()
    \brief Standard Constructor
  */
  /*!
    \fn g_gg::g_gg(Sudakov_Tools * _tools)
    \brief Constructor, taking a class to deliver running alpha_s
  */

  //##############################################################################

  /*!
    \class g_qq
    \brief represents gluon to quark + anti-quark splitting function
 
    This class represents the splitting function for a \f$g \to q \ \bar q \f$
    splitting:
    \f[
       P_{g \to q \bar q} (z) = T_R [z^2 + (1-z)^2]
    \f]

 */
  /*!
    \fn g_qq::g_qq(ATOOLS::Flavour quarkflavour)
    \brief Constructor, taking the quark flavour as argument (obsolete)
  */
  /*!
    \fn g_qq::g_qq(ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools)
    \brief Constructor, takng the quark flavour and a class to
           deliver running alpha_qed 
  */

  //##############################################################################

  /*!
    \class q_gq
    \brief quark to qluon + quark splitting function 
           (only used in Initial State Shower)


    This class represents the splitting function for a \f$q \to g \ q \f$
    splitting:
    \f[
       P_{q \to g q} (z) = CF \frac{1 + (1-z)^2}{z}
    \f]

  */
  /*!
    \fn q_gq::q_gq(ATOOLS::Flavour quarkflavour)
    \brief Constructor, taking the quark flavour as argument (obsolete)
  */
  /*!
    \fn q_gq::q_gq(ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools)
    \brief Constructor, takng the quark flavour and a class to
           deliver running alpha_qed 
  */


} // namespace APACIC
#endif








