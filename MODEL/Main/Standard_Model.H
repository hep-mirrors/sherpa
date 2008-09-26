#ifndef Standard_Model_h
#define Standard_Model_h

#include "Model_Base.H"

namespace MODEL {
  class Standard_Model : public Model_Base {
    int  m_trivialextension;
  private :
    int  m_ewscheme,m_ckmorder;
    void ParticleInit();
    void FillSpectrum();
    void FixEWParameters();
    void FixCKM();
  public :
    Standard_Model(std::string,std::string,bool,const int =0);
    ~Standard_Model()                   {}
    bool FillDecay(ATOOLS::Decay_Table *) { return false; }
  };
  /*!
    \class Standard_Model
    \brief contains the Standard_Model (no surprise)
    
    There are different options implemented how to define the elctrowek sector of the SM, all of them
    steered through m_ewscheme that can be read in. 
  */
  /*!
    \var int Standard_Model::m_trivialextension
    \brief To steer and select some trivial extensions of the SM.

    In Sherpa, a number of trivial extensions of the Standard Model have been implemented,
    labelled by m_trivialextension:
    -0 : No extension
    -1 : Non-standard top-bottom-W couplings with a left- and a right-handed part.
         This includes two parameters:
	 - \f$\kappa_{tbW}\f$, the relative strength of the \f$tbW\f$-coupling w.r.t.
	   the SM value,
	 - and \f$\theta_{tbW}\f$ an angle to steer the ratio of left-handed (\f$\cos\theta_{tbW}\f$)
           and right-handed (\f$\sin\theta_{tbW}\f$) coupling. 
    -2 : A fourth generation with a slightly enlarged coupling sector (due to new
         quark and lepton mixing terms like in CKM.
  */
  /*!
    \var int Standard_Model::m_ewscheme
    This is used to define the way, the electroweak parameters are fixed by some input; the default
    value of m_ewscheme, if not specified differently in the input data, is 0:
    -0 : all are given explicitly.
    -1 : \f$\alpha_{QED},\; M_W,\; M_Z\f$ and \f$M_H\f$
    -2 : \f$\alpha_{QED},\; \sin\theta_W,\; v_{vev}\f$ and \f$M_H\f$
    Note that \f$\alpha_{QED}\f$ is given in the Thomson limit.

    The default values of the parameters above are:
    \f$1/\alpha_{QED}=137.03599976,\f$, \f$M_W\f$, \f$M_H\f$,  and \f$M_Z\f$ as given by the particle 
    data file, \f$\sin^2\theta_W=0.23,\; v_{vev}=246\f$ GeV, and \f$\lambda_H=0.47591\f$.
    \todo Define & test more schemes
  */
  /*!
    \var int Standard_Model::m_ckmorder
    The order in the Cabbibo angle up to which the CKM matrix is evaluated. Note that the input
    parameters are given in the Wolfenstein parametrization
  */
  /*!
    \fn Standard_Model::Standard_Model(std::string,std::string)
    The constructor initializes the parameters of this class and the pointers to
    the searchable maps.
  */
  /*!
    \fn void Standard_Model::FillSpectrum()
    This method is called by the constructor and, through FixEWParameters and FixCKM defines
    all parameters in the electroweak sector before the Yukawa couplings of the fermions
    are read in and set in the corresponding flavour data. The running coupling \f$\alpha_s\f$
    is also initialized in this method.
  */
  /*!
    \fn void Standard_Model::FixEWParameters()
    According to the value of Standard_Model::m_ewscheme in this method the electroweak 
    parameters are set. These include: The masses of the gauge bosons and the Higgs boson,
    \f$M_W,\; M_Z\f$, and \f$M_H\f$, the electroweak (Weinberg-) angle, i.e.
    \f$\sin^2\theta_W\f$ and \f$\cos^2\theta_W\f$, the vacuum expectation value of the
    Higgs potential, \f$v_{vev}\f$, and the Higgs self-coupling, \f$\lambda_H\f$. Also,
    the Fermi constant \f$G_F\f$ is calculated and the running QED coupling is initialised.
  */
  /*!
    \fn void Standard_Model::FixCKM()
    In this method the four parameters of the CKM matrix in the Wolfenstein parametrization
    are read in and the CKM matrix is set up to a certain order in the Cabbibo-angle, as
    given by m_ckmorder, per default zero, i.e. a diagonal CKM matrix. The default values for the 
    paraemters are \f$\lambda=0.22,\; A=0.8,\; \eta=\rho=0.5\f$.
  */
}

#endif


