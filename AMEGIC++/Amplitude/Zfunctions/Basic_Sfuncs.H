#ifndef AMEGIC_Amplitude_Zfunctions_Basic_Sfuncs_H
#define AMEGIC_Amplitude_Zfunctions_Basic_Sfuncs_H


#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Spinor.H"
#include "AMEGIC++/Main/Pol_Info.H"
#include "AMEGIC++/Amplitude/Pfunc.H"
#include <list>

#define masslessskip 30
#define massiveskip 60
 
namespace AMEGIC {

  //!Container for explicit 4-vectors (momentums and polarization vectors)
  class Momfunc {
  public:
    //!Number of arguments
    int       argnum;
    /*!
      List of Arguments.
      For external momentums:
      - [0] index of the external particle
      
      For propagator momentums:
      - [0] index in AMEGIC::Basic_Sfuncs::Momlist
      - [1]..[argnum] list of external particle indices for calculation

      For polarization vectors
      - [0] index in AMEGIC::Basic_Sfuncs::Momlist
      - [1] index of the corresponding momentum
    */
    int*      arg;
    //!real part of the vector
    ATOOLS::Vec4D     mom;
    //!imaginary part of the vector
    ATOOLS::Vec4D     mom_img;
    //! polarization angle for external linear polarized vector bosons 
    double    angle;
    //! mass of the particle, described by the vector (polarizations only)
    double    mass;
    //! type of the vector
    mt::momtype type;
    //! kf code
    kf_code kfc;
    Momfunc() {arg = 0;argnum = 0;angle=0.;mass=0.;kfc=0;}
    Momfunc(const Momfunc& m) {
      arg = 0;
      argnum = 0;
      *this = m;
    }

    ~Momfunc() {if (arg) delete[] arg;}
      

    Momfunc& operator=(const Momfunc& m) {
      if (this!=&m) {
	argnum = m.argnum;
	if (arg) delete[] arg;
	if (argnum>0) {
	  arg = new int[argnum];
	  for (short int i=0;i<argnum;i++) arg[i] = m.arg[i];
	}
	mom     = m.mom;
	mom_img = m.mom_img;
	mass    = m.mass;
	angle   = m.angle;
	type    = m.type;
	kfc     = m.kfc;
      }
      return *this;
    }
  };

  class Basic_Sfuncs {
    std::vector<Momfunc> Momlist;
    int  momcount;
    ATOOLS::Flavour* fl;
    int  nmom,nvec;
    ATOOLS::Vec4D* p, m_k1, m_k2, m_k3;
    std::vector<ATOOLS::Vec4D>* p_epol;
    int* b;
    Complex* _eta;
    Complex* _mu;
    Complex** _S0;
    Complex** _S1;
    int** calc_st;
    int  k0_n;
    int  m_precalc;
    int  InitializeMomlist();
    void CalcMomlist();
    void CalcS(int i, int j);
    void PrecalcS();
    inline static int R1() { return ATOOLS::Spinor<double>::R1(); }
    inline static int R2() { return ATOOLS::Spinor<double>::R2(); }
    inline static int R3() { return ATOOLS::Spinor<double>::R3(); }
    inline ATOOLS::Vec4D K1() const { return m_k1; }
    inline ATOOLS::Vec4D K2() const { return m_k2; }
    inline ATOOLS::Vec4D K3() const { return m_k3; }
    inline double C1(const ATOOLS::Vec4D &p) const { return -m_k1*p; }
    inline double C2(const ATOOLS::Vec4D &p) const { return -m_k2*p; }
    inline double C3(const ATOOLS::Vec4D &p) const { return -m_k3*p; }
  public:
    Basic_Sfuncs(int,int, ATOOLS::Flavour*,int*);
    Basic_Sfuncs(int,int, ATOOLS::Flavour*,int*,std::string,std::string);
    ~Basic_Sfuncs();
    void Output(std::string name);

    int    BuildMomlist(Pfunc_List&); 
    int    BuildTensorPolarisations(int); 
    int    BuildPolarisations(int,char,double angle=0.);
    int    BuildPolarisations(int,ATOOLS::Flavour);

    void   PrintMomlist();
    int    GetMomNumber(Pfunc*);
    void   PropPolarisation(int,Pfunc_List&,std::vector<int>&);
    int    GetPolNumber(int,int,double,int check=0);

    void   Setk0(int);
    ATOOLS::Vec4D Getk0();
    ATOOLS::Vec4D Getk1();
    void   Initialize();
    int    CalcEtaMu(ATOOLS::Vec4D*);//setS
    void   InitGaugeTest(double);//ResetS_GT
    void   SetEPol(std::vector<ATOOLS::Vec4D>* epol) {p_epol=epol;}

    inline Complex Mu(int i)  { return (i>0) ? _mu[i] : -_mu[-i]; }
    inline Complex Eta(int i) { return _eta[ATOOLS::iabs(i)]; }
    inline Complex S0d(int i,int j) { return _S0[ATOOLS::iabs(i)][ATOOLS::iabs(j)]; }
    inline Complex S1d(int i,int j) { return _S1[ATOOLS::iabs(i)][ATOOLS::iabs(j)]; }
    inline Complex S0(int i,int j) {
      i=ATOOLS::iabs(i);j=ATOOLS::iabs(j);
      if (!calc_st[i][j] && !m_precalc) CalcS(i,j); 
      return _S0[i][j];
    }
    inline Complex S1(int i,int j) {
      i=ATOOLS::iabs(i);j=ATOOLS::iabs(j);
      if (!calc_st[i][j] && !m_precalc) CalcS(i,j); 
      return _S1[i][j];
    }
    Complex CalcS(ATOOLS::Vec4D& m, ATOOLS::Vec4D& m1);
    std::pair<Complex, Complex> GetS(ATOOLS::Vec4D v, int i);

    double Norm(int,int);//N

    inline int  GetNmomenta() { return nmom; }
    inline int  MomlistSize() { return Momlist.size(); }
    inline int  Sign(int i)   { return b[ATOOLS::iabs(i)]; }
    inline ATOOLS::Flavour GetFlavour(int i)   { return fl[ATOOLS::iabs(i)]; }
    inline ATOOLS::Vec4D& Momentum(int i)     { return Momlist[i].mom; }
    inline ATOOLS::Vec4D& MomentumImg(int i) { return Momlist[i].mom_img; }
    inline bool IsComplex(const int i) {
      return ( Momlist[i].type==mt::p_p || Momlist[i].type==mt::p_m || 
	       Momlist[i].type==mt::p_l || Momlist[i].type==mt::p_s || Momlist[i].type==mt::p_si);
    } 
    void StartPrecalc();
    bool IsMomSum(int,int,int);
  };

  std::ostream& operator<<(std::ostream& os, const Momfunc& mf);
  std::istream& operator>>(std::istream& is, Momfunc& mf);


  /*!
    \class Basic_Sfuncs
    \brief Calculation of S-functions.

    This class generates a list of all four-momentums and polarization vectors 
    and determines them from the given four-momentums of external particles.
    It also calculates and administrates the \f$\eta\f$- and \f$\mu\f$-functions
    as well as the basic spinor products \f$S(+;p_1,p_2)\f$ and \f$S(-;p_1,p_2)\f$.
  */
  /*!
    \var std::vector<Momfunc> Basic_Sfuncs::Momlist
    List of all momentums and polarization vectors.
  */
  /*!
    \var ATOOLS::Flavour* Basic_Sfuncs::fl
    Array of Flavours for external particles.
  */
  /*!
    \var int Basic_Sfuncs::nmom
    Number of external particles
  */
  /*!
    \var int Basic_Sfuncs::nvec
    Numbers of vectors for external particles (nmom + extra vectors for the old gauge boson treatment)
  */
  /*!
    \var ATOOLS::Vec4D* Basic_Sfuncs::p
    Array of momentums for external particles
  */
  /*!
    \var int* Basic_Sfuncs::b
    Signs for the external particles:
    - +1 for outgoing
    - -1 for incoming
  */
  /*!
    \var Complex* Basic_Sfuncs::_eta
    Array of \f$\eta\f$'s for each vector in Momlist. 
    It is calculated in Basic_Sfuncs::CalcEtaMu(ATOOLS::Vec4D*).
  */
  /*!
    \var Complex* Basic_Sfuncs::_mu
    Array of \f$\mu\f$'s for each vector in Momlist. 
    It is calculated in Basic_Sfuncs::CalcEtaMu(ATOOLS::Vec4D*).
  */
  /*!
    \var Complex** Basic_Sfuncs::_S0
    Array of \f$S(+,p_i,p_j)\f$, calculated in Basic_Sfuncs::CalcS(int i, int j) when first used.
  */
  /*!
    \var Complex** Basic_Sfuncs::_S1
    Array of \f$S(-,p_i,p_j)\f$, calculated in Basic_Sfuncs::CalcS(int i, int j) when first used.
   */
  /*!
    \var int** Basic_Sfuncs::calc_st
    Array to keep track over already calculated S-functions.
   */
  /*!
    \var int  Basic_Sfuncs::momcount
    Number of vectors in Basic_Sfuncs::Momlist.
   */
  /*!
    \var int  Basic_Sfuncs::k0_n
    Indicates what set of vectors \f$k_0\f$ and \f$k_1\f$ is used to define the spinor basis.
    See Basic_Sfuncs::Setk0(int) for details.
  */

  /*!
    \fn int  Basic_Sfuncs::Basic_Sfuncs::InitializeMomlist()
    Initializes momentums of external particles in Basic_Sfuncs::Momlist.
  */
  /*!
    \fn void Basic_Sfuncs::CalcMomlist()
    Calculates all vectors in Basic_Sfuncs::Momlist.

    <table border>    
      <tr>
        <td> </td>
        <td> <B>type</B> (see class mt)</td>
	<td> </td>
      </tr>
     <tr>
        <td>Momentums and Popagators: </td>
        <td>mom</td>
        <td>\f$p_i\f$</td>
      </tr>
      <tr>
        <td> </td>
        <td>prop</td>
        <td>\f$\sum_i p_i\f$</td>
      </tr>
      <tr>
        <td> </td>
        <td>cmprop</td>
        <td>\f$p_0+p_1\f$</td>
      </tr>
      <tr>
        <td>Polarizations: </td>
        <td>p_m/p_p</td>
        <td>\f$\frac{1}{\sqrt{2}\sqrt{p_x^2+p_y^2}}\left(0,\frac{p_xp_z}{|\vec{p}|}\mp ip_y,
	       \frac{p_yp_z}{|\vec{p}|}\pm ip_x,-\frac{p_x^2+p_y^2}{|\vec{p}|}\right)\f$</td>
      </tr>
      <tr>
        <td> </td>
        <td>p_lh</td>
        <td>\f$\frac{1}{\sqrt{p_x^2+p_y^2}}\left(0,\frac{p_xp_z}{|\vec{p}|},
	       \frac{p_yp_z}{|\vec{p}|},-\frac{p_x^2+p_y^2}{|\vec{p}|}\right)\f$</td>
      </tr>
      <tr>
        <td> </td>
        <td>p_lv</td>
        <td>\f$\frac{1}{\sqrt{p_x^2+p_y^2}}\left(0,p_y,-p_x,0\right)\f$</td>
      </tr>
      <tr>
        <td> </td>
        <td>p_l</td>
        <td>\f$\frac{1}{\sqrt{p^2}}\left(|\vec{p}|,p_0\frac{\vec{p}}{\vec{p}|}\right)\f$</td>
      </tr>
      <tr>
        <td> </td>
        <td>p_s</td>
        <td>\f$\sqrt{\frac{p^2-(m^2+im\Gamma)}{p^2 (m^2+im\Gamma)}}p\f$</td>
      </tr>
       <tr>
        <td> </td>
        <td>p_si</td>
        <td>\f$\sqrt{-\frac{1}{p^2}}p\f$</td>
      </tr>
    </table>
  */
  /*!
    \fn void Basic_Sfuncs::CalcS(int i, int j)
    Calculates the basic spinor products \f$S(+,p_i,p_j)\f$ and \f$S(-,p_i,p_j)\f$. 
    i,j are indices of Basic_Sfuncs::Momlist. 

    S-Functions are defined as products of massless spinor:
    \f[
    S(+,p_1,p_2)=\bar{u}(p_1,+)u(p_2,-),
    \f] 
    \f[
    S(-,p_1,p_2)=\bar{u}(p_1,-)u(p_2,+).
    \f]
 
    We introduce a dimensionless chiral spinor basis 
    \f[
    w(k_0,\lambda)\bar w(k_0,\lambda) = \frac{1+\lambda\gamma_5}{2}k\!\!\!/_0
    \f] 
    and
    \f[
    w(k_0,\lambda) = \lambda k\!\!\!/_1 w(k_0,-\lambda)\,,
    \f] 
    where \f$k_1 k_1=-1\f$ and \f$k_0 k_1=0\f$.
    
    Now a (massless) spinor with an arbitrary momentum \f$p\f$ can be defined by
    \f[
    u(p,\lambda)=\frac{p\!\!\!/}{\sqrt{2pk_0}}w(k_0,-\lambda).
    \f]

    Using this the S-Function can be calculated:
    \f[
    S(+;p_1,p_2) = 
    2\frac{(p_1\cdot k_0)(p_2\cdot k_1)-
           (p_1\cdot k_1)(p_2\cdot k_0)+
           i\epsilon_{\mu\nu\rho\sigma}
           p_1^\mu p_2^\nu k_0^\rho k_1^\sigma}{\eta_1\eta_2}\,, 
    \f]
    \f[
    S(-;p_1,p_2) = -S(+;p_1,p_2)^*\,.
    \f]

    With the default choice of \f$k_0\f$ and \f$k_1\f$ (see Basic_Sfuncs::Setk0), that 
    expression can be simplified to
    \f$
    S(\pm,p_1,p_2)=\left(\pm p_1^y+\frac{i}{\sqrt{2}}(p_1^z-p_1^x)\right)
                   \frac{\eta_2}{\eta_1} - (1 \leftrightarrow 2)\,.
    \f$
  */
  /*!
    \fn void Basic_Sfuncs::PrecalcS()
    Precalculation of all used spinor products. 
    This method is used after calling Basic_Sfuncs::StartPrecalc().
  */
 
  /*!
    \fn int    Basic_Sfuncs::BuildMomlist(Pfunc_List&)
    Adds propagators and internal polarization vectors to Basic_Sfuncs::Momlist.
  */
  /*!
    \fn int    Basic_Sfuncs::BuildTensorPolarisations(int)
    Adds polarization vectors, neccessary to calculate polarization tensors f
    or external spin-2 particles to Basic_Sfuncs::Momlist.
  */
  /*!
    \fn int    Basic_Sfuncs::BuildPolarisations(int,char,double angle=0.)
    Initializes polarizations for external vector bosons in Basic_Sfuncs::Momlist.
  */
  /*!
    \fn int    Basic_Sfuncs::BuildPolarisations(int,ATOOLS::Flavour)
    Initializes polarizations for vector boson and spin-2 propagators in Basic_Sfuncs::Momlist.
  */

  /*!
    \fn void   Basic_Sfuncs::PrintMomlist()
    Prints Basic_Sfuncs::Momlist.
  */
  /*!
    \fn int    Basic_Sfuncs::GetMomNumber(Pfunc*)
    Checks Basic_Sfuncs::Momlist if a propagator already exists. 
    Returns the index or -1.
  */
  /*!
    \fn void   Basic_Sfuncs::PropPolarisation(int,Pfunc_List&,std::vector<int>&)
    In the last argument a list of corresponding polarizaton vector types (mt) for a 
    propagator is returned.
  */
  /*!
    \fn int    Basic_Sfuncs::GetPolNumber(int momindex,int sign, double mass,int check=0)
    Checks Basic_Sfuncs::Momlist if a polarization vector already exists.
    Polarizations are characterized a propagator or external momentum , a type (mt)
    and a mass of the corresponding boson if the calculation depends on it.
    Returns the index or -1.
  */
  /*!
    \fn void   Basic_Sfuncs::Setk0(int)
    Set the auxiliary vectors \f$k_0\f$ and \f$k_1\f$, defining the spinor basis.

    The following sets are available:
    - 0:  \f$k_0=(1,\sqrt(1/2),0,\sqrt(1/2))\;\;\;\;\;k_1=(0,0,1,0)\f$  (default)
    - 1:  \f$k_0=(1,0,\sqrt(1/2),\sqrt(1/2))\;\;\;\;\;k_1=(0,1,0,0)\f$
    - 2:  \f$k_0=(1,\sqrt(1/2),\sqrt(1/2),0)\;\;\;\;\;k_1=(0,0,0,1)\f$

    This function can be used for internal tests.
  */
  /*!
    \fn ATOOLS::Vec4D Getk0()
    Returns the currently set k0 as a 4-vector.
  */
  /*!
    \fn void   Basic_Sfuncs::Initialize()
    Memory allocation for Basic_Sfuncs::_mu, Basic_Sfuncs::_eta, Basic_Sfuncs::_S0,
    Basic_Sfuncs::_S1 and Basic_Sfuncs::calc_st.
  */
  /*!
    \fn int    Basic_Sfuncs::CalcEtaMu(ATOOLS::Vec4D*)
    Calculation of \f$\eta\f$ and \f$\mu\f$ for every vector in Basic_Sfuncs::Momlist.

    \f$\eta_i=\sqrt{2p_ik_0}\f$ 

    \f$\mu_i=\pm\frac{m_i}{\eta_i}\f$ for external particles/antiparticles, and
    
    \f$\mu_i=\frac{\sqrt{p_i}}{\eta_i}\f$ for propagators and polarizations. 

    If Basic_Sfuncs::m_precalc is set, all basic spinor products, enabled in 
    Basic_Sfuncs::calc_st are calculated. 
  */
  /*!
    \fn void   Basic_Sfuncs::InitGaugeTest(double theta)
    Initializes the internal gauge test for external massless vector bosons (photons or gluons).
    The test can only performed for unpolarized beams with a circular polarization basis.

    External massless gauge bosons obey the completeness relation in the light-cone gauge
    \f[
    \sum_{\lambda=\pm}\epsilon_{\mu}(p,\lambda)\epsilon^{*}_{\nu}(p,\lambda) =
    -\eta_{\mu\nu} + \frac{p_\mu q_\nu+p_\nu q_\mu}{p\cdot q}
    \f]
    where \f$q\f$ is any four-vector not alligned with \f$p\f$ and \f$q^2=0\f$.
    
    By default \f$q\f$ is choosen to be \f$q=(p_0,-\vec{p})\f$, where
    \f[
    \epsilon(p,\pm) = \frac{1}{\sqrt{2}}\left(0,
    \cos\theta\cos\varphi \mp i\sin\varphi,\cos\theta\sin\varphi \pm i\cos\varphi,
    -\sin\theta \right)    
    \f]
    is a possible choice for the polarization vectors of a particle with momentum
    \f$p_{\mu} = \left(p_0,p_0\sin\theta\cos\varphi,
               p_0\sin\theta\sin\varphi,p_0\cos\theta\right)\f$.

    This function redefines the polarization vectors to yield a completeness relation, 
    where spacial part of \f$q\f$ has an arbitrary angle \f$\theta_0\f$ with respect to \f$-\vec{p}\f$.
    The circular polarization vectors are now
    \f[
    \epsilon(p,\pm)=\frac{1}{{\cal N}}\times
    \left[
    \left(
    \begin{array}{l}
    \cos{\frac{\theta_0}{2}}\cos{\frac{\theta}{2}}+\sin{\frac{\theta_0}{2}}\sin{\frac{\theta}{2}}\cos\varphi
    \\
    \sin{\frac{\theta_0}{2}}\cos{\frac{\theta}{2}}+\cos{\frac{\theta_0}{2}}\sin{\frac{\theta}{2}}\cos\varphi
    \\
    \cos{\frac{\theta_0}{2}}\sin{\frac{\theta}{2}}\sin\varphi
    \\
    \cos{\frac{\theta_0}{2}}\cos{\frac{\theta}{2}}-\sin{\frac{\theta_0}{2}}\sin{\frac{\theta}{2}}\cos\varphi
    \end{array}
    \right) 
    \mp i
    \left(
    \begin{array}{l}
    \sin{\frac{\theta_0}{2}}\sin{\frac{\theta}{2}}\sin\varphi                                         
    \\
    \cos{\frac{\theta_0}{2}}\sin{\frac{\theta}{2}}\sin\varphi
    \\
    \sin{\frac{\theta_0}{2}}\cos{\frac{\theta}{2}}-
    \cos{\frac{\theta_0}{2}}\sin{\frac{\theta}{2}}\cos\varphi
    \\
    -\sin{\frac{\theta_0}{2}}\sin{\frac{\theta}{2}}\sin\varphi
    \end{array}
    \right)\right]
    \f]
    where
    \f[
    {\cal N} = \sqrt{1-\sin\theta\cos\varphi\sin{\theta_0}-\cos\theta\cos{\theta_0}}\,.
    \f]
  */
  /*!
    \fn inline Complex Basic_Sfuncs::Mu(int i) 
    Returns \f$\mu_i\f$.
  */
  /*!
    \fn inline Complex Basic_Sfuncs::Eta(int i) 
    Returns \f$\eta_i\f$.
  */
  /*!
    \fn inline Complex Basic_Sfuncs::S0(int i,int j) 
    Checks Basic_Sfuncs::calc_st if \f$S(+,p_i,p_j)\f$ is already calculated.
    If not Basic_Sfuncs::CalcS is started.
    Returns the complex value. 
  */
  /*!
    \fn inline Complex Basic_Sfuncs::S1(int i,int j) 
    Checks Basic_Sfuncs::calc_st if \f$S(-,p_i,p_j)\f$ is already calculated.
    If not Basic_Sfuncs::CalcS is started.
    Returns the complex value. 
  */
  /*!
    \fn inline Complex Basic_Sfuncs::S0d(int i,int j) 
    Returns \f$S(+,p_i,p_j)\f$. Basic_Sfuncs::StartPrecalc must be performed first.
  */
  /*!
    \fn inline Complex Basic_Sfuncs::S1d(int i,int j) 
    Returns \f$S(-,p_i,p_j)\f$. Basic_Sfuncs::StartPrecalc must be performed first.
  */
  /*!
    \fn std::pair<Complex, Complex> Basic_Sfunc::GetS(ATOOLS::Vec4D v, int i)
    Returns the results of both S-functions S+ and S- of the vectors v and the i-th
    vector in the momlist. S+ is stored in the first entry of the pair, S- is stored
    in the second entry.
  */
  /*!
    \fn double Basic_Sfuncs::Norm(int,int)
  */

  /*!
    \fn int  Basic_Sfuncs::GetNmomenta()
    Returns the number of external particles.
  */
  /*!
    \fn inline int  Basic_Sfuncs::Sign(int i)
    Returns a sign for external particles (see Basic_Sfuncs::b).
  */
  /*!
    \fn inline ATOOLS::Flavour  Basic_Sfuncs::GetFlavour(int i)
    Returns the flavour of an external particle.
  */
  /*!
    \fn inline ATOOLS::Vec4D  Basic_Sfuncs::Momentum(int i)
    Returns the real part of a momentum in Basic_Sfuncs::Momlist.
  */
  /*!
    \fn inline ATOOLS::Vec4D  Basic_Sfuncs::MomentumImg(int i) 
    Returns the imaginary part of a momentum in Basic_Sfuncs::Momlist.
  */
  /*!
    \fn inline bool Basic_Sfuncs::IsComplex(int i)
    Returns true, if the vector in Basic_Sfuncs::Momlist[i] may become complex.
    This is the case for mt::p_m, mt::p_p, mt::p_l, mt::p_s and mt::p_si.
  */
  /*!
    \fn void Basic_Sfuncs::StartPrecalc()
    Activates the precalculation modus for the basic spinor products.
    All used spinor products are now calculated in advance.
  */

  /*!
    \fn bool Basic_Sfuncs::IsMomSum(int x,int y,int z) 
    Checks if the momentum in Basic_Sfuncs::Momlist[x] is a sum of the other momentums
    or proportional to it.
    The method is used to check if some Basic_Xfunc are identical zero.
  */
}
#endif
