#ifndef BEAM_Main_Laser_Backscattering_H
#define BEAM_Main_Laser_Backscattering_H

#include "BEAM/Main/Beam_Base.H"

namespace BEAM {  
  class Laser_Backscattering : public Beam_Base {
    double m_energyL,m_polarisationL;
    int    m_mode;
    bool   m_angles,m_pol;
    double m_Ebounds[2];
    double m_rho2,m_delta;
    double m_nonlin1,m_nonlin2,m_xi;
    double m_xe, m_xmax, m_xmax2;
    double m_upper,m_peak;

    double m_yfix,m_yden;
    int    m_ysteps;

    double m_totalC,m_total2,m_totalE;

    double m_polar;

    double Compton(double,double,double,double &);
    double TwoPhotons(double,double,double,double &);
    double Rescattering(double,double,double,double &);
    double SimpleCompton(double,double,double);
    double Polarisation(double,double,double,double);
  public:
    Laser_Backscattering(const ATOOLS::Flavour,const double,const double,
			 const double,const double,const int,const int,const int,const int);
    ~Laser_Backscattering();

    Beam_Base     * Copy();
    void            PrintSpectra(std::string,int=0);
    bool            CalculateWeight(double,double);
    double          Weight(ATOOLS::Flavour=ATOOLS::Flavour(kf_none));
    ATOOLS::Vec4D   OutMomentum();

    double          Polarisation()       { return m_polar; }
    void            SetMode(int _mode)   { m_mode = _mode; }
    void            SetPol(bool on)      { m_pol  = on; }
    bool            On()                 { return 1;  }
    int             Mode()               { return m_mode; }
    bool            Pol()                { return m_pol; }
    double          Exponent()           { return 0.5; }
    double          Xmax()               { return m_upper; }
    double          Peak()               { return m_peak; }
    bool            PolarisationOn() ;
  };


  /*!
    \file 
    \brief contains the class BEAM::Laser_Backscattering.
  */
  /*!
    \class Laser_Backscattering
    \brief This class is for Laser back scattering off electrons.

    In this class the 
    <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?rawcmd=FIND+eprint+hep-ex%2F0207021">
    CompAZ</A> parametrization is implemented. This is a phenomenological parametrization, including 
    three components, each of which is modified with some additional non linear terms. These 
    components are:
    - Compton scattering,
    - Rescattering processes, i.e. interactions involving two electrons, and
    - Two photon processes.

    They are all based on a simple Compton backscattering spectrum implemented in
    Laser_Backscattering::SimpleCompton(double,double,double).
  */
  /*!
    \var int Laser_Backscattering::m_mode
    The mode for the CompAZ spectrum:
    - 0 all contributions.
    - 1 direct contribution of Compton process only.
    - 2 scattering with 2 photons only.
    - 3 scattering off secondary electrons.
    .
    another option is to use just the simple compton spectrum
    - -1 for comparison only.
    .
  */
  /*!
    \var bool Laser_Backscattering::m_angles
    Collinear photons (0) or angles taken into account (1). At the moment  only option 0,
    i.e. collinear photons, is supported. 

    \todo Enable m_angles = 1.
  */
  /*!
    \var bool Laser_Backscattering::m_pol
    Whether degree of polarization should be calculated or not. This happens, if the incoming 
    electron beam has a polarization degree different from 0 or if the Laser photons are polarized.
  */
  /*! 
    \var double Laser_Backscattering::m_Ebounds[2]
    The electron energies bounds for this parametrization to be resaonable.
    At the moment they are set to 50 GeV < E < 500 GeV.
  */
  /*!
    \var double Laser_Backscattering::m_energyL
    The energy of the laser-photons, \f$\omega_L\f$.
  */
  /*!
    \var double Laser_Backscattering::m_polarisationL
    The polarization degree of the laser-photons, \f$\lambda_L\f$.
  */
  /*!
    \var double Laser_Backscattering::m_xe
    The characteristic parameter of the process, given by:
    \f[
    \xi = \frac{4E_e\omega_L}{m_e^2}\,.
    \f]
    If nonlinear corrections are included, this is modified to read
    \f[
    \xi = \frac{4E_e\omega_L}{m_e^2(1+\chi)}\,.
    \f]
  */
  /*!
    \var double Laser_Backscattering::m_xmax
    The maximal energy fraction \f$ x_{\rm max}\f$ the backscattered photons can obtain from Compton
    scattering. It is given by
    \f[
    x_{\rm max} = \frac{\xi}{1+\xi}
    \f]
  */
  /*!
    \var double Laser_Backscattering::m_xmax2
    The maximal energy fraction for backscattered photons from the two photon process, 
    \f$ x_{\rm max}^{(2)}\f$.
    It is given by  
    \f[
    x_{\rm max}^{(2)} = \frac{2\xi}{1+2\xi}\,.
    \f]
  */
  /*!
    \var double Laser_Backscattering::m_rho2
    A parameter for the damping of the Compton pieces of the spectrum, \f$\rho^2\f$.
  */
  /*!
    \var double Laser_Backscattering::m_delta
    A parameter for the damping of the two photon piece of the spectrum, \f$\delta\f$.
  */
  /*!
    \var double Laser_Backscattering::m_nonlin1
    Parameter for nonlinear corrections \f$\nu_1\f$.
  */
  /*!
    \var double Laser_Backscattering::m_nonlin2  
    Parameter for nonlinear corrections \f$\nu_2\f$.
  */
  /*!
    \var double Laser_Backscattering::m_xi
    The effective nonlinearity parameter \f$\chi\f$. It is given by
    \f[
    \chi = \nu_1 + \nu_2 E_e\,.
    \f]
  */
  /*! 
    \var double Laser_Backscattering::m_upper
    Depending on the mode, this is either \f$ x_{\rm max}\f$ or \f$ x_{\rm max}^{(2)}\f$.
  */
  /*!
    \var double Laser_Backscattering::m_peak
    The peak position of the spectrum, \f$ x_{\rm max}\f$.
  */
  /*!
    \var double Laser_Backscattering::m_yfix
    A parameter for the evaluation of the rescattering-convolution, \f$y_{\rm fix} = 1/(1-\xi)\f$.
  */
  /*!
    \var double Laser_Backscattering::m_yden
    A parameter for the evaluation of the rescattering-convolution, \f$y_{\rm den} = \log(1-\xi)\f$.
  */
  /*!
    \var int Laser_Backscattering::m_ysteps
    The number of steps in the convolution for the rescattering.
  */
  /*!
    \var double Laser_Backscattering::m_totalC
    The norm of the Compton contribution.
  */
  /*!
    \var double Laser_Backscattering::m_total2
    The norm of the two-photon contribution.
  */
  /*!
    \var double Laser_Backscattering::m_totalE
    The norm of the rescatering contribution.
  */
  /*!
    \var double Laser_Backscattering::m_polar
    The polarization degree.
  */
  /*!
    \fn Laser_Backscattering::Laser_Backscattering(const ATOOLS::Flavour,const double,const double,
                                                   const double,const double,const int,
						   const int,const int,bool &);
    This is the constructor of the class. It initializes a large variety of internal parameters.
    The vectors are constructed by calling the default constructor of the Beam_Base.
  */
  /*!
    \fn Beam_Base * Laser_Backscattering::Copy()
    A method to fully copy the class, including all parameter settings.
  */
  /*!
    \fn void Laser_Backscattering::PrintSpectra(std::string)
    When called, this method outputs the spectra and its contributions to a file specified by
    the string variable.
  */
  /*!
    \fn bool Laser_Backscattering::CalculateWeight(double,double);
    Calculates the weight dependent on hte energy fraction \f$x\f$ of the photon w.r.t the
    incoming lepton and in accordance with the mode given by Laser_Backscattering::m_mode.
  */
  /*!
    \fn double Laser_Backscattering::Weight(ATOOLS::Flavour=ATOOLS::Flavour(kf_none));
    After checking that the outgoing particle indeed is a photon, the calculated weight
    is returned. In principle, from there one could also return an electron energy distribution,
    but this has not been implmeneted yet.
  */
  /*!
    \fn double Laser_Backscattering::SimpleCompton(double,double,double).
    Evaluates and returns the normalized simple Compton backscattering spectrum
    \f[
    {\cal F}_C(x,z,\Lambda) = \frac{1}{{\cal N}(z,\Lambda)}\times
                            \left[1-x  + \frac{1}{1-x} - \frac{4z}{1+z} + \frac{4z^2}{(1+z)^2}
                            -\Lambda \frac{2x-x^2}{1-x}\cdot\left(\frac{2z}{1+z}-1\right)\right]\,,
    \f]
    where the norm reads
    \f[
    {\cal N}_C(z,\Lambda) = \frac{z^3+18z^2+32z+16}{2z(1+z)^2}
                          +\frac{z^2-4z-8}{z^2}\log(1+z)
			  -\Lambda\left[2+\frac{z^2}{(1+z)^2}-\frac{2+z}{z}\log(1+z)\right]\,,
    \f]
    for \f$ x\in [0,z]\f$. In most theory papers this form alone is used, then \f$ z = \xi\f$.
  */
  /*!
    \var double Laser_Backscattering::Compton(double,double,double,double &)
    This is the - eventually modified - Compton piece of the spectrum. It is given
    by the simple Compton piece above and some corrections and reads
    \f[
    {\cal F}_1(x,\lambda_e,\lambda_L,{\cal P}) = \frac{1}{{\cal N}_1}\cdot 
                   \exp\left[-\frac{\rho^2}{8}\left(\frac{\xi}{x}-\xi-1\right)\right]
                   \cdot {\cal F}_C(x,\xi,\lambda_e\cdot\lambda_L)\,.
    \f]
    The range for \f$ x\f$ is \f$ x\in [0,x_{\rm max}]\f$.
  */
  /*!
    \var double Laser_Backscattering::TwoPhotons(double,double,double,double &)
    This is the - eventually modified - two photon piece of the spectrum. It is given
    by the simple Compton piece above and some corrections and reads
    \f[
    {\cal F}_2(x,\lambda_e,\lambda_L,{\cal P}) = \frac{1}{{\cal N}_2}\cdot 
                   \exp\left[-\frac{\rho^2}{8}\left(\frac{2\xi}{x}-2\xi-1\right)\right]
		   \cdot \left(\frac{2\xi}{x}-2\xi-1\right)^\delta
                   \cdot {\cal F}_C(x,2\xi,\lambda_e\cdot\lambda_L)\,.
    \f]
    The range for \f$ x\f$ is \f$ x\in [0,x_{\rm max}^{(2)}]\f$.
  */
  /*!
    \var double Laser_Backscattering::Rescattering(double,double,double,double &)
    For the rescattering piece of the spectrum a convolution of two simple Compton
    spectra has to be performed. 
    The range for \f$ x\f$ is \f$ x\in [0,x_{\rm max}]\f$.
  */
  /*!
    \var double Laser_Backscattering::Polarization(double,double,double,double)
    The polarization piece of the spectrum yielding the degree of polarization for the outgoing
    photon. It is given by
    \f[
    {\cal P}_C(x,z,\lambda_e,\lambda_L) = \frac{1}{\cal N}_P
    \lambda_e \left[\frac{x}{1-x}\left(1+(1-x)\sqrt{\frac{2x}{z(1-x)}-1}\right)\right] -
    \lambda_L \left[\left(1-x+\frac{1}{1-x}\right)\left(\frac{2x}{z(1-x)}-1\right)\right]\,,
    \f]
    where the norm is given by
    \f[
    {\cal N}_P = 1-x+\frac{1}{1-x}+\frac{4x}{z(1-x)}\left(\frac{x}{z(1-x)}-1\right)
                 -\lambda_e\lambda_L\frac{x(2-x)}{1-x}\left(\frac{2x}{z(1-x)}-1\right)\,.
    \f]
    This expression is valid only for 
    \f[
    x \in \left[0,\frac{z}{1+z}\right]\,.
    \f]
  */
}


#endif
