#ifndef Laser_Backscattering_H
#define Laser_Backscattering_H

#include "Beam_Base.H"

namespace BEAM {  
  class Laser_Backscattering : public Beam_Base {
    /*!
      The mode for the spectrum:
      0 = all contributions.
      1 = direct Compton process only.
      2 = scattering with 2 photons only.
      3 = scattering off secondary electrons.
    */
    int mode;
    /*!
      Collinear photons (0) or angles taken into account (1)
    */
    bool angles;
    /*!
      Whether degree of polarization should be calculated or not. 
    */
    bool pol;
    /*! 
      The electron energies bounds for this approximation to be resaonable.
      At the moment they are set to 50 GeV < E < 500 GeV.
    */
    double E,Ebounds[2];
    /*!
      The energy of the laser-photons.
    */
    double omegaL;
    double polE,polL;
    double rho2,delta;
    double nonlin1,nonlin2,xi;
    /*!
      The maximal energy fraction the backscattered photons can obtain from Compton
      scattering :          xmax. With xe = 4E*omega_L/m_e^2 it is xe/(xe+1)
      The maximal energy fraction for backscattered photons from the two photon process.
      It is given by  xmax2 = 2.*xe/(1.+2.*xe).
    */
    double xe, xmax, xmax2, Emax;
    double upper,peak;

    double yfix,yden;
    int    ysteps;

    double totalC,total2,totalE;

    double polar;

    double Compton(double,double,double,double &);
    double TwoPhotons(double,double,double,double &);
    double Rescattering(double,double,double,double &);
    double SimpleCompton(double x,double z,double pol);
    double Polarization(double x,double z,double pole,double poll);
  public:
    Laser_Backscattering(APHYTOOLS::Flavour,double);
    ~Laser_Backscattering();
    void   PrintSpectra(std::string);

    bool   CalculateWeight(double,double);
    double Weight(APHYTOOLS::Flavour);


    //! Simple access
    double PolarisationDegree()  { return polar; }

    void      SetMode(int _mode) { mode = _mode; }
    void      SetPol(bool on)    { pol  = on; }

    bool      On()               { return 1;  }
    int       Mode()             { return mode; }
    bool      Pol()              { return pol; }

    double    Exponent()         { return 2.; }
    double    Xmax()             { return upper; }
    double    Peak()             { return peak; }
  };
}


#endif
