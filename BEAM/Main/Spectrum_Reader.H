#ifndef BEAM_Main_Spectrum_Reader_H
#define BEAM_Main_Spectrum_Reader_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace ATOOLS {
  class Function_Base;
}

namespace BEAM {
  class Spectrum_Point {
    double m_x1,m_x2,m_y,m_pol;
    unsigned int m_id;
  public:
    Spectrum_Point(double x1,double x2,double y, double pol):
      m_x1(x1), m_x2(x2), m_y(y), m_pol(pol), m_id((unsigned int)-1) {};
    Spectrum_Point(double x1,double x2, unsigned int id) :
      m_x1(x1), m_x2(x2), m_y(0.), m_pol(0.), m_id(id) {};
    double X1() const  { return m_x1; }
    double X2() const  { return m_x2; }
    double Y()  const  { return m_y; }
    double Pol() const  { return m_pol; }
    unsigned int Id() const { return m_id; }
  };

  typedef  std::vector<ATOOLS::Function_Base *> Spectrum_Function_List;
  typedef  std::vector<Spectrum_Point> Spectrum_Data_List;

  class Spectrum_Reader : public Beam_Base {
    std::string m_fname;
    Spectrum_Function_List m_spectrum_funcs;
    Spectrum_Data_List m_spectrum_histo;
    void ReadFromFile();
    ATOOLS::Function_Base * GetFunc(MyStrStream &,const double );
    double m_energy_laser, m_polarisation_laser;
    int m_mode;
    double m_polar;
    double m_upper, m_peak;
    ATOOLS::Flavour  m_remnant;
    double m_x;
  public:
    Spectrum_Reader(const ATOOLS::Flavour,const double,const double,
		    const double, const double, const std::string fname,const int);
    Beam_Base * Copy();     

    ATOOLS::Vec4D   OutMomentum();
    ATOOLS::Flavour Remnant();

    bool        CalculateWeight(const double,const double);
    double      Weight(ATOOLS::Flavour=ATOOLS::Flavour(kf_none));
    double          Polarisation()       { return m_polar; }
    bool        On() { return 1; }
    double          Exponent()           { return 0.5; }
    double          Xmax()               { return m_upper; }
    double          Peak()               { return m_peak; }
    void            PrintSpectra(std::string,int=0);
  };
  /*!
    \file 
    \brief contains the class BEAM::Spectrum_Reader.
  */
  /*!
    \class Spectrum_Reader
    \brief This class is for as simple spectrum read from a file

    The idea is to give a simple interface for returning a in princple
    arbitrary spectrum given by the user. The main application will
    be photon photon laser backscatering, where realistic spectra 
    functions exist only for certain collider configuration (specific
    choice of energy of laser and electron beam).

    There are two main possibilities, to specify a spectrum.
     -# give a normalized histogram for the range 0. to 1.
     -# give a set of function and ranges, that piecewise 
        approximate the spectrum
    .
  */
  /*!
    \fn Spectrum_Reader::Spectrum_Reader(const ATOOLS::Flavour,const double,const double,const int)
    This is the constructor of the class. Vectors are constructed through a call to the 
    constructor of the Beam_Base.
  */
  /*!
    \fn Beam_Base * Spectrum_Reader::Copy()
    Returns an exact copy of this class.
  */  
  /*!
    \fn 
    bool Spectrum_Reader::CalculateWeight(const double,const double)
  */  
  /*!
    \fn double Spectrum_Reader::Weight(ATOOLS::Flavour=ATOOLS::Flavour(kf_none))
  */  
}
#endif
