#ifndef HADRONS_PS_Library_ResonanceFlavour_H
#define HADRONS_PS_Library_ResonanceFlavour_H

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Histogram.H"

namespace HADRONS {
   
  

  class SimpleResonanceFlavour {

    protected:
      std::string         m_name;
      double              m_mass, m_width, m_mass2;

    public:
      SimpleResonanceFlavour() {}
      SimpleResonanceFlavour( std::string name, double _mass, double _width );
      ~SimpleResonanceFlavour() {}
       
      double      Mass() { return m_mass; }
      double      Mass2() { return m_mass2; }
      double      Width() { return m_width; }
      double      MassWidth() { return m_mass*m_width; }
      std::string Name() { return m_name; }
  };

  class ResonanceFlavour : public SimpleResonanceFlavour {

    protected:
      kf_code    m_kfc;
      int                 m_running;
      std::string         m_path;
      ATOOLS::Histogram * p_hist;
      double              m_G_at_m2;
      int                 m_body;

      double TwoBodyResonanceMassWidth( double s, double m1, double m2 ); 
      double TwoBodyResonanceMassWidth( double s, double m ); 
      double TwoBodyResonanceMassWidth_12( double s, double m ); 
      double ThreeBodyResonanceMassWidth( double s );

      ATOOLS::Histogram * CreateGHistogram( ResonanceFlavour res1, ResonanceFlavour res2, double beta, kf_code out ); 
      double      GetValueOfG( double s ); 
      double      IntegralG( double Q2, ResonanceFlavour res1, ResonanceFlavour res2, double beta, kf_code out );
     
    public:
      static double Lambda( double a, double b, double c ); 
      static double Sqrt_Lambda( double a, double b, double c); 
	public:
      ResonanceFlavour() {}
      ResonanceFlavour( kf_code _kfc, double _mass, double _width, int _run, std::string path="" );
      ~ResonanceFlavour() { if(p_hist) delete p_hist; }

      void        InitialiseThreeBodyResonance( ResonanceFlavour &res1, ResonanceFlavour &res2, double beta=0. );
      void        InitialiseThreeBodyResonance( ResonanceFlavour &res1 );

      Complex     BreitWigner( double s ); 
      Complex     BreitWignerAlt( double s ); 
      double      OffShellMassWidth( double s ); 
      int         Running() { return m_running; }

  };

}

#endif
  
