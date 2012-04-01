#ifndef SHRIMPS_Tools_MinBias_Parameters_H
#define SHRIMPS_Tools_MinBias_Parameters_H

#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <map>
#include <string>




namespace SHRIMPS {
  struct run_mode {
    enum code {
      xsecs_only                = 0,
      elastic_events            = 1,
      single_diffractive_events = 2,
      double_diffractive_events = 3,
      quasi_elastic_events      = 4,
      inelastic_events          = 5,
      all_min_bias              = 10,
      underlying_event          = 20,
      unknown                   = -1
    };
  };
  std::ostream & operator<<(std::ostream & s, 
			    const run_mode::code & runmode);

  struct weight_mode {
    enum code {
      unweighted = 0,
      weighted   = 1,
      unknown    = -1
    };
  };
  std::ostream & operator<<(std::ostream & s, 
			    const weight_mode::code & weightmode);

  /*! 
    \struct ff_form
    \brief A struct used to select the form of the form factor in the class
    SHRIMPS::Form_Factor.  Available options up to now are a Gaussian form and
    a modified dipole form, supplemented with an additional exponent factor 
    to guarantee convergence of the Fourier transform.
  */
  struct ff_form {
    enum code {
      Gauss  = 2,
      dipole = 1
    };
  };
  std::ostream & operator<<(std::ostream & s, 
			    const ff_form::code & ffform);

  struct absorption {
    enum code {
      exponential = 2,
      factorial   = 1
    };
  };
  std::ostream & operator<<(std::ostream & s, 
			    const absorption::code & absorb);

  struct ladder_weight {
    enum code {
      IntervalOnly,
      ReggeDiffusion,
      Regge
    };
  };
  std::ostream & operator<<(std::ostream & s, 
			    const ladder_weight::code & wt);





  struct ktform {
    enum code {
      IR0      = 2,
      smooth = 1,
      frozen = 0,
      cut    = -1
    };
  };
  struct ordering {
    enum code {
      ao_phys  = 6,
      rap_phys = 4,
      ao_keep  = 3,
      ao       = 2,
      keep     = 1,
      rap_only = 0
    };
  };
  struct resc_ktmin {
    enum code {
      on  = 1,
      off = 0
    };
  };
  struct resc_nosing {
    enum code {
      on  = 1,
      off = 0
    };
  };
  struct resc_over_sing {
    enum code {
      on  = 1,
      off = 0
    };
  };
  struct resc_mode {
    enum code {
      on  = 1,
      off = 0
    };
  };
  struct spat_mode {
    enum code {
      expo = 1,
      off  = 0
    };
  };
  struct reconn_mode {
    enum code {
      run = 2,
      fix = 1,
      off = 0
    };
  };


  typedef std::map<std::string,double> StringDoubleMap;
  typedef StringDoubleMap::iterator    SDiter;

  class MinBias_Parameters {
  private:
    run_mode::code       m_runmode;
    weight_mode::code    m_weightmode;
    ff_form::code        m_ffform;
    absorption::code     m_absorp;
    ladder_weight::code  m_ladderweight;
    ktform::code         m_ktform;
    ordering::code       m_ordering;
    resc_ktmin::code     m_resc_ktmin;
    resc_nosing::code    m_resc_nosing;
    resc_over_sing::code m_resc_over_sing;
    resc_mode::code      m_rescmode;
    spat_mode::code      m_spacemode;
    reconn_mode::code    m_reconnmode;
    StringDoubleMap      m_params;
    MODEL::asform::code  m_as_form;
  public:
    MinBias_Parameters();
    ~MinBias_Parameters();

    void   Init(ATOOLS::Data_Reader * dr);
    double operator()(std::string keyword);
    void   Set(std::string keyword, double value) {
      if (m_params.find(keyword)!=m_params.end()) m_params[keyword] = value;
      else {
	std::cerr<<"Error in MinBias_Parameters::Set: '"
		 <<keyword<<"' not found.\n";
      }
    }

    run_mode::code       RunMode()      const { return m_runmode; }
    weight_mode::code    WeightMode()   const { return m_weightmode; }
    ff_form::code        FF_Form()      const { return m_ffform; }  
    absorption::code     Absorption()   const { return m_absorp; }  
    ladder_weight::code  LadderWeight() const { return m_ladderweight; }
    ktform::code         KTForm()       const { return m_ktform; }
    ordering::code       Ordering()     const { return m_ordering; }
    resc_mode::code      RescMode()     const { return m_rescmode; }
    spat_mode::code      SpaceMode()    const { return m_spacemode; }
    resc_ktmin::code     RescKTMin()    const { return m_resc_ktmin; }
    resc_nosing::code    RescNoSing()   const { return m_resc_nosing; }
    resc_over_sing::code RescOverSing() const { return m_resc_over_sing; }
    reconn_mode::code    ReconnMode()   const { return m_reconnmode; }
    MODEL::asform::code  AsForm()       const { return m_as_form; }
  };
  
  extern MinBias_Parameters MBpars;
}

#endif
