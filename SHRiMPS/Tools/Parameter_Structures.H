#ifndef SHRIMPS_Tools_Parameter_Structures_H
#define SHRIMPS_Tools_Parameter_Structures_H

#include <ostream>

namespace SHRIMPS {
  struct run_mode {
    enum code {
      test                      = 0,
      xsecs_only                = 1,
      inelastic_events          = 10,
      elastic_events            = 11,
      soft_diffractive_events   = 12,
      single_diffractive_events = 13,
      double_diffractive_events = 14,
      quasi_elastic_events      = 15,
      all_min_bias              = 20,
      underlying_event          = 100,
      unknown                   = -1
    };
  };

  struct weight_mode {
    enum code {
      unweighted = 0,
      weighted   = 1,
      unknown    = -1
    };
  };

  struct ff_form {
    enum code {
      Gauss  = 2,
      dipole = 1
    };
  };

  struct absorption {
    enum code {
      truncatedExp = 3,
      exponential = 2,
      factorial   = 1
    };
  };

  struct ladder_weight {
    enum code {
      IntervalOnly,
      ReggeDiffusion,
      Regge
    };
  };

  struct ktform {
    enum code {
      IR0    = 2,
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
      props = 2,
      on    = 1,
      off   = 0,
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
  struct reconn_mode {
    enum code {
      run = 2,
      fix = 1,
      off = 0
    };
  };


  struct Run_Parameters {
    run_mode::code    runmode;
    weight_mode::code weightmode;
  public:
    Run_Parameters() {}
  };

  struct FormFactor_Parameters {
  public:
    ff_form::code form;
    size_t number;
    double norm, beta02, Lambda2, kappa, xi, bmax, accu;
    size_t bsteps;
    
    FormFactor_Parameters() {}
    FormFactor_Parameters(const FormFactor_Parameters & params) :
      form(params.form),number(-1),
      norm(params.norm),beta02(params.beta02), 
      Lambda2(params.Lambda2),kappa(params.kappa),xi(params.xi), 
      bmax(params.bmax), accu(params.accu), bsteps(params.bsteps) {}
  };

  struct Eikonal_Parameters {
  public:
    absorption::code absorp;
    double originalY, cutoffY, Ymax, lambda, Delta, beta02;
    double bmax, accu;

    Eikonal_Parameters() {}
    Eikonal_Parameters(const Eikonal_Parameters & params) :
      absorp(params.absorp),
      originalY(params.originalY), cutoffY(params.cutoffY), Ymax(params.Ymax),  
      lambda(params.lambda),Delta(params.Delta),beta02(params.beta02), 
      bmax(params.bmax), accu(params.accu) {}
  };

  struct Ladder_Parameters {
  public:
    double Q02, Qas2;
    Ladder_Parameters() {}
  };

  struct ShowerLink_Parameters {
  public:
    double CEKT2min;
    double KT2min;
    ShowerLink_Parameters() {}
  };
}

#endif
