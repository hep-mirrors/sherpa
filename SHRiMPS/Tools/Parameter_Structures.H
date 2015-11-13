#ifndef SHRIMPS_Tools_Parameter_Structures_H
#define SHRIMPS_Tools_Parameter_Structures_H

#include <ostream>

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
  // std::ostream & 
  // operator<<(std::ostream & s,const run_mode::code & runmode) {
  //   switch (runmode) {
  //   case run_mode::xsecs_only:
  //     s<<" XSecs only";
  //     break;
  //   case run_mode::elastic_events:
  //     s<<" Elastic events";
  //     break;
  //   case run_mode::single_diffractive_events:
  //     s<<" Single-diffractive events";
  //     break;
  //   case run_mode::quasi_elastic_events:
  //     s<<" Quasi-elastic events";
  //     break;
  //   case run_mode::inelastic_events:
  //     s<<" Inelastic events";
  //     break;
  //   case run_mode::all_min_bias:
  //     s<<" Minimum Bias";
  //     break;
  //   case run_mode::underlying_event:
  //     s<<" Underlying event";
  //     break;
  //   default:
  //     s<<" Unknown mode";
  //     break;
  //   }
  //   return s;
  // };

  struct weight_mode {
    enum code {
      unweighted = 0,
      weighted   = 1,
      unknown    = -1
    };
  };
  // std::ostream & operator<<(std::ostream & s, 
  // 			    const weight_mode::code & weightmode) {
  //   switch (weightmode) {
  //   case weight_mode::weighted:
  //     s<<" Weighted events";
  //     break;
  //   case weight_mode::unweighted:
  //     s<<" Unweighted events";
  //     break;
  //   default:
  //     s<<" Unknown mode";
  //     break;
  //   }
  //   return s;
  // };

  struct ff_form {
    enum code {
      Gauss  = 2,
      dipole = 1
    };
  };
  // std::ostream & operator<<(std::ostream & s, 
  // 			    const ff_form::code & ffform) {
  //   if (ffform==ff_form::Gauss)       s<<"Gaussian";
  //   else if (ffform==ff_form::dipole) s<<"Dipole";
  //   else                           s<<"unknown";
  //   return s;
  // };

  struct absorption {
    enum code {
      exponential = 2,
      factorial   = 1
    };
  };
  // std::ostream & operator<<(std::ostream & s, 
  // 			    const absorption::code & absorb) {
  //   if (absorb==absorption::exponential)    s<<"exponential";
  //   else if (absorb==absorption::factorial) s<<"factorial";
  //   else                               s<<"unknown";
  //   return s;
  // };

  struct ladder_weight {
    enum code {
      IntervalOnly,
      ReggeDiffusion,
      Regge
    };
  };
  // std::ostream & operator<<(std::ostream & s, 
  // 			    const ladder_weight::code & wt) {
  //   if (wt==ladder_weight::IntervalOnly)        s<<"y interval only ";
  //   else if (wt==ladder_weight::ReggeDiffusion) s<<"Regge factor + diffusion ";
  //   else if (wt==ladder_weight::Regge)          s<<"Regge factor ";
  //   else                                    s<<"unknown";
  //   return s;
  // }

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


  struct FormFactor_Parameters {
  public:
    ff_form::code form;
    int    number;
    double norm, beta0, Lambda2, kappa, xi, bmax, accu;

    FormFactor_Parameters() {}
    FormFactor_Parameters(const FormFactor_Parameters & params) :
      form(params.form),number(-1),
      norm(params.norm),beta0(params.beta0), 
      Lambda2(params.Lambda2),kappa(params.kappa),xi(params.xi), 
      bmax(params.bmax), accu(params.accu) {}
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
}

#endif
