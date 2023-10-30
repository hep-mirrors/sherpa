#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/RealReal.H"

#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "EXTAMP/External_ME_Interface.H"
#include "MODEL/Main/Running_AlphaQED.H"

using namespace YFS;
using namespace MODEL;

std::ofstream rr_out, out_ps_rr;


RealReal::RealReal(const PHASIC::Process_Info& pi)  {
   /* Load RealReal ME */
   PHASIC::Process_Info real_pi(pi);
   real_pi.m_fi.m_nlotype = ATOOLS::nlo_type::loop;
   real_pi.m_mincpl[0] = pi.m_mincpl[0];
   real_pi.m_maxcpl[0] = pi.m_maxcpl[0];
   real_pi.m_mincpl[1] = pi.m_mincpl[1];
   real_pi.m_maxcpl[1] = pi.m_maxcpl[1];
   p_real_me = PHASIC::Tree_ME2_Base::GetME2(real_pi);;
   if (!p_real_me)  {
    msg_Error()<<real_pi;
    THROW(not_implemented, "Couldn't find real ME for this process.");
  }
   MODEL::s_model->GetCouplings(m_cpls);
   /* Load color-correlated ME. TODO: orders */
   PHASIC::External_ME_Args args(real_pi.m_ii.GetExternal(),
                                 real_pi.m_fi.GetExternal(),
                                 real_pi.m_maxcpl);
   p_real_me->SetCouplings(m_cpls);
   for(auto f: args.m_inflavs) m_flavs.push_back(f);
   for(auto f: args.m_outflavs) m_flavs.push_back(f);
   // p_real_me->SetSubType(sbt::qed);
   m_sym  = ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
   m_sym *= ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
   // m_factor = p_real_me->AlphaQED()/m_sym;
   // PRINT_VAR(m_sym);
   m_factor = m_rescale_alpha/m_sym;
   ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
   if(m_check_rr){
    if(FileExists("recola-real-real.txt")) Remove("recola-real-real.txt");
    if(FileExists("ps-points.yaml")) Remove("ps-points.yaml");
    rr_out.open("recola-real-real.txt", std::ios_base::app); // append instead of overwrite
    out_ps_rr.open("ps-points.yaml",std::ios_base::app);
    out_ps_rr<<"MOMENTA:"<<std::endl;
  }
} 

RealReal::~RealReal() {

}

double RealReal::Calc_R(const ATOOLS::Vec4D_Vector& p)
  {
    if(m_check_rr){
      out_ps_rr<<std::setprecision(15)<<"  - ["<<std::endl;
      int j=0;
      for(auto k: p){
        out_ps_rr<<"      [";
        if(m_flavs[j].IsAnti()) out_ps_rr<<"-"<<m_flavs[j].Kfcode()<<", ";
        else out_ps_rr<<m_flavs[j].Kfcode()<<", ";
        for(int i=0; i<4; i++){
          if(i!=3) out_ps_rr<<k[i]<<",";
          else out_ps_rr<<k[i];
        }
        out_ps_rr<<"],"<<std::endl;
        j++;
      }
      out_ps_rr<<"    ]"<<std::endl;
  }
    // double R = p_real_me->ME_Finite();
    double R = p_real_me->Calc(p);
    if(m_check_rr) rr_out<<std::setprecision(15)<<R/m_sym<<std::endl;
    return R*m_factor;
  }