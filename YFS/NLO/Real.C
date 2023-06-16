#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/Real.H"

#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "EXTAMP/External_ME_Interface.H"
#include "MODEL/Main/Running_AlphaQED.H"

using namespace YFS;
using namespace MODEL;

std::ofstream real_out, out_ps;


Real::Real(const PHASIC::Process_Info& pi)  {
   /* Load Real ME */
   PHASIC::Process_Info real_pi(pi);
   real_pi.m_mincpl[0] = pi.m_mincpl[0];
   real_pi.m_maxcpl[0] = pi.m_maxcpl[0];
   real_pi.m_mincpl[1] = pi.m_mincpl[1];
   real_pi.m_maxcpl[1] = pi.m_maxcpl[1];
   p_real_me =  PHASIC::Tree_ME2_Base::GetME2(real_pi);
   if (!p_real_me)  THROW(not_implemented, "Couldn't find real ME for this process.");
   MODEL::s_model->GetCouplings(m_cpls);
   /* Load color-correlated ME. TODO: orders */
   PHASIC::External_ME_Args args(real_pi.m_ii.GetExternal(),
                                 real_pi.m_fi.GetExternal(),
                                 real_pi.m_maxcpl);
   p_real_me->SetCouplings(m_cpls);
   Flavour_Vector born_flavs;
   for (int i = 0; i < args.m_outflavs.size()-1; ++i) born_flavs.push_back(args.m_outflavs[i]);
   m_sym  = ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
   m_sym *= ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
   double bornsym = ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
   bornsym*= ATOOLS::Flavour::FSSymmetryFactor(born_flavs);
   // m_sym/=bornsym;
   ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
   // m_factor = m_sym*p_real_me->AlphaQED()/2/M_PI;
  m_factor = 1./2/M_PI;///m_sym;
  if(m_check_real){
    if(FileExists("recola-real.txt")) Remove("recola-real.txt");
    if(FileExists("ps-points.yaml")) Remove("ps-points.yaml");
    real_out.open("recola-real.txt", std::ios_base::app); // append instead of overwrite
    out_ps.open("ps-points.yaml",std::ios_base::app);
    out_ps<<"MOMENTA:"<<std::endl;
  }

}

Real::~Real() {
}

double Real::Calc_R(const ATOOLS::Vec4D_Vector& p)
  {
    if(m_check_real){
      out_ps<<std::setprecision(15)<<"  - ["<<std::endl;
      int j=0;
      for(auto k: p){
        out_ps<<"      [";
        if(j==0) out_ps<<"11, ";
        if(j==1) out_ps<<"-11, ";
        if(j==2) out_ps<<"13, ";
        if(j==3) out_ps<<"-13, ";
        if(j==4) out_ps<<"22, ";
        for(int i=0; i<4; i++){
          if(i!=3) out_ps<<k[i]<<",";
          else out_ps<<k[i];
        }
        out_ps<<"],"<<std::endl;
        j++;
      }
      out_ps<<"    ]"<<std::endl;
  }
    
    double R = p_real_me->Calc(p);
    if(m_check_real) real_out<<std::setprecision(15)<<R/m_sym<<std::endl;
    return R*m_rescale_alpha*m_factor;///2/M_PI;
  }