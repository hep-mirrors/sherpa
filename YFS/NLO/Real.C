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
   p_real_me =  PHASIC::Tree_ME2_Base::GetME2(pi);
   if (!p_real_me)  THROW(not_implemented, "Couldn't find real ME for this process.");
   MODEL::s_model->GetCouplings(m_cpls);
   PHASIC::External_ME_Args args(pi.m_ii.GetExternal(),
                                 pi.m_fi.GetExternal(),
                                 pi.m_maxcpl);
   p_real_me->SetCouplings(m_cpls);
   Flavour_Vector born_flavs;
   m_sym =  ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
   m_sym *= ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
   for(auto f: args.m_inflavs) m_flavs.push_back(f);
   for(auto f: args.m_outflavs) m_flavs.push_back(f);
   m_factor = m_rescale_alpha/m_sym;

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
        if(m_flavs[j].IsAnti()) out_ps<<"-"<<m_flavs[j].Kfcode()<<", ";
        else out_ps<<m_flavs[j].Kfcode()<<", ";
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
    return m_factor*R;
  }