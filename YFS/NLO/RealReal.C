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
using namespace PHASIC;

std::ofstream rr_out, out_ps_rr, out_mom_rr;


RealReal::RealReal(const PHASIC::Process_Info& pi)  {
  p_real_me = NULL;
  p_rrproc = NULL;
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  std::string gen = s["RR_Generator"].SetDefault("").Get<std::string>();
  m_check = s["Compare_RR"].SetDefault(0).Get<bool>();
  m_writemom = s["Write_RR_Momenta"].SetDefault(0).Get<bool>();
  m_nmom = s["N_RR_Momenta"].SetDefault(100).Get<int>();/* Load RealReal ME */
  PHASIC::Process_Info real_pi(pi);
  for(auto f: pi.ExtractFlavours()) m_flavs.push_back(f);
  if(m_check && gen=="") THROW(fatal_error, "Need two generators to compare.");
  if(gen!=""){
    PHASIC::External_ME_Args args(pi.m_ii.GetExternal(),
                               pi.m_fi.GetExternal(),
                               pi.m_maxcpl,
                               gen);
     p_real_me =  PHASIC::Tree_ME2_Base::GetME2(args);
     if (!p_real_me)  {
      msg_Error()<<real_pi;
      THROW(not_implemented, "Couldn't find Real-Real ME for this process.");
    }
     MODEL::s_model->GetCouplings(m_cpls);
     p_real_me->SetCouplings(m_cpls);
     for(auto f: args.m_inflavs) m_flavs.push_back(f);
     for(auto f: args.m_outflavs) m_flavs.push_back(f);
     // p_real_me->SetSubType(sbt::qed);
     m_sym  = ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
     m_sym *= ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
     double cplfac(1.0);
     cplfac *= pow(p_real_me->AlphaQED(),pi.m_maxcpl[1]);
     m_factor =  1./m_sym;
     // m_factor = 1./m_sym;
     if(m_check_rr){
      if(FileExists("recola-real-real.txt")) Remove("recola-real-real.txt");
      if(FileExists("ps-points.yaml")) Remove("ps-points.yaml");
      rr_out.open("recola-real-real.txt", std::ios_base::app); // append instead of overwrite
      out_ps_rr.open("ps-points.yaml",std::ios_base::app);
      out_ps_rr<<"MOMENTA:"<<std::endl;
    }
  }
  // if(m_writemom){
  //   m_fill=0;
  //   std::string filename="Momenta";
  //   std::string MEfilename="ME";
  //   for(auto f: m_flavs) {
  //     filename+="_";
  //     MEfilename+="_";
  //     filename+=f.IDName();
  //     MEfilename+=f.IDName();
  //   }
  //   filename+=".yaml";
  //   MEfilename+=".yaml";
  //   if(FileExists(filename)) Remove(filename);
  //   if(FileExists(MEfilename)) Remove(MEfilename);
  //   out_mom_rr.open(filename, std::ios_base::app);
  //   real_out.open(MEfilename, std::ios_base::app);
  //   out_mom_rr<<"MOMENTA:"<<std::endl;
  //   real_out<<"ME:"<<std::endl;
  // }
  // if(m_check){
  //   if (!ATOOLS::DirectoryExists("./Real_Histogram")) ATOOLS::MakeDir("./Real_Histogram");
  //   m_histograms1d["RRME_Dev"] = new Histogram(0,-1e-6, 1e-6, 100 );
  //   m_histograms1d["RRME_DevWide"] = new Histogram(0,-1, 1, 100 );
  // }
} 

RealReal::~RealReal() {

}

double RealReal::Calc_R(const ATOOLS::Vec4D_Vector& p){
  m_failcut = false;
  if(!p_rrproc->Trigger(p)) {
    m_failcut = true;
    return 0;
  }
  double external_real;
  if(p_real_me) {
      if(!m_check) return Calc_External(p);
      external_real = Calc_External(p);
  }
  p_ampl=CreateAmplitude(p);
  int rmode = 130;
  Weights_Map iR = p_rrproc->Differential(*p_ampl, Variations_Mode::nominal_only,rmode);
  if(p_ampl) p_ampl->Delete();
  if(m_check){
    double ratio = iR.Nominal()/external_real;
    // if(!IsEqual(ratio,1.,1e-4)){
    msg_Out()<<std::setprecision(15)<<"ratio = "<<ratio<<std::endl;
    msg_Out()<<std::setprecision(15)<<"external_real = "<<external_real<<std::endl;
       // m_histograms1d["RealME_Dev"]->Insert(1.-ratio);
       // m_histograms1d["RealME_DevWide"]->Insert(1.-ratio);
      // for (int i = 0; i < p.size(); ++i)
      // {
      //   msg_Out()<<"Flavour = "<<p_rrproc->Flavours()[i]<<std::endl
      //            <<"Momentum = "<<p[i]<<std::endl
      //            <<"P.PPerP() = "<<p[i].PPerp()<<std::endl
      //            <<"###############################################"<<std::endl;
      // }
    // }
  }
  return iR.Nominal();
}

double RealReal::Calc_External(const ATOOLS::Vec4D_Vector& p)
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
    if(m_check_rr) rr_out<<std::setprecision(15)<<R*m_factor<<std::endl;
    return R*m_factor;
  }


Cluster_Amplitude *RealReal::CreateAmplitude(const ATOOLS::Vec4D_Vector &p) const
{
  Cluster_Amplitude *ampl = Cluster_Amplitude::New();
  ampl->SetNIn(p.size());
  ampl->SetMS(p_rrproc->Generator());
  ampl->SetOrderQCD(p_rrproc->MaxOrder(0));
  ampl->SetMuF2(100);
  ampl->SetMuR2(100);
  ampl->SetMuQ2(100);
  ampl->SetMu2(100);
  for (size_t i(1);i<p_rrproc->MaxOrders().size();++i)
    ampl->SetOrderEW(ampl->OrderEW()+p_rrproc->MaxOrder(i));
  Int_Vector ci(p.size(), 0), cj(p.size(), 0);
  for (size_t i = 0; i < p.size(); ++i) {
    ampl->CreateLeg(p[i], p_rrproc->Flavours()[i]);
  }
  ampl->SetProc(p_rrproc);
  return ampl;
}
