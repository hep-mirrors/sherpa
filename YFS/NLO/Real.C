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
#include "PHASIC++/Main/Phase_Space_Point.H"

using namespace YFS;
using namespace MODEL;
using namespace PHASIC;

std::ofstream real_out, out_ps, out_mom;

Real::Real(const PHASIC::Process_Info& pi)  {
   /* Load Real ME */
   p_real_me = NULL;
   p_realproc = NULL;
   Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
   std::string gen = s["Real_Generator"].SetDefault("Comix").Get<std::string>();
   m_check = s["Compare_Real"].SetDefault(0).Get<bool>();
   m_writemom = s["Write_Real_Momenta"].SetDefault(0).Get<bool>();
   m_nmom = s["N_Real_Momenta"].SetDefault(100).Get<int>();
   for(auto f: pi.ExtractFlavours()) m_flavs.push_back(f);
   if(m_check && gen=="") THROW(fatal_error, "Need two generators to compare.");
   if(gen!="Comix" && gen != "Amegic"){
     PHASIC::External_ME_Args args(pi.m_ii.GetExternal(),
                                   pi.m_fi.GetExternal(),
                                   pi.m_maxcpl,
                                   gen);
     p_real_me =  PHASIC::Tree_ME2_Base::GetME2(args);
     if (!p_real_me)  THROW(not_implemented, "Couldn't find real ME for this process.");
     MODEL::s_model->GetCouplings(m_cpls);
     p_real_me->SetCouplings(m_cpls);
     m_sym =  ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
     m_sym *= ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
     m_factor = 1./m_sym;
    }
    if(m_check_real){
      std::string filename=gen;
      for(auto f: m_flavs) {
          filename+="_";
          filename+=f.IDName();
      }
      filename+="_";
      filename+=gen;
      filename+="_";
      if(FileExists(filename+"real.txt")) Remove(filename+"-real.txt");
      if(FileExists(filename+"ps-points.yaml")) Remove(filename+"-ps-points.yaml");
      if(FileExists(filename+"real.yaml")) Remove(filename+"real.yaml");
      real_out.open(filename+"real.yaml", std::ios_base::app); // append instead of overwrite
      out_ps.open(filename+"ps-points.yaml",std::ios_base::app);
      out_ps<<"MOMENTA:"<<std::endl;
  }
  // if(m_writemom){
  //   m_fill=0;
  //   std::string filename="Momenta";
  //   std::string MEfilename="ME";
  //   MEfilename+="_";
  //   MEfilename+=gen;
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
  //   out_mom.open(filename, std::ios_base::app);
  //   real_out.open(MEfilename, std::ios_base::app);
  //   out_mom<<"MOMENTA:"<<std::endl;
  //   real_out<<"ME:"<<std::endl;
  // }
  m_npoints = 0;
  m_nbad = 0;
  m_maxdev = 0.;
}

Real::~Real() {
  if(m_check && m_npoints>0){
    msg_Out()<<ATOOLS::om::bold<<ATOOLS::om::blue
             <<"###############################################"<<std::endl
             <<"Real ME comparison summary:"<<ATOOLS::om::reset<<std::endl
             <<"  points checked    = "<<m_npoints<<std::endl
             <<"  points mismatched = "<<m_nbad<<" ("
             <<std::setprecision(3)<<100.*m_nbad/m_npoints<<"%)"<<std::endl
             <<"  max |1-ratio|     = "<<std::setprecision(6)<<m_maxdev<<std::endl;
    if(m_nbad>0)
      msg_Out()<<ATOOLS::om::bold<<ATOOLS::om::red
               <<"WARNING: "<<m_nbad<<" / "<<m_npoints
               <<" points disagreed beyond tolerance!"<<ATOOLS::om::reset<<std::endl;
    else
      msg_Out()<<ATOOLS::om::bold<<ATOOLS::om::green
               <<"All points agreed within tolerance."<<ATOOLS::om::reset<<std::endl;
    msg_Out()<<ATOOLS::om::bold<<ATOOLS::om::blue
             <<"###############################################"
             <<ATOOLS::om::reset<<std::endl;
  }
}

double Real::Calc_R(const ATOOLS::Vec4D_Vector& p)
  {
    double external_real;
    m_failcut = false;
    if(m_nlocuts && !p_realproc->Trigger(p)) {
      m_failcut = true;
      msg_Debugging()<<"Rejecting real event, failed cuts"<<std::endl;
      return 0;
    }
    if(p_real_me) {
      if(!m_check) return Calc_External(p);
      external_real = Calc_External(p);
    }
    p_ampl=CreateAmplitude(p);
    // rmode bits: 128=GeneratePoint(), 2=SetFixedScale(ampl scales), 1=disable
    // the selector before Trigger() so it can't reject the point - the
    // m_nlocuts check above already applies cuts explicitly, so Differential()
    // doesn't need to (and shouldn't) re-run the selector itself.
    int rmode =  2 + 1;
    Weights_Map iR = p_realproc->Differential(*p_ampl, Variations_Mode::nominal_only,rmode);
    if(iR.Nominal()==0) {
      if(p_ampl) p_ampl->Delete();
      if(m_check) msg_Out()<<"Real is 0"<<std::endl;
      return 0;
    }
    if(m_writemom && m_fill < m_nmom){
      out_ps<<std::setprecision(20)<<"  - ["<<std::endl;
      real_out<<std::setprecision(20)<<""<<m_fill<<":"<<std::endl;
      real_out<<std::setprecision(20)<<"  value: "<< (p_real_me ? external_real : iR.Nominal())<<std::endl;
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
      m_fill++;
    } 
    // double ratio = iR.Nominal()/(m_factor*R);
    if(p_ampl) p_ampl->Delete();
    if(m_check){
      const double tol = 1e-6;
      double ratio = iR.Nominal()/external_real;
      double dev = std::abs(1.-ratio);
      m_npoints++;
      if(dev > m_maxdev) m_maxdev = dev;
      if(dev > tol){
        m_nbad++;
        double softe(-1.), softpt(-1.);
        for(size_t i(0);i<m_flavs.size() && i<p.size();++i){
          if(m_flavs[i].IsPhoton() && (softe<0. || p[i][0]<softe)){
            softe = p[i][0];
            softpt = p[i].PPerp();
          }
        }
        msg_Out()<<ATOOLS::om::bold<<ATOOLS::om::red<<"WARNING: "<<ATOOLS::om::reset
                 <<ATOOLS::om::red<<"Real ME mismatch: ratio = "
                 <<std::setprecision(10)<<ratio<<", |1-ratio| = "<<dev;
        if(softe>=0.) msg_Out()<<", softest photon: E = "<<softe<<", pT = "<<softpt;
        msg_Out()<<ATOOLS::om::reset<<std::endl;
      }
      else {
        msg_Debugging()<<ATOOLS::om::green<<"Real ME OK: ratio = "
                        <<std::setprecision(10)<<ratio<<ATOOLS::om::reset<<std::endl;
      }
    }
    return iR.Nominal();
  }

double Real::Calc_External(const ATOOLS::Vec4D_Vector &p){
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
  if(m_check_real) {
    real_out<<std::setprecision(20)<<""<<m_fill<<":"<<std::endl;
    real_out<<std::setprecision(20)<<"  value: "<< (R)<<std::endl;
    m_fill++;
  }
  // if(m_writemom && m_fill < m_nmom) real_out<<std::setprecision(15)<<R/m_sym<<std::endl;
  return R*m_factor;
}


Cluster_Amplitude *Real::CreateAmplitude(const ATOOLS::Vec4D_Vector &p) const
{
  Cluster_Amplitude *ampl = Cluster_Amplitude::New();
  ampl->SetNIn(p.size());
  ampl->SetMS(p_realproc->Generator());
  ampl->SetOrderQCD(p_realproc->MaxOrder(0));
  ampl->SetMuF2(100);
  ampl->SetMuR2(100);
  ampl->SetMuQ2(100);
  ampl->SetMu2(100);
  for (size_t i(1);i<p_realproc->MaxOrders().size();++i)
    ampl->SetOrderEW(ampl->OrderEW()+p_realproc->MaxOrder(i));
  Int_Vector ci(p.size(), 0), cj(p.size(), 0);
  for (size_t i = 0; i < p.size(); ++i) {
    ampl->CreateLeg(p[i], p_realproc->Flavours()[i]);
  }
  ampl->SetProc(p_realproc);
  return ampl;
}
