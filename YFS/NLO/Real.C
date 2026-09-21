#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/Real.H"
#include <iostream>
#include <limits>
#include <cmath>
#include "METOOLS/Currents/Cancel_Probe.H"
#include "ATOOLS/Phys/Spinor.H"

#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "EXTAMP/External_ME_Interface.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Main/Phase_Space_Point.H"
#include "ATOOLS/YAML/yaml-cpp/yaml.h"

using namespace YFS;
using namespace MODEL;
using namespace PHASIC;

std::ofstream real_out, out_ps, out_mom;
std::ofstream wr_ps, wr_val;

Real::Real(const PHASIC::Process_Info& pi)  {
   /* Load Real ME */
   p_real_me = NULL;
   p_real_me2 = NULL;
   p_realproc = NULL;
   Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
   std::string gen = s["Real_Generator"].SetDefault("Comix").Get<std::string>();
   m_gen = gen;
   // optional second EXTERNAL generator - if set, Compare_Real checks this
   // against Real_Generator directly (bypassing p_realproc/Comix entirely)
   // instead of comparing Real_Generator against the internal ME.
   std::string gen2 = s["Real_Generator2"].SetDefault("").Get<std::string>();
   m_check = s["Compare_Real"].SetDefault(0).Get<bool>();
   m_writemom = s["Write_Real_Momenta"].SetDefault(0).Get<bool>();
   m_nmom = s["N_Real_Momenta"].SetDefault(100).Get<int>();
   m_fill = 0;
   const std::string readmom = s["Read_Real_Momenta"].SetDefault("").Get<std::string>();
   for(auto f: pi.ExtractFlavours()) m_flavs.push_back(f);
   if (!readmom.empty()) ReadMomenta(readmom);
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
    if(gen2!="" && gen2!="Comix" && gen2!="Amegic"){
      if(!p_real_me) THROW(fatal_error, "Real_Generator2 requires an external Real_Generator too.");
      PHASIC::External_ME_Args args2(pi.m_ii.GetExternal(),
                                     pi.m_fi.GetExternal(),
                                     pi.m_maxcpl,
                                     gen2);
      p_real_me2 = PHASIC::Tree_ME2_Base::GetME2(args2);
      if (!p_real_me2) THROW(not_implemented, "Couldn't find real ME for this process (generator 2).");
      p_real_me2->SetCouplings(m_cpls);
    }
    if(m_check_real){
      // pi.m_megenerator is the INTERNAL ME_Generator (Comix/Amegic) that
      // computes the amplitude being checked here; `gen` above is the
      // EXTERNAL Real_Generator it is compared against. Tag with both, so
      // runs made with different ME_Generator settings do not silently
      // overwrite each other's diagnostic output.
      const std::string megen(pi.m_megenerator.empty()?"auto":pi.m_megenerator);
      std::string filename=megen+"_vs_"+gen;
      for(auto f: m_flavs) {
          filename+="_";
          filename+=f.IDName();
      }
      filename+="_";
      if(FileExists(filename+"real.txt")) Remove(filename+"-real.txt");
      if(FileExists(filename+"ps-points.yaml")) Remove(filename+"-ps-points.yaml");
      if(FileExists(filename+"real.yaml")) Remove(filename+"real.yaml");
      real_out.open(filename+"real.yaml", std::ios_base::app); // append instead of overwrite
      out_ps.open(filename+"ps-points.yaml",std::ios_base::app);
      out_ps<<"MOMENTA:"<<std::endl;
  }
  // Write_Real_Momenta: dump up to N_Real_Momenta phase-space points (and
  // the value computed for each) so they can be replayed later with
  // Read_Real_Momenta -- for regression testing, or to isolate one specific
  // point (e.g. one SHERPA_SOFT_SCAN or a bad @@@ ROC line flagged) without
  // depending on the random sequence reproducing it again. Deliberately its
  // own files, independent of Compare_Real/m_check_real: the two features are
  // orthogonal, and m_check_real's out_ps/real_out are only opened when
  // Compare_Real is also on.
  if (m_writemom) {
    const std::string megen(pi.m_megenerator.empty()?"auto":pi.m_megenerator);
    std::string filename = "Real_events_" + megen;
    for (auto f: m_flavs) { filename += "_"; filename += f.IDName(); }
    if (FileExists(filename+"-momenta.yaml")) Remove(filename+"-momenta.yaml");
    if (FileExists(filename+"-values.yaml")) Remove(filename+"-values.yaml");
    wr_ps.open(filename+"-momenta.yaml", std::ios_base::app);
    wr_val.open(filename+"-values.yaml", std::ios_base::app);
    wr_ps<<"MOMENTA:"<<std::endl;
  }
  p_cmp = new ME_Compare(m_check, "Real_Histogram", "Real");
}

Real::~Real() {
  delete p_cmp;
}

double Real::Calc_R(const ATOOLS::Vec4D_Vector& p)
  {
    // Read_Real_Momenta: on the first call, replay every point read in by
    // ReadMomenta() instead of the one the phase-space generator produced,
    // then end the run. `replaying` guards against the recursive calls
    // ReplayMomenta makes back into Calc_R re-triggering this.
    static bool replaying(false);
    if (!m_replaypoints.empty() && !replaying) {
      replaying = true;
      ReplayMomenta();
      THROW(normal_exit,"Replay done.");
    }
    double external_real;
    m_failcut = false;
    if(m_nlocuts && !p_realproc->Trigger(p)) {
      m_failcut = true;
      msg_Debugging()<<"Rejecting real event, failed cuts"<<std::endl;
      return 0;
    }
    if(p_real_me) {
      external_real = Calc_External(p);
      if(p_real_me2){
        // compare two EXTERNAL generators directly, bypassing p_realproc/Comix
        double external_real2 = p_real_me2->Calc(p)*m_factor
                                *ExternalFormFactor(p,m_flavs);
        if(m_check) p_cmp->CheckAgreement(p, external_real2, external_real,
                                           m_flavs, p_realproc->NIn());
        return external_real;
      }
      if(!m_check) {
        // Same feature as below, reached here because an external
        // Real_Generator with Compare_Real off returns before ever building
        // the internal amplitude -- Write_Real_Momenta should still capture
        // the point and the value actually returned.
        if(m_writemom && m_fill < m_nmom) WriteMomentum(p, external_real);
        return external_real;
      }
    }
    METOOLS::Cancel_Probe::Reset();
    p_ampl=CreateAmplitude(p);
    // rmode bits: 128=GeneratePoint(), 2=SetFixedScale(ampl scales), 1=disable
    // the selector before Trigger() so it can't reject the point - the
    // m_nlocuts check above already applies cuts explicitly, so Differential()
    // doesn't need to (and shouldn't) re-run the selector itself.
    int rmode = 128 +  2 + 1;
    Weights_Map iR = p_realproc->Differential(*p_ampl, Variations_Mode::nominal_only,rmode);
    // gauge-flip error estimator: re-evaluate the same point with the spinor
    // reference axis rotated. The result is gauge invariant analytically, so
    // the difference between the two evaluations is a direct probe of the
    // rounding error of THIS point - unlike the cancellation heuristic it
    // needs no calibration against a reference amplitude.
    double gaugedev(-1.0);
    if (m_check && iR.Nominal()!=0.0) {
      int sd(ATOOLS::Spinor<double>::DefaultGauge());
      ATOOLS::Spinor<double>::SetGauge(sd>0?sd-1:sd+1);
      Weights_Map iR2 = p_realproc->Differential
        (*p_ampl, Variations_Mode::nominal_only,rmode);
      ATOOLS::Spinor<double>::ResetGauge();
      gaugedev=std::abs(iR2.Nominal()/iR.Nominal()-1.0);
    }
    if(iR.Nominal()==0) {
      if(p_ampl) p_ampl->Delete();
      if(m_check) msg_Out()<<"Real is 0"<<std::endl;
      return 0;
    }
    if(m_writemom && m_fill < m_nmom) WriteMomentum(p, p_real_me ? external_real : iR.Nominal());

    // double ratio = iR.Nominal()/(m_factor*R);
    if(p_ampl) p_ampl->Delete();
    // Input-conditioning probe: nudge every momentum component by one ulp and
    // re-evaluate BOTH generators. This measures the sensitivity of each to the
    // last representable digit of its input. If both move together the phase
    // space point itself is under-determined here; if only one moves, that one
    // is doing something the other is not.
    static const bool ulpchk(getenv("SHERPA_ULP_CHECK")!=NULL);
    if (ulpchk && iR.Nominal()!=0.0 && external_real!=0.0) {
      Vec4D_Vector pu(p);
      for (size_t j(0);j<pu.size();++j)
        for (int c(0);c<4;++c)
          pu[j][c]=std::nextafter(pu[j][c],
                                  std::numeric_limits<double>::infinity());
      const double eu(Calc_External(pu));
      Cluster_Amplitude *au(CreateAmplitude(pu));
      const double cu(p_realproc->Differential
                      (*au,Variations_Mode::nominal_only,rmode).Nominal());
      au->Delete();
      std::cerr<<"@@@ ULP comix="<<std::abs(cu/iR.Nominal()-1.0)
               <<" ol="<<std::abs(eu/external_real-1.0)
               <<" refdev="<<std::abs(iR.Nominal()/external_real-1.0)
               <<std::endl;
    }
    // Rotation-invariance self-check: 
    double rotdev(-1.0);
    static const bool rotchk(getenv("SHERPA_ROT_CHECK")!=NULL);
    if (rotchk && iR.Nominal()!=0.0) {
      const double ca(cos(0.6)), sa(sin(0.6));
      Vec4D_Vector pr(p);
      for (size_t j(0);j<pr.size();++j)
        pr[j]=Vec4D(p[j][0],ca*p[j][1]+sa*p[j][3],p[j][2],
                    -sa*p[j][1]+ca*p[j][3]);
      Cluster_Amplitude *ar(CreateAmplitude(pr));
      const double r(p_realproc->Differential
                     (*ar,Variations_Mode::nominal_only,rmode).Nominal());
      ar->Delete();
      rotdev=std::abs(r/iR.Nominal()-1.0);
    }
    if (m_check && external_real!=0.0) {
      std::cerr<<"@@@ ROT2 rotdev="<<rotdev<<std::endl;
    }
    if (m_check && external_real!=0.0) {
      std::cerr<<"@@@ ROC ratio="<<iR.Nominal()/external_real
               <<" dev="<<std::abs(iR.Nominal()/external_real-1.0)
               <<" cancel="<<METOOLS::Cancel_Probe::s_worst
               <<" gauge="<<gaugedev<<std::endl;
    }
    if (getenv("SHERPA_SOFT_SCAN")) {
      static bool done(false);
      // only bother with a point that actually disagrees, so the scan starts
      // inside the region under investigation
      if (!done && external_real!=0.0 &&
          std::abs(iR.Nominal()/external_real-1.0)>
          atof(getenv("SHERPA_SOFT_SCAN"))) {
        done=true;
        SoftScan(p);
        THROW(normal_exit,"SoftScan done.");
      }
    }
    if(m_check) p_cmp->CheckAgreement(p, iR.Nominal(), external_real,
                                       m_flavs, p_realproc->NIn());
    return iR.Nominal();
  }

void Real::SoftScan(const ATOOLS::Vec4D_Vector &p)
{
  using namespace ATOOLS;
  const size_t n(p.size());
  size_t ig(n-1);
  for (size_t i(0);i<n;++i) if (m_flavs[i].Kfcode()==22) ig=i;
  const Vec4D p1(p[0]), p2(p[1]), k0(p[ig]);
  const double me(m_flavs[0].Mass()), Eb(p1[0]), Eg(k0[0]);
  if (me<=0.0 || Eg<=0.0) { msg_Error()<<"SoftScan: bad point\n"; return; }

  // hard system: everything outgoing except the photon. Its configuration is
  // held fixed in its own rest frame, so only the photon direction varies.
  std::vector<Vec4D> hrest;
  std::vector<size_t> hidx;
  for (size_t i(2);i<n;++i) if (i!=ig) { hrest.push_back(p[i]); hidx.push_back(i); }
  const Vec4D Q0(p1+p2-k0);
  const double sQ0(Q0.Abs2());
  Poincare cms0(Q0);
  for (size_t j(0);j<hrest.size();++j) cms0.Boost(hrest[j]);

  const double mth(me/Eb);
  std::cerr<<"@@@ SCANHEAD Eg="<<Eg<<" me="<<me<<" Eb="<<Eb
           <<" m/E="<<mth<<" sqrtQ0="<<sqrt(sQ0)<<" nhard="<<hrest.size()
           <<std::endl;
  for (int it(0);it<25;++it) {
    const double th(mth*pow(10.0,-1.0+7.0*it/24.0));  // 0.1 .. 1e6 x m/E
    const Vec4D k(Eg,Eg*sin(th),0.0,Eg*cos(th));
    const Vec4D Q(p1+p2-k);
    const double sQ(Q.Abs2());
    if (sQ<=0.0) continue;
    const double sc(sqrt(sQ/sQ0));
    Vec4D_Vector pn(p);
    pn[ig]=k;
    Poincare cms(Q);
    for (size_t j(0);j<hrest.size();++j) {
      Vec4D q(hrest[j]*sc);
      cms.BoostBack(q);
      pn[hidx[j]]=q;
    }
    // eikonal of the two charged initial-state legs, mass terms included
    const double pk1(p1*k), pk2(p2*k), p12(p1*p2);
    const double S(2.0*p12/(pk1*pk2)-me*me/(pk1*pk1)-me*me/(pk2*pk2));
    const double ext(p_real_me?Calc_External(pn):0.0);
    Cluster_Amplitude *a(CreateAmplitude(pn));
    const double cx(p_realproc->Differential
                    (*a,Variations_Mode::nominal_only,128+2+1).Nominal());
    { double dmax(0.0);
      for (size_t j(0);j<pn.size();++j) {
        const Vec4D back(j<p_realproc->NIn()?-a->Leg(j)->Mom():a->Leg(j)->Mom());
        for (int c(0);c<4;++c)
          dmax=Max(dmax,std::abs(back[c]-pn[j][c]));
      }
      if (dmax>0.0)
        std::cerr<<"@@@ MOMCHANGE thr="<<th/mth<<" max|dp|="<<dmax<<std::endl; }
    a->Delete();
    // Same event, rigidly rotated off the beam axis. The physics is invariant,
    // but the light-cone decomposition p[0] +- p[3] is no longer a
    // cancellation for the beams, so any failure that is an artefact of the
    // axis-aligned numerics has to move. A failure that survives is structural.
    const double ca(cos(0.6)), sa(sin(0.6));
    Vec4D_Vector pr(pn);
    for (size_t j(0);j<pr.size();++j)
      pr[j]=Vec4D(pn[j][0],ca*pn[j][1]+sa*pn[j][3],pn[j][2],
                  -sa*pn[j][1]+ca*pn[j][3]);
    const double extr(p_real_me?Calc_External(pr):0.0);
    Cluster_Amplitude *ar(CreateAmplitude(pr));
    const double cxr(p_realproc->Differential
                     (*ar,Variations_Mode::nominal_only,128+2+1).Nominal());
    ar->Delete();
    std::cerr<<std::setprecision(12)<<"@@@ ROT thr="<<th/mth
             <<" int="<<cx<<" introt="<<cxr<<" extrot/ext="<<(ext?extr/ext:0.0)
             <<std::setprecision(6)<<std::endl;
    std::cerr<<std::setprecision(12)
             <<"@@@ SCAN thr="<<th/mth<<" th="<<th<<" S="<<S
             <<" int="<<cx<<" ext="<<ext
             <<" int/S="<<cx/S<<" ext/S="<<ext/S
             <<std::setprecision(6)<<std::endl;
  }
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
   double R = p_real_me->Calc(p)*ExternalFormFactor(p,m_flavs);
  if(m_check_real) {
    real_out<<std::setprecision(20)<<""<<m_fill<<":"<<std::endl;
    real_out<<std::setprecision(20)<<"  value: "<< (R)<<std::endl;
    m_fill++;
  }
  // if(m_writemom && m_fill < m_nmom) real_out<<std::setprecision(15)<<R/m_sym<<std::endl;
  return R*m_factor;
}


void Real::WriteMomentum(const ATOOLS::Vec4D_Vector &p, double value) {
  wr_ps<<std::setprecision(20)<<"  - ["<<std::endl;
  wr_val<<std::setprecision(20)<<""<<m_fill<<":"<<std::endl;
  wr_val<<std::setprecision(20)<<"  value: "<<value<<std::endl;
  int j=0;
  for(auto k: p){
    wr_ps<<"      [";
    if(m_flavs[j].IsAnti()) wr_ps<<"-"<<m_flavs[j].Kfcode()<<", ";
    else wr_ps<<m_flavs[j].Kfcode()<<", ";
    for(int i=0; i<4; i++){
      if(i!=3) wr_ps<<k[i]<<",";
      else wr_ps<<k[i];
    }
    wr_ps<<"],"<<std::endl;
    j++;
  }
  wr_ps<<"    ]"<<std::endl;
  m_fill++;
}

std::string Real::ActiveGenName() const {
  if (p_real_me) return m_gen;
  if (p_realproc && p_realproc->Generator()) return p_realproc->Generator()->Name();
  return "none";
}

void Real::ReadMomenta(const std::string &path) {
  SHERPA_YAML::Node root;
  try { root = SHERPA_YAML::LoadFile(path); }
  catch (const SHERPA_YAML::Exception &e) {
    THROW(fatal_error, "Read_Real_Momenta: cannot parse '"+path+"': "+e.what());
  }
  // Same MOMENTA: schema Write_Real_Momenta produces below, and that
  // Compare_Real's ps-points.yaml dump already used: a list of events, each
  // a list of [pdgid, E, px, py, pz] legs in m_flavs order (incoming legs
  // positive, matching this argument to Calc_R -- CreateAmplitude negates
  // them itself on the way in).
  const SHERPA_YAML::Node events = root["MOMENTA"];
  if (!events || !events.IsSequence())
    THROW(fatal_error, "Read_Real_Momenta: '"+path+"' has no MOMENTA: list.");
  for (size_t i(0);i<events.size();++i) {
    const SHERPA_YAML::Node ev = events[i];
    if (!ev.IsSequence() || ev.size()!=m_flavs.size()) {
      msg_Error()<<METHOD<<"(): "<<path<<": event "<<i<<" has "
                 <<(ev.IsSequence()?ev.size():0)<<" legs, expected "
                 <<m_flavs.size()<<" -- skipped.\n";
      continue;
    }
    Vec4D_Vector pt(m_flavs.size());
    for (size_t j(0);j<ev.size();++j) {
      const SHERPA_YAML::Node leg = ev[j];
      // [pdgid, x0, x1, x2, x3] with x0..x3 = (E,px,py,pz), Vec4D's own
      // component order -- WriteMomentum writes k[0..3] verbatim, so this is
      // a direct positional read, not a reordering.
      pt[j] = Vec4D(leg[1].as<double>(), leg[2].as<double>(),
                    leg[3].as<double>(), leg[4].as<double>());
    }
    m_replaypoints.push_back(pt);
  }
  msg_Info()<<METHOD<<"(): read "<<m_replaypoints.size()
            <<" phase-space point(s) from "<<path<<".\n";
}

void Real::ReplayMomenta() {
  msg_Out()<<"@@@ REPLAY start n="<<m_replaypoints.size()<<std::endl;
  for (size_t i(0);i<m_replaypoints.size();++i) {
    const double v(Calc_R(m_replaypoints[i]));
    msg_Out()<<"@@@ REPLAY i="<<i<<" value="<<v<<std::endl;
  }
  msg_Out()<<"@@@ REPLAY done"<<std::endl;
}

Cluster_Amplitude *Real::CreateAmplitude(const ATOOLS::Vec4D_Vector &p) const
{
  Cluster_Amplitude *ampl = Cluster_Amplitude::New();
  // NIn is the number of INCOMING legs, not the leg count. Comix crosses
  // in exactly the first NIn momenta (Single_Process.C: p[i] = i<NIn ?
  // -Leg(i)->Mom() : Leg(i)->Mom()), so passing p.size() here negates the
  // final state as well - which is equivalent to flipping the sign of every
  // mass relative to the momenta, and shows up wherever a mass term matters.
  ampl->SetNIn(p_realproc->NIn());
  ampl->SetMS(p_realproc->Generator());
  ampl->SetOrderQCD(p_realproc->MaxOrder(0));
  ampl->SetMuF2(100);
  ampl->SetMuR2(100);
  ampl->SetMuQ2(100);
  ampl->SetMu2(100);
  for (size_t i(1);i<p_realproc->MaxOrders().size();++i)
    ampl->SetOrderEW(ampl->OrderEW()+p_realproc->MaxOrder(i));
  Int_Vector ci(p.size(), 0), cj(p.size(), 0);
  // Incoming legs are STORED negated - Comix undoes that when it crosses them
  // in. Storing them positive and setting NIn to the leg count (as this did)
  // flips every momentum instead: momentum is still conserved, so the result
  // stays finite and looks right away from any singularity, but p -> -p at
  // fixed m reverses each momentum relative to its mass. The error is then
  // invisible except where a mass term matters - i.e. inside the dead cone
  // theta <~ m/E, which is exactly where the real ME disagreed with OpenLoops.
  const size_t nin(p_realproc->NIn());
  for (size_t i = 0; i < p.size(); ++i) {
    ampl->CreateLeg(i<nin?-p[i]:p[i], p_realproc->Flavours()[i]);
  }
  ampl->SetProc(p_realproc);
  return ampl;
}
