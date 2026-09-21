#include "YFS/NLO/Real_Correction.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Spinor.H"
#include "METOOLS/Currents/Cancel_Probe.H"
#include "METOOLS/Main/Spin_Structure.H"

#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Point.H"
#include "EXTAMP/External_ME_Interface.H"
#include "MODEL/Main/Running_AlphaQED.H"

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <map>
#include <mutex>

using namespace YFS;
using namespace MODEL;
using namespace PHASIC;
using namespace ATOOLS;

/*!
  The runcard keys and labels for a given photon multiplicity, GENERATED
  rather than tabulated, so that an arbitrary number of photons needs no new
  row and no new class.

  The scheme reproduces the historical n=1 and n=2 names exactly:

      n=1   Real_Generator    Compare_Real   Real_Histogram        ROC
      n=2   RR_Generator      Compare_RR     RR_Histogram          RRC
      n=3   RRR_Generator     Compare_RRR    RRR_Histogram         RRRC
      n=4   RRRR_Generator    ...            ...                   RRRRC

  so existing runcards keep working and a higher multiplicity is configured the
  way the reader would guess.
*/
const Real_Config &Real_Correction::ConfigFor(size_t nphotons)
{
  static std::map<size_t, Real_Config> s_cfg;
  static std::mutex s_mtx;
  std::lock_guard<std::mutex> lock(s_mtx);

  if (nphotons == 0)
    THROW(fatal_error, "A real correction needs at least one photon.");

  std::map<size_t, Real_Config>::const_iterator it(s_cfg.find(nphotons));
  if (it != s_cfg.end()) return it->second;

  Real_Config c;
  if (nphotons == 1) {
    c.m_generator = "Real_Generator";
    c.m_generator2 = "Real_Generator2";
    c.m_compare   = "Compare_Real";
    c.m_writemom  = "Write_Real_Momenta";
    c.m_nmom      = "N_Real_Momenta";
    c.m_histdir   = "Real_Histogram";
    c.m_label     = "Real";
    c.m_tag       = "ROC";
    c.m_what      = "real";
  } else {
    const std::string pre(nphotons, 'R');      // RR, RRR, RRRR, ...
    std::string label, what("real");
    for (size_t i(0); i < nphotons; ++i) label += "Real";
    for (size_t i(1); i < nphotons; ++i) what += "-real";
    c.m_generator  = pre + "_Generator";
    c.m_generator2 = pre + "_Generator2";
    c.m_compare    = "Compare_" + pre;
    c.m_writemom   = "Write_" + pre + "_Momenta";
    c.m_nmom       = "N_" + pre + "_Momenta";
    c.m_histdir    = pre + "_Histogram";
    c.m_label      = label;
    c.m_tag        = pre + "C";
    c.m_what       = what;
  }
  return s_cfg.insert(std::make_pair(nphotons, c)).first->second;
}

/*!
  Build one real correction of the given photon multiplicity.

  This is the merge of the former Real:: and RealReal:: constructors, which
  differed only in the five settings keys, the histogram labels and the text
  of the error messages - all of which now come from the Real_Config row.
*/
Real_Correction::Real_Correction(const PHASIC::Process_Info &pi, size_t nphotons)
  : p_proc(nullptr), p_ampl(nullptr),
    m_sym(1.), m_factor(1.), m_check(false), m_writemom(false), m_nmom(100),
    m_fill(0), p_cmp(nullptr), m_nphotons(nphotons),
    p_cfg(&ConfigFor(nphotons)), m_dump(0)
{
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  const std::string gen
    (s[p_cfg->m_generator].SetDefault("Comix").Get<std::string>());
  m_gen = gen;
  // optional second EXTERNAL generator - if set, the comparison checks this
  // against the first directly (bypassing p_proc/Comix entirely) instead of
  // comparing the first against the internal ME.
  const std::string gen2
    (s[p_cfg->m_generator2].SetDefault("").Get<std::string>());
  m_check    = s[p_cfg->m_compare].SetDefault(0).Get<bool>();
  m_writemom = s[p_cfg->m_writemom].SetDefault(0).Get<bool>();
  m_nmom     = s[p_cfg->m_nmom].SetDefault(100).Get<int>();
  // Check_Real / Check_RR only exist for the two multiplicities that predate
  // this class; higher ones have no dump flag of their own, and the momenta
  // dumps they drive are a cross-generator debugging aid that has no external
  // generator to compare against up there anyway.
  m_dump     = (m_nphotons == 1 ? m_check_real
                                : (m_nphotons == 2 ? m_check_rr : 0));

  // Flavours come from the process info for EVERY multiplicity. RealReal used
  // to fill m_flavs from the external generator's argument list, inside the
  // branch that builds it - so with an internal generator its m_flavs stayed
  // empty and every m_flavs lookup in the diagnostics was out of bounds.
  for (const auto &f : pi.ExtractFlavours()) m_flavs.push_back(f);

  if (m_check && gen == "")
    THROW(fatal_error, "Need two generators to compare.");

  if (gen != "Comix" && gen != "Amegic") {
    PHASIC::External_ME_Args args(pi.m_ii.GetExternal(),
                                  pi.m_fi.GetExternal(),
                                  pi.m_maxcpl,
                                  gen);
    p_real_me.reset(PHASIC::Tree_ME2_Base::GetME2(args));
    if (!p_real_me) {
      msg_Error()<<pi;
      THROW(not_implemented, std::string("Couldn't find ")+p_cfg->m_what
                             +" ME for this process.");
    }
    MODEL::s_model->GetCouplings(m_cpls);
    p_real_me->SetCouplings(m_cpls);
    m_sym  = ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
    m_sym *= ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
    m_factor = 1./m_sym;
  }

  if (gen2 != "" && gen2 != "Comix" && gen2 != "Amegic") {
    if (!p_real_me)
      THROW(fatal_error, std::string(p_cfg->m_generator2)+" requires an "
                         "external "+p_cfg->m_generator+" too.");
    PHASIC::External_ME_Args args2(pi.m_ii.GetExternal(),
                                   pi.m_fi.GetExternal(),
                                   pi.m_maxcpl,
                                   gen2);
    p_real_me2.reset(PHASIC::Tree_ME2_Base::GetME2(args2));
    if (!p_real_me2)
      THROW(not_implemented, std::string("Couldn't find ")+p_cfg->m_what
                             +" ME for this process (generator 2).");
    p_real_me2->SetCouplings(m_cpls);
  }

  if (m_dump) {
    std::string filename(gen);
    for (const auto &f : m_flavs) { filename += "_"; filename += f.IDName(); }
    filename += "_"; filename += gen; filename += "_";
    if (FileExists(filename+"me.yaml"))        Remove(filename+"me.yaml");
    if (FileExists(filename+"ps-points.yaml")) Remove(filename+"ps-points.yaml");
    m_me_out.open((filename+"me.yaml").c_str(), std::ios_base::app);
    m_ps_out.open((filename+"ps-points.yaml").c_str(), std::ios_base::app);
    m_ps_out<<"MOMENTA:"<<std::endl;
  }

  p_cmp = std::make_unique<ME_Compare>(m_check, p_cfg->m_histdir, p_cfg->m_label);
}

Real_Correction::~Real_Correction() = default;

std::string Real_Correction::ActiveGenName() const {
  if (p_real_me) return m_gen;
  if (p_proc && p_proc->Generator()) return p_proc->Generator()->Name();
  return "none";
}

double Real_Correction::Calc_R(const ATOOLS::Vec4D_Vector& p)
  {
    // Zero-initialised: with no external generator configured (p_real_me
    // null) nothing below assigns it, yet the m_check comparisons read it.
    double external_real(0.);
    m_failcut = false;
    if(m_nlocuts && !p_proc->Trigger(p)) {
      m_failcut = true;
      msg_Debugging()<<"Rejecting real event, failed cuts"<<std::endl;
      return 0;
    }
    if(p_real_me) {
      external_real = Calc_External(p);
      if(p_real_me2){
        // compare two EXTERNAL generators directly, bypassing p_proc/Comix
        double external_real2 = p_real_me2->Calc(p)*m_factor
                                *ExternalFormFactor(p,m_flavs);
        if(m_check) p_cmp->CheckAgreement(p, external_real2, external_real,
                                           m_flavs, p_proc->NIn());
        return external_real;
      }
      if(!m_check) return external_real;
    }
    METOOLS::Cancel_Probe::Reset();
    p_ampl=CreateAmplitude(p);
    // rmode bits: 128=GeneratePoint(), 2=SetFixedScale(ampl scales), 1=disable
    // the selector before Trigger() so it can't reject the point - the
    // m_nlocuts check above already applies cuts explicitly, so Differential()
    // doesn't need to (and shouldn't) re-run the selector itself.
    int rmode = 128 +  2 + 1;
    Weights_Map iR = p_proc->Differential(*p_ampl, Variations_Mode::nominal_only,rmode);
    /*
      Comix's per-helicity amplitudes for this same point.

      FillAmplitudes() must follow Differential(), which is what actually
      evaluates the amplitude (see AddOns/EWSud/Comix_Interface.C:90-96).
      Spin_Amplitudes is helicity-indexed with SumSquare() summing |A|^2 over
      all combinations, so the ratio to the Differential is the spin average,
      flux and symmetry factors - a CONSTANT. Constant across phase space is
      the validation; drift means the two are not the same object.

      The point of this is CEEX: Ceex_Base hand-codes ~600 lines of spinor
      products and coupling assembly to build exactly these amplitudes, and
      every CEEX bug found so far lived in that layer. If Comix can supply
      them, that layer goes away.
    */
    static const bool campchk(getenv("SHERPA_COMIX_AMPS")!=NULL);
    if (m_keepamps) {
      std::vector<std::vector<Complex> > cols;
      m_spinamps.clear();
      p_proc->FillAmplitudes(m_spinamps, cols);
    }
    if (campchk && iR.Nominal()!=0.0) {
      std::vector<METOOLS::Spin_Amplitudes> amps;
      std::vector<std::vector<Complex> > cols;
      p_proc->FillAmplitudes(amps, cols);
      double ss(0.);
      size_t nhel(0);
      for (size_t a(0); a < amps.size(); ++a) { ss += amps[a].SumSquare(); nhel += amps[a].size(); }
      std::cerr<<"@@@ CAMP n="<<m_nphotons<<" nampl="<<amps.size()
               <<" nhel="<<nhel<<" sumsq="<<ss
               <<" diff="<<iR.Nominal()
               <<" ratio="<<(ss!=0.? iR.Nominal()/ss : 0.)<<std::endl;
    }
    // gauge-flip error estimator: re-evaluate the same point with the spinor
    // reference axis rotated. The result is gauge invariant analytically, so
    // the difference between the two evaluations is a direct probe of the
    // rounding error of THIS point - unlike the cancellation heuristic it
    // needs no calibration against a reference amplitude.
    double gaugedev(-1.0);
    if (m_check && iR.Nominal()!=0.0) {
      int sd(ATOOLS::Spinor<double>::DefaultGauge());
      ATOOLS::Spinor<double>::SetGauge(sd>0?sd-1:sd+1);
      Weights_Map iR2 = p_proc->Differential
        (*p_ampl, Variations_Mode::nominal_only,rmode);
      ATOOLS::Spinor<double>::ResetGauge();
      gaugedev=std::abs(iR2.Nominal()/iR.Nominal()-1.0);
    }
    if(iR.Nominal()==0) {
      if(p_ampl) p_ampl->Delete();
      if(m_check) msg_Out()<<"Real is 0"<<std::endl;
      return 0;
    }
    if(m_writemom && m_fill < m_nmom){
      m_ps_out<<std::setprecision(20)<<"  - ["<<std::endl;
      m_me_out<<std::setprecision(20)<<""<<m_fill<<":"<<std::endl;
      m_me_out<<std::setprecision(20)<<"  value: "<< (p_real_me ? external_real : iR.Nominal())<<std::endl;
      int j=0;
      for(auto k: p){
        m_ps_out<<"      [";
        if(m_flavs[j].IsAnti()) m_ps_out<<"-"<<m_flavs[j].Kfcode()<<", ";
        else m_ps_out<<m_flavs[j].Kfcode()<<", ";
        for(int i=0; i<4; i++){
          if(i!=3) m_ps_out<<k[i]<<",";
          else m_ps_out<<k[i];
        }
        m_ps_out<<"],"<<std::endl;
        j++;
      }
      m_ps_out<<"    ]"<<std::endl;
      m_fill++;
    } 
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
      const double cu(p_proc->Differential
                      (*au,Variations_Mode::nominal_only,rmode).Nominal());
      au->Delete();
      std::cerr<<"@@@ "<<p_cfg->m_tag<<"-ULP comix="<<std::abs(cu/iR.Nominal()-1.0)
               <<" ol="<<std::abs(eu/external_real-1.0)
               <<" refdev="<<std::abs(iR.Nominal()/external_real-1.0)
               <<std::endl;
    }
    // Rotation-invariance self-check: the ME is a Lorentz scalar, so a rigid
    // rotation must leave it unchanged. Any difference is the calculation's own
    // conditioning error on THIS point, measured against a symmetry the exact
    // answer obeys - no reference amplitude needed. Unlike the cancellation
    // heuristic this is not a proxy for the error, it IS an error estimate.
    double rotdev(-1.0);
    static const bool rotchk(getenv("SHERPA_ROT_CHECK")!=NULL);
    if (rotchk && iR.Nominal()!=0.0) {
      const double ca(cos(0.6)), sa(sin(0.6));
      Vec4D_Vector pr(p);
      for (size_t j(0);j<pr.size();++j)
        pr[j]=Vec4D(p[j][0],ca*p[j][1]+sa*p[j][3],p[j][2],
                    -sa*p[j][1]+ca*p[j][3]);
      Cluster_Amplitude *ar(CreateAmplitude(pr));
      const double r(p_proc->Differential
                     (*ar,Variations_Mode::nominal_only,rmode).Nominal());
      ar->Delete();
      rotdev=std::abs(r/iR.Nominal()-1.0);
    }
    if (m_check && external_real!=0.0) {
      std::cerr<<"@@@ "<<p_cfg->m_tag<<"-ROT2 rotdev="<<rotdev<<std::endl;
    }
    if (m_check && external_real!=0.0) {
      std::cerr<<"@@@ "<<p_cfg->m_tag<<" ratio="<<iR.Nominal()/external_real
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
                                       m_flavs, p_proc->NIn());
    return iR.Nominal();
  }

double Real_Correction::Calc_External(const ATOOLS::Vec4D_Vector &p){
  if(m_dump){
      m_ps_out<<std::setprecision(15)<<"  - ["<<std::endl;
      int j=0;
      for(auto k: p){
        m_ps_out<<"      [";
        if(m_flavs[j].IsAnti()) m_ps_out<<"-"<<m_flavs[j].Kfcode()<<", ";
        else m_ps_out<<m_flavs[j].Kfcode()<<", ";
        for(int i=0; i<4; i++){
          if(i!=3) m_ps_out<<k[i]<<",";
          else m_ps_out<<k[i];
        }
        m_ps_out<<"],"<<std::endl;
        j++;
      }
      m_ps_out<<"    ]"<<std::endl;
  }
   double R = p_real_me->Calc(p)*ExternalFormFactor(p,m_flavs);
  if(m_dump) {
    m_me_out<<std::setprecision(20)<<""<<m_fill<<":"<<std::endl;
    m_me_out<<std::setprecision(20)<<"  value: "<< (R)<<std::endl;
    m_fill++;
  }
  // if(m_writemom && m_fill < m_nmom) m_me_out<<std::setprecision(15)<<R/m_sym<<std::endl;
  return R*m_factor;
}

Cluster_Amplitude *Real_Correction::CreateAmplitude(const ATOOLS::Vec4D_Vector &p) const
{
  Cluster_Amplitude *ampl = Cluster_Amplitude::New();
  // NIn is the number of INCOMING legs, not the leg count. Comix crosses
  // in exactly the first NIn momenta (Single_Process.C: p[i] = i<NIn ?
  // -Leg(i)->Mom() : Leg(i)->Mom()), so passing p.size() here negates the
  // final state as well - which is equivalent to flipping the sign of every
  // mass relative to the momenta, and shows up wherever a mass term matters.
  ampl->SetNIn(p_proc->NIn());
  ampl->SetMS(p_proc->Generator());
  ampl->SetOrderQCD(p_proc->MaxOrder(0));
  ampl->SetMuF2(100);
  ampl->SetMuR2(100);
  ampl->SetMuQ2(100);
  ampl->SetMu2(100);
  for (size_t i(1);i<p_proc->MaxOrders().size();++i)
    ampl->SetOrderEW(ampl->OrderEW()+p_proc->MaxOrder(i));
  Int_Vector ci(p.size(), 0), cj(p.size(), 0);
  // Incoming legs are STORED negated - Comix undoes that when it crosses them
  // in. Storing them positive and setting NIn to the leg count (as this did)
  // flips every momentum instead: momentum is still conserved, so the result
  // stays finite and looks right away from any singularity, but p -> -p at
  // fixed m reverses each momentum relative to its mass. The error is then
  // invisible except where a mass term matters - i.e. inside the dead cone
  // theta <~ m/E, which is exactly where the real ME disagreed with OpenLoops.
  const size_t nin(p_proc->NIn());
  for (size_t i = 0; i < p.size(); ++i) {
    ampl->CreateLeg(i<nin?-p[i]:p[i], p_proc->Flavours()[i]);
  }
  ampl->SetProc(p_proc);
  return ampl;
}

void Real_Correction::SoftScan(const ATOOLS::Vec4D_Vector &p)
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
    const double cx(p_proc->Differential
                    (*a,Variations_Mode::nominal_only,128+2+1).Nominal());
    { double dmax(0.0);
      for (size_t j(0);j<pn.size();++j) {
        const Vec4D back(j<p_proc->NIn()?-a->Leg(j)->Mom():a->Leg(j)->Mom());
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
    const double cxr(p_proc->Differential
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
