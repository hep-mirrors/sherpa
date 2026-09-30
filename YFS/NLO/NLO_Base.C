#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "YFS/NLO/NLO_Base.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include <functional>
#include <map>
#include <cstdlib>
#include <iostream>
#include "YFS/NLO/Virtual.H"
#include "YFS/NLO/VirtualVirtual.H"
#include "YFS/NLO/Photon_Counterterm.H"
#include "MODEL/Main/Model_Base.H"
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include "YFS/NLO/NLO_Base_Internal.H"

using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;


// Lambda (Kaellen function) now lives once in YFS/Tools/Dipole.H.


NLO_Base::NLO_Base() {
  p_yfsFormFact = std::make_unique<YFS::YFS_Form_Factor>();
  p_nlodipoles = std::make_unique<YFS::Define_Dipoles>();
  // p_real/p_virt/p_realvirt/p_realreal/p_vv are default-null in the header.
  // p_realreal was the one the list here forgot, so m_rrtool read an
  // uninitialised pointer whenever SetProviders had not run yet.
  m_evts = 0;
  m_recola_evts = 0;
  m_realtool = 0;
  m_realvirt = 0;
  m_looptool = 0;
  m_rrtool = 0;
  m_vvtool = 0;
  m_zeroRV = 0;
  m_zeroRR = 0;
  m_nonZeroRR = 0;
  m_zeroV = 0;
  m_nonZeroRV=0;
  m_real_hard1 = 0.;
  m_rv_hard1 = 0.;
  m_rr_hard2 = 0.;
  m_real_hard2 = 0.;
  m_rv_hard2 = 0.;
  m_zero_real_amp = 0;
  m_ceex_done = false;
  m_softRV = 0;
  m_softRR = 0;
  m_rvUnstable = 0;
  m_rvHiC = 0;
  m_rvBlowup = 0;
  m_rvBlowupRtree0 = 0;
  m_rvBlowupHiC = 0;
  m_rvBlowupSoft = 0;
  m_rvBlowupHardWide = 0;
  BookHistograms();
  if (m_check_poles == 1) {
    if (!ATOOLS::DirectoryExists(m_debugDIR_NLO))
      ATOOLS::MakeDir(m_debugDIR_NLO);
    m_histograms1d["SinglePoleCD"] = std::make_unique<Histogram>(0, 0, 25, 25);
    m_histograms1d["SinglePoleVV"] = std::make_unique<Histogram>(0, 0, 25, 25);
    m_histograms1d["DoublePoleVV"] = std::make_unique<Histogram>(0, 0, 25, 25);
    m_histograms1d["OneLoopEpsLP"] = std::make_unique<Histogram>(0, -1.5, -0.5, 50);
    m_histograms1d["OneLoopEpsYFS"] = std::make_unique<Histogram>(0, -1.5, -0.5, 50);
    m_histograms1d["RealLoopEpsLP"] = std::make_unique<Histogram>(0, -5, 0.0, 50);
    m_histograms1d["RealLoopEpsYFS"] = std::make_unique<Histogram>(0, -5, 0.0, 50);
    m_histograms1d["relativediff"] = std::make_unique<Histogram>(0, -20., -5.0, 50);
    m_histograms1d["RVSinglePoleCD"] = std::make_unique<Histogram>(0, 0, 25, 25);
    m_histograms2d["REAL_SUB"] =
        std::make_unique<Histogram_2D>(0, 0, sqrt(m_s) / 2., 200, 0, 2 * M_PI, 20);
    m_histograms2d["REAL"] =
        std::make_unique<Histogram_2D>(0, 0, sqrt(m_s) / 2., 200, 0, 2 * M_PI, 20);
  }
  if (m_rv_cancel_hist) {
    if (!ATOOLS::DirectoryExists(m_debugDIR_NLO))
      ATOOLS::MakeDir(m_debugDIR_NLO);
    m_histograms1d["RV_tot_by_logC_w"] = std::make_unique<Histogram>(0, -16., 0., 80);
    m_histograms1d["RV_tot_by_logC_n"] = std::make_unique<Histogram>(0, -16., 0., 80);
    m_histograms1d["RV_tot_by_Efrac_w"] = std::make_unique<Histogram>(0, 0., 0.5, 100);
    m_histograms1d["RV_tot_by_Efrac_n"] = std::make_unique<Histogram>(0, 0., 0.5, 100);
    m_histograms1d["RV_MEstab_all"] = std::make_unique<Histogram>(0, -2., 40., 84);
    m_histograms1d["RV_MEstab_hardwide"] = std::make_unique<Histogram>(0, -2., 40., 84);
    m_histograms1d["RV_tot_by_MEstab_w"] = std::make_unique<Histogram>(0, -2., 40., 84);
  }
}

NLO_Base::~NLO_Base() {
  WriteHistograms();
  msg_Out()<<"Total zero V: "<<m_zeroV<<std::endl;
  msg_Out()<<"Total zero RV: "<<m_zeroRV<<std::endl;
  msg_Out()<<"Total zero RR: "<<m_zeroRR<<std::endl;
  msg_Out()<<"Total non-zero RR: "<<m_nonZeroRR<<std::endl;
  msg_Out()<<"Total non-zero RV: "<<m_nonZeroRV<<std::endl;
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    int gbuf[3] = {m_softRV, m_rvUnstable, m_softRR};
    mpi->Allreduce(gbuf, 3, MPI_INT, MPI_SUM);
    m_softRV = gbuf[0];
    m_rvUnstable = gbuf[1];
    m_softRR = gbuf[2];
  }
#endif
  msg_Out()<<"Total soft RV skipped: "<<m_softRV<<std::endl;
  msg_Out()<<"Total unstable-ME RV skipped (RV_ME_MAX_RATIO): "<<m_rvUnstable<<std::endl;
  if (RVMode() == rvmode::remainder && m_realvirt) {
#ifdef USING__MPI
    if (mpi->Size() > 1) {
      int rb[2] = {m_rvPoleFail, m_rvNoVirt};
      mpi->Allreduce(rb, 2, MPI_INT, MPI_SUM);
      m_rvPoleFail = rb[0];
      m_rvNoVirt = rb[1];
    }
#endif
    msg_Out()<<"RV_MODE 1: pole mismatches > 1e-6: "<<m_rvPoleFail
             <<", photons without a Born virtual: "<<m_rvNoVirt<<std::endl;
    msg_Out()<<"RV_MODE 1 (this rank): loop frame "<<RVLoopFrame()
             <<", two-frame checks "<<m_rvLoopChecked<<", dropped as unstable "
             <<m_rvLoopUnstable<<std::endl;
  }
  msg_Out()<<"Total soft RR pairs skipped: "<<m_softRR<<std::endl;
  if (m_rv_cancel_hist) {
#ifdef USING__MPI
    if (mpi->Size() > 1) {
      int buf[6] = {m_rvHiC,      m_rvBlowup,    m_rvBlowupRtree0,
                    m_rvBlowupHiC, m_rvBlowupSoft, m_rvBlowupHardWide};
      mpi->Allreduce(buf, 6, MPI_INT, MPI_SUM);
      m_rvHiC = buf[0];
      m_rvBlowup = buf[1];
      m_rvBlowupRtree0 = buf[2];
      m_rvBlowupHiC = buf[3];
      m_rvBlowupSoft = buf[4];
      m_rvBlowupHardWide = buf[5];
    }
#endif
    msg_Out()<<"RV photons with C>=1 (subtraction not cancelling): "<<m_rvHiC<<std::endl;
    msg_Out()<<"RV blow-ups |tot|>1e3*|Born|: "<<m_rvBlowup
             <<"  (of these: rtree==0: "<<m_rvBlowupRtree0
             <<", C>=1: "<<m_rvBlowupHiC
             <<", soft E/sqrt(s)<0.01: "<<m_rvBlowupSoft
             <<", HARD WIDE-ANGLE: "<<m_rvBlowupHardWide<<")"<<std::endl;
    if (m_rvBlowup>0)
      msg_Out()<<"  -> hard wide-angle fraction of blow-ups: "
               <<(100.*m_rvBlowupHardWide/m_rvBlowup)<<"%"
               <<" (if >0, instability is NOT confined to soft/collinear)"<<std::endl;
    if (m_rvBlowup>0)
      msg_Out()<<"  -> rtree==0 fraction of blow-ups: "
               <<(100.*m_rvBlowupRtree0/m_rvBlowup)<<"%"
               <<" (mechanism confirmed if ~100%)"<<std::endl;
  }
  msg_Out()<<"Total zero real amplitudes: "<<m_zero_real_amp<<std::endl;
  msg_Out()<<"Total events : "<<m_evts<<std::endl;
  ReportLoopHelicity();
}

void NLO_Base::SetProviders(YFS::Virtual *virt, YFS::Real *real,
                            YFS::RealVirtual *realvirt, YFS::RealReal *realreal,
                            YFS::VirtualVirtual *vv) {
  p_virt     = virt;
  p_real     = real;
  p_realvirt = realvirt;
  p_realreal = realreal;
  p_vv       = vv;
  m_looptool = (p_virt     != NULL);
  m_realtool = (p_real     != NULL);
  m_realvirt = (p_realvirt != NULL);
  m_rrtool   = (p_realreal != NULL);
  m_vvtool   = (p_vv       != NULL);
}

void NLO_Base::Init(Flavour_Vector &flavs, Vec4D_Vector &plab,
                    Vec4D_Vector &born) {
  m_rawbeta.clear();
  m_flavs = flavs;
  m_plab = plab;
  m_bornMomenta = born;
}

// ======================================================================
//  General-n real corrections: one subset recursion for every multiplicity
// ======================================================================

YFS::Real_Correction *NLO_Base::RealProvider(size_t nphotons) const
{
  if (nphotons == 1) return p_real;
  if (nphotons == 2) return p_realreal;
  if (nphotons < m_realprov.size()) return m_realprov[nphotons];
  return NULL;
}

void NLO_Base::SetRealProvider(size_t nphotons, YFS::Real_Correction *prov)
{
  if (m_realprov.size() <= nphotons) m_realprov.resize(nphotons+1, NULL);
  m_realprov[nphotons] = prov;
}

size_t NLO_Base::MaxRealPhotons() const
{
  size_t n(0);
  if (p_real)     n = 1;
  if (p_realreal) n = 2;
  for (size_t i(m_realprov.size()); i-- > 3; )
    if (m_realprov[i] != NULL) { n = Max(n, i); break; }
  return n;
}

size_t NLO_Base::RequestedMaxRealPhotons()
{
  static const int n
    (ATOOLS::Settings::GetMainSettings()["YFS"]["NLO_MAX_PHOTONS"]
     .SetDefault(2).Get<int>());
  return n < 1 ? 1 : (size_t)n;
}
