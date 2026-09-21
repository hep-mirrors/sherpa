/*!
  \file Ceex_Partitions.C

  The ISR/FSR partition sum. Calculate() walks all 2^n assignments of the
  photons to the initial or final line and accumulates the betas coherently;
  that coherence IS initial-final interference in CEEX.

*/

#include "YFS/CEEX/Ceex_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "EXTAMP/External_ME_Interface.H"
#include "PHASIC++/Process/External_ME_Args.H"

using namespace YFS;


void Ceex_Base::CalculateSfactors() {
  m_Sfac_ini.clear();
  m_Sfac_fin.clear();
  const Complex qratio(m_qe != 0. ? -m_qf/m_qe : 0., 0.);
  for (size_t i(0); i < m_allphotons.size(); ++i) {
    m_Sfac_ini.push_back(Sfactor(m_bornmomenta[0], m_bornmomenta[1],
                                 m_allphotons[i], m_PhoHel[i]));
    m_Sfac_fin.push_back(qratio * Sfactor(m_pceex[2], m_pceex[3],
                                          m_allphotons[i], m_PhoHel[i]));
  }
}


void Ceex_Base::PartitionStart(int &last) {
  const bool fsr(HasFSR());
  m_isrflag.assign(m_allphotons.size(), fsr ? 0 : 1);
  last = fsr ? 0 : 1;
}


void Ceex_Base::PartitionPlus(int &last) {
  const size_t n(m_isrflag.size());
  if (n == 0) { last = 2; return; }
  if (n == 1) last = 1;
  ++m_isrflag[0];
  for (size_t i(0); i < n; ++i)
    if (m_isrflag[i] == 2) {
      m_isrflag[i] = 0;
      if (i + 1 < n) {
        ++m_isrflag[i+1];
        if (m_isrflag[n-1] == 2) last = 2;
      } else {
        last = 2;   // carried off the top photon: every partition is done
      }
    }
}


void Ceex_Base::MakePhotonHel() {
  // CEEX_PIN_PHOTON_HEL pins the helicity for validation. The random number is
  // still drawn, so the RNG stream - and therefore the phase-space sequence - is
  // identical to an unpinned run. That is what lets two runs be summed point by
  // point to recover the helicity sum.
  static const char *pin(getenv("CEEX_PIN_PHOTON_HEL"));
  static const int pinned(pin ? atoi(pin) : 0);
  m_PhoHel.clear();
  for (size_t i(0); i < m_allphotons.size(); ++i) {
    const int h(ran->Get() < 0.5 ? 1 : -1);
    m_PhoHel.push_back(pinned ? pinned : h);
  }
}


void Ceex_Base::Calculate() {
  // Lab copies first: everything below boosts into the CEEX frame.
  m_allphotons_lab = m_isrphotons;
  m_allphotons_lab.insert(m_allphotons_lab.end(),
                          m_fsrphotons.begin(), m_fsrphotons.end());
  for (size_t i(0); i < m_isrphotons.size(); ++i) m_cms.Boost(m_isrphotons[i]);
  for (size_t i(0); i < m_fsrphotons.size(); ++i) m_cms.Boost(m_fsrphotons[i]);
  for (size_t i(0); i < m_bornmomenta.size(); ++i) m_cms.Boost(m_bornmomenta[i]);

  m_justdumped = false;
  m_rhocrud = 0.0;
  m_beta10 = 0.0;
  m_beta01 = 0.0;
  m_beta00 = 0.0;
  BuildCeexMomenta();
  ZerAmplit();
  m_allphotons = m_isrphotons;
  m_allphotons.insert(m_allphotons.end(),
                      m_fsrphotons.begin(), m_fsrphotons.end());
  MakePhotonHel();
  CalculateSfactors();
  m_spincache.resize(1 + 4*m_allphotons.size());
  m_spinvalid.assign(m_spincache.size(), 0);
  m_realphot.assign(m_allphotons.size(), Amplitude());
  m_svarQ = (m_pceex[2] + m_pceex[3]).Abs2();
  m_sQ    = m_svarQ;
  MakePropT(m_pceex);

  
  int last(0), nparts(0);
  
  static const size_t maxphot(Settings::GetMainSettings()["CEEX"]["MAX_PARTITION_PHOTONS"]
                              .SetDefault(100).Get<int>());
  if (m_allphotons.size() > maxphot) {
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Error()<<METHOD<<"(): "<<m_allphotons.size()<<" photons would need "
                 <<"2^"<<m_allphotons.size()<<" partitions. Falling back to "
                 <<"the all-ISR partition for events above CEEX: "
                 <<"MAX_PARTITION_PHOTONS = "<<maxphot<<", which DROPS the "
                 <<"initial-final interference on them. Reported once."
                 <<std::endl;
    }
    m_isrflag.assign(m_allphotons.size(), 1);
    last = 1;
  } else {
    PartitionStart(last);
  }
  if (m_allphotons.size() > m_maxnphot) m_maxnphot = m_allphotons.size();
  for (;;) {
    ++nparts;
    Vec4D PX(m_pceex[0] + m_pceex[1]);
    Complex sProd(1., 0.);
    std::vector<Complex> Sactu(m_allphotons.size(), Complex(1., 0.));
    for (size_t j(0); j < m_allphotons.size(); ++j) {
      Sactu[j] = m_isrflag[j] ? m_Sfac_ini[j] : m_Sfac_fin[j];
      sProd   *= Sactu[j];
      if (m_isrflag[j]) PX -= m_allphotons[j];
    }
    const double svarX(PX.Abs2());
    if (m_checkxs) {
      const Complex pz(1./Complex(svarX - sqr(m_MZ), m_gZ*svarX/m_MZ));
      const double a(std::abs(pz));
      if (nparts == 1) { m_pzmin = m_pzmax = a; }
      else { m_pzmin = Min(m_pzmin, a); m_pzmax = Max(m_pzmax, a); }
    }
    if (svarX <= 0. || m_svarQ <= 0.) {
      if (last) break;
      PartitionPlus(last);
      if (last == 2) break;
      continue;
    }
    m_sp   = svarX;
    m_Sprod = sProd;
    m_cfac = sProd * Complex(svarX/m_svarQ, 0.);
    MakeProp();
    // Electroweak form factors at THIS partition's scale, and the scattering
    // angle the WW/ZZ boxes depend on. No-op unless CEEX: WEAK is set.
    MakeEWFF(m_sp, cos(m_pceex[2].Theta()));
    MakeBoxMandelstams(PX);

    InfraredSubtractedME_0_0();
    InfraredSubtractedME_0_1();
    for (size_t j(0); j < m_allphotons.size(); ++j) {
      if (m_checkxs && m_isrflag[j]) {
        CheckSoftFactor(m_allphotons[j]);
        CheckUVDiagonality(m_allphotons[j]);
      }
      if (m_isrflag[j]) {
        InfraredSubtractedME_1_0(m_allphotons[j], m_PhoHel[j], sProd, Sactu[j],
                                 1 + 4*(int)j, (int)j);
      } else {
        // CKine = (q1+q2+k)^2/(q1+q2)^2, KKMC's svarX1/svarQ.
        const double CKine((m_pceex[2] + m_pceex[3]
                            + m_allphotons[j]).Abs2() / m_svarQ);
        InfraredSubtractedME_1_0_FSR(m_allphotons[j], m_PhoHel[j],
                                     sProd, Sactu[j], CKine, 3 + 4*(int)j, (int)j);
      }
    }

    if (last == 1) break;
    PartitionPlus(last);
    if (last == 2) break;
  }
  m_nparts = nparts;
  if (m_checkxs && nparts > 1 && m_pzmin > 0.) {
    const size_t n(Min(m_allphotons.size(), size_t(15)));
    m_pzspread[n] = Max(m_pzspread[n], m_pzmax/m_pzmin);
    m_pzsum[n] += m_pzmax/m_pzmin;
    ++m_pzn[n];
    const double ecm((m_momenta[0]+m_momenta[1]).Mass());
    for (int c(0); c < 4; ++c) {
      const double xcut(m_pzxcut[c]);
      Vec4D Xhard(m_pceex[0] + m_pceex[1]), Xall(Xhard);
      size_t nsoft(0);
      for (size_t j(0); j < m_allphotons.size(); ++j) {
        const double x(2.*m_allphotons[j][0]/ecm);
        if (x > xcut) { Xhard -= m_allphotons[j]; Xall -= m_allphotons[j]; }
        else { Xall -= m_allphotons[j]; ++nsoft; }
      }
      if (nsoft == 0) continue;
      const double sh(Xhard.Abs2()), sa(Xall.Abs2());
      if (sh <= 0. || sa <= 0.) continue;
      const double a(std::abs(1./Complex(sh - sqr(m_MZ), m_gZ*sh/m_MZ)));
      const double b(std::abs(1./Complex(sa - sqr(m_MZ), m_gZ*sa/m_MZ)));
      const double r(a > b ? a/b : b/a);
      m_pzsoft[c] = Max(m_pzsoft[c], r);
      m_pzsoftsum[c] += r; ++m_pzsoftn[c];
      m_pzsoftsav[c] += double(1u << Min(nsoft, size_t(20)));
    }
  }
 
  {
    const size_t n(m_allphotons.size());
    const int want(!HasFSR() || n > maxphot ? 1
                   : (n < 31 ? (1 << n) : nparts));
    if (nparts != want) {
      ++m_partbad;
      static bool warned(false);
      if (!warned) {
        warned = true;
        msg_Error()<<METHOD<<"(): partition counter enumerated "<<nparts
                   <<" partitions for "<<n<<" photons, expected "<<want
                   <<". Reported once."<<std::endl;
      }
    }
    ++m_partn;
  }

  MakeRho();

  static const double dumpxmin(getenv("SHERPA_CEEX_DUMP_XMIN") ?
                               atof(getenv("SHERPA_CEEX_DUMP_XMIN")) : 0.);
  static const size_t dumpn(getenv("SHERPA_CEEX_DUMP_NPHOT") ?
                            atoi(getenv("SHERPA_CEEX_DUMP_NPHOT")) : 1);
  const double xgam(!m_allphotons.empty() && m_momenta.size() >= 2 ?
                    2.*m_allphotons[0][0]/(m_momenta[0]+m_momenta[1]).Mass() : 0.);
  if (m_checkxs && !m_ceexdumped && m_allphotons.size() == dumpn &&
      xgam > dumpxmin && m_momenta.size() >= 6 && m_result0 != 0.) {
    
    Vec4D bal(m_momenta[0] + m_momenta[1] - m_momenta[4] - m_momenta[5]);
    for (size_t i(0); i < m_allphotons.size(); ++i) bal -= m_allphotons[i];
    const double worst(Max(Max(dabs(bal[0]), dabs(bal[1])),
                           Max(dabs(bal[2]), dabs(bal[3]))));
    if (worst > 1e-6) {
      msg_Error()<<METHOD<<"(): CEEX dump skipped, momentum imbalance "
                 <<worst<<std::endl;
      return;
    }
    m_ceexdumped = true;
    m_justdumped = true;
    std::ofstream pt("ceex_point.dat");
    pt << std::setprecision(17);
    const int idx[5] = {0, 1, 4, 5, -1};
    for (int n(0); n < 4; ++n)
      pt << m_momenta[idx[n]][0] << " " << m_momenta[idx[n]][1] << " "
         << m_momenta[idx[n]][2] << " " << m_momenta[idx[n]][3] << "\n";
    
    pt << m_allphotons[0][0] << " " << m_allphotons[0][1] << " "
       << m_allphotons[0][2] << " " << m_allphotons[0][3] << "\n";
    pt << "NPHOT " << m_allphotons.size() << "\n";
    for (size_t i(0); i < m_allphotons.size(); ++i)
      pt << m_allphotons[i][0] << " " << m_allphotons[i][1] << " "
         << m_allphotons[i][2] << " " << m_allphotons[i][3] << " "
         << m_PhoHel[i] << "\n";
    
    {
      Amplitude bref;
      BornAmplitude(m_pceex, bref);
      msg_Error() << std::setprecision(10);
      msg_Error() << "@@@ SHMASS nflav=" << m_flavs.size()
                << " m1=" << (m_flavs.size()>0? m_flavs[0].Mass():-1.)
                << " m2=" << (m_flavs.size()>1? m_flavs[1].Mass():-1.)
                << " m3=" << (m_flavs.size()>2? m_flavs[2].Mass():-1.)
                << " m4=" << (m_flavs.size()>3? m_flavs[3].Mass():-1.)
                << " mass_I=" << m_mass_I << " mass_F=" << m_mass_F
                << std::endl;
      
      msg_Error() << std::setprecision(16) << "@@@ SHPMASS";
      for (size_t i = 0; i < m_pceex.size(); ++i)
        msg_Error() << " p" << i << "=" << m_pceex[i].Mass();
      msg_Error() << std::setprecision(10) << std::endl;
      msg_Error() << "@@@ SHCPL sin2tw=" << m_sin2tw << " norm=" << m_norm
                << " ve=" << m_ve << " vf=" << m_vf
                << " ae=" << m_ae << " af=" << m_af
                << " qe=" << m_qe << " qf=" << m_qf
                << " CZ1=" << CouplingZ(1., 1) << " CZ0=" << CouplingZ(1., 0)
                << " CG=" << CouplingG() << std::endl;
      
      msg_Error() << std::setprecision(16)
                << "@@@ SHVIRT hasfsr=" << (HasFSR() ? 1 : 0)
                << " s=" << m_Sc.real() << " t=" << m_Tc.real()
                << " u=" << m_Uc.real() << " sQ=" << m_sQ
                << " deltI=(" << m_vertexI.real() << "," << m_vertexI.imag() << ")"
                << " deltF=(" << m_vertexF.real() << "," << m_vertexF.imag() << ")"
                << " BoxGGtu=(" << m_BoxGGtu.real() << "," << m_BoxGGtu.imag() << ")"
                << " BoxGZtu=(" << m_BoxGZtu.real() << "," << m_BoxGZtu.imag() << ")"
                << " BoxGGut=(" << m_BoxGGut.real() << "," << m_BoxGGut.imag() << ")"
                << " BoxGZut=(" << m_BoxGZut.real() << "," << m_BoxGZut.imag() << ")"
                << " rawGG=(" << m_rawboxGG.real() << "," << m_rawboxGG.imag() << ")"
                << " rawGZ=(" << m_rawboxGZ.real() << "," << m_rawboxGZ.imag() << ")"
                << " rawSub=(" << m_rawboxSub.real() << "," << m_rawboxSub.imag() << ")"
                << " coef=(" << m_rawcoef.real() << "," << m_rawcoef.imag() << ")"
                << " alpha=" << m_alpha << " MZ=" << m_MZ << " GZ=" << m_gZ
                << std::setprecision(10) << std::endl;
      msg_Error() << "@@@ SHPROP sp=" << m_sp
                << " propG=(" << m_propG.real() << "," << m_propG.imag() << ")"
                << " propZ=(" << m_propZ.real() << "," << m_propZ.imag() << ")"
                << std::endl;
      for (int j = 0; j <= 1; ++j)
        msg_Error() << "@@@ SHFF j=" << j
                  << " TC=(" << m_TC[j].real() << "," << m_TC[j].imag() << ")"
                  << " UC=(" << m_UC[j].real() << "," << m_UC[j].imag() << ")"
                  << std::endl;
      
      msg_Error() << "@@@ SHPHEL n=" << m_PhoHel.size()
                << " hel0=" << (m_PhoHel.empty()? 0 : m_PhoHel[0])
                << " nisr=" << m_isrphotons.size()
                << " nfsr=" << m_fsrphotons.size()
                << " nparts=" << m_nparts
                << std::endl;
      for (int h = 0; h <= 1; ++h) {
        const Complex sf(Sfactor(m_pceex[0], m_pceex[1], m_allphotons[0],
                                 h == 0 ? 1 : -1));
        msg_Error() << "@@@ SHSINI h=" << h
                  << " S=(" << sf.real() << "," << sf.imag() << ")"
                  << " |S|=" << std::abs(sf) << std::endl;
      }
      for (int j1 = 0; j1 <= 1; ++j1)
       for (int j2 = 0; j2 <= 1; ++j2)
        for (int j3 = 0; j3 <= 1; ++j3)
         for (int j4 = 0; j4 <= 1; ++j4) {
           const Complex tt(m_Tamp[j1][j2][j3][j4]);
           const Complex uu(m_Uamp[j1][j2][j3][j4]);
           const Complex bc(bref.m_A[j1][j2][j3][j4]);
           if (std::abs(tt)==0. && std::abs(uu)==0. && std::abs(bc)==0.) continue;
           msg_Error() << "@@@ SHSPIN " << j1 << j2 << j3 << j4
                     << " TT=(" << tt.real() << "," << tt.imag() << ")"
                     << " UU=(" << uu.real() << "," << uu.imag() << ")"
                     << " Born=(" << bc.real() << "," << bc.imag() << ")"
                     << std::endl;
         }
    }

    
    msg_Error() << std::setprecision(10);
    for (int j1 = 0; j1 <= 1; ++j1)
      for (int j2 = 0; j2 <= 1; ++j2)
        for (int j3 = 0; j3 <= 1; ++j3)
          for (int j4 = 0; j4 <= 1; ++j4) {
            const Complex a0(m_AmpExpo0.m_A[j1][j2][j3][j4]);
            const Complex a1(m_AmpExpo1.m_A[j1][j2][j3][j4]);
            if (std::abs(a0) == 0. && std::abs(a1) == 0.) continue;
            const Complex bo(m_snapBorn.m_A[j1][j2][j3][j4]);
            const Complex vi(m_snapVirt.m_A[j1][j2][j3][j4]);
            const Complex re(m_snapReal.m_A[j1][j2][j3][j4]);
            msg_Error() << "@@@ SHPART " << j1 << j2 << j3 << j4
                      << " born=(" << bo.real() << "," << bo.imag() << ")"
                      << " virt=(" << vi.real() << "," << vi.imag() << ")"
                      << " real=(" << re.real() << "," << re.imag() << ")"
                      << " virt/born=" << (std::abs(bo)>0.? std::abs(vi)/std::abs(bo):0.)
                      << " real/born=" << (std::abs(bo)>0.? std::abs(re)/std::abs(bo):0.)
                      << std::endl;
            msg_Error() << "@@@ SHAMP " << j1 << j2 << j3 << j4
                      << " A0=(" << a0.real() << "," << a0.imag() << ")"
                      << " A1=(" << a1.real() << "," << a1.imag() << ")"
                      << " A1/A0=(" << (std::abs(a0)>0.? (a1/a0).real() : 0.)
                      << "," << (std::abs(a0)>0.? (a1/a0).imag() : 0.) << ")"
                      << std::endl;
          }
    msg_Out() << "=== Sherpa CEEX point written to ceex_point.dat ===\n"
              << "  rho0 (O(alpha^0)) = " << m_result0 << "\n"
              << "  rho1 (O(alpha^1)) = " << m_result << "\n"
              << "  rho1/rho0 - 1     = " << (m_result/m_result0 - 1.) << "\n"
              << "  compare against the KKMC harness's\n"
              << "    (RhoExp1 - RhoExp0)/RhoExp0 on the same point.\n"
              << "  sqrt(s) of this point = "
              << (m_momenta[0]+m_momenta[1]).Mass()
              << "   (momentum balance " << worst << ")\n"
              << "  m_sp (propagator s used)  = " << m_sp << "\n"
              << "  (P-k)^2 = s-channel momentum for ISR = "
              << (m_momenta[4]+m_momenta[5]).Abs2()
              << "   <-- CEEX eq.(one-photon): B(X) with X = P-k_1\n"
              << "  photon x = 2E/sqrt(s) = " << xgam
              << ",  E_gamma = " << m_allphotons[0][0] << " GeV\n";
  }

  static const bool cxchk(getenv("SHERPA_CEEX_COMIX")!=NULL);
  if (cxchk && m_bornmomenta.size() >= 4 && p_bornproc) {
    msg_Debugging()<<METHOD<<"(): entering Comix Born comparison\n";
    Vec4D_Vector bp(m_bornmomenta);
    Poincare bcms(bp[0] + bp[1]);
    for (size_t i(0); i < bp.size(); ++i) bcms.Boost(bp[i]);
    Amplitude cx, hand;
    const double sp_save(m_sp);
    m_sp = (bp[2] + bp[3]).Abs2();
    MakeProp();
    BornAmplitude(bp, hand);
    const bool ok(ComixBornAmplitude(bp, cx));
    m_sp = sp_save;
    MakeProp();
    if (ok) {
      double sh(0.), sc(0.), worst(0.);
      for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
          for (int c = 0; c <= 1; ++c)
            for (int d = 0; d <= 1; ++d) {
              const Complex H(m_e*m_e*hand.m_A[a][b][c][d]), C(cx.m_A[a][b][c][d]);
              sh += std::norm(H); sc += std::norm(C);
              const double den(std::abs(H) + std::abs(C));
              if (den > 0.) worst = Max(worst, std::abs(H - C)/den);
            }
      std::cerr<<"@@@ CEEXCMP sum_hand="<<sh<<" sum_comix="<<sc
               <<" ratio="<<(sc!=0.? sh/sc : 0.)
               <<" worst_elem_reldiff="<<worst
               <<" t="<<m_tinv<<" cth="<<cos(bp[2].Theta())<<std::endl;
      static int nrow(0);
      if (nrow++ < 3)
        for (int a = 0; a <= 1; ++a)
          for (int b = 0; b <= 1; ++b)
            for (int c = 0; c <= 1; ++c)
              for (int d = 0; d <= 1; ++d) {
                const Complex H(m_e*m_e*hand.m_A[a][b][c][d]), C(cx.m_A[a][b][c][d]);
                std::cerr<<"@@@ CEEXHEL "<<a<<b<<c<<d
                         <<" absH="<<std::abs(H)<<" absC="<<std::abs(C)
                         <<" H=("<<H.real()<<","<<H.imag()<<")"
                         <<" C=("<<C.real()<<","<<C.imag()<<")"<<std::endl;
              }
    }
  }

  if (m_checkxs && m_bornmomenta.size() >= 4) {
    Vec4D_Vector bp(m_bornmomenta);
    Poincare bcms(bp[0] + bp[1]);
    for (size_t i(0); i < bp.size(); ++i) bcms.Boost(bp[i]);
    const double sp_save(m_sp);
    m_sp = (bp[2] + bp[3]).Abs2();
    MakeProp();
    MakePropT(bp);
    Amplitude b;
    BornAmplitude(bp, b);
    double sum(0.);
    for (int j1 = 0; j1 <= 1; ++j1)
      for (int j2 = 0; j2 <= 1; ++j2)
        for (int j3 = 0; j3 <= 1; ++j3)
          for (int j4 = 0; j4 <= 1; ++j4) {
            const Complex a(m_e * m_e * b.m_A[j1][j2][j3][j4]);
            sum += std::real(a * conj(a));
          }
    m_sp = sp_save;
    MakeProp();
    MakePropT(m_pceex);
    AccumulateBornShape(cos(bp[2].Theta()), sum / 4.);
  }
}
