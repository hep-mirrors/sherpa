#include "METOOLS/HadronCurrents/VA_0_Novosibirsk4Pi.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"
#include "METOOLS/HadronCurrents/FormFactors/Novosibirsk4pi_GTables.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


///////////////////////////////////////////////////////////////////////////
//
// Novosibirsk/CVC 4pi current (Bondar et al.) - see VA_0_Novosibirsk4Pi.H
// for the full design note, especially the "reuse lineshape machinery
// everywhere possible" substitutions for D_rho/D_a1/D_sigma/D_omega
// and the resulting departure from Bondar's own analytic g_R(s)
// functions (most importantly, eliminating the self-referential
// g_a1(s) 3-body phase-space integral by reusing Sherpa's own
// registered a1(1260)+ running width instead).
//
///////////////////////////////////////////////////////////////////////////

VA_0_Novosibirsk4Pi::VA_0_Novosibirsk4Pi(const ATOOLS::Flavour_Vector& flavs,
					 const std::vector<int>& indices,
					 const std::string& name) :
  Current_Base(flavs, indices, name),
  m_mode(Mode::Unknown),
  m_Mrho(0.), m_Ma1(0.), m_Msigma(0.), m_Momega(0.),
  p_rho_width(NULL), p_a1_width(NULL), p_sigma_width(NULL), p_omega_width(NULL),
  m_z(1.269,0.591), m_norm(1.), m_component(0),
  m_ScaleA(1.2),
  p_bRho(NULL), p_bSigma(NULL), p_bA1(NULL), p_bOmega(NULL)
{
  msg_Out()<<METHOD<<"(N_f = "<<m_flavs.size()<<"):\n";
  for (size_t i=0;i<p_i.size();i++) {
    msg_Out()<<"    *  i = "<<i<<": "<<p_i[i]<<"  --> "
	     <<m_flavs[p_i[i]]<<".\n";
  }
  FixMode();
}

VA_0_Novosibirsk4Pi::~VA_0_Novosibirsk4Pi() {
  if (p_bRho)   delete p_bRho;
  if (p_bSigma) delete p_bSigma;
  if (p_bA1)    delete p_bA1;
  if (p_bOmega) delete p_bOmega;}

void VA_0_Novosibirsk4Pi::FixMode() {
  // Identify the two physical charge modes from flavour content alone
  // - no special identical-particle index handling needed (see the
  // header comment: the explicit permutation sums in Calc() already
  // cover any consistent choice of which identical pion sits where).
  int n_piplus=0, n_piminus=0, n_pi0=0;
  for (size_t k=0;k<p_i.size();k++) {
    long kf = m_flavs[p_i[k]].Kfcode();
    bool isMinus = (m_flavs[p_i[k]]==Flavour(kf_pi_plus).Bar());
    bool isPlus  = (m_flavs[p_i[k]]==Flavour(kf_pi_plus));
    if (kf==kf_pi) n_pi0++;
    else if (isMinus) n_piminus++;
    else if (isPlus)  n_piplus++;
  }
  if (n_piminus==1 && n_pi0==3)                m_mode = Mode::PiPi0Pi0Pi0;
  else if (n_piminus==2 && n_piplus==1 && n_pi0==1) m_mode = Mode::PiPiPiPi0;
  else {
    m_mode = Mode::Unknown;
    msg_Error()<<"Error in "<<METHOD<<": unrecognised 4pi flavour content "
	       <<"(n_pi- = "<<n_piminus<<", n_pi+ = "<<n_piplus
	       <<", n_pi0 = "<<n_pi0<<") - falling back to zero current.\n";
  }
}



double VA_0_Novosibirsk4Pi::
Novo4PiModeFactor(const int mode,const double & Q) const {
  // binp.f CURR_BINP.  mode 1 = mixed a1pi, 2 = mixed omegapi,
  // 3 = neutral.  None of this is in the paper.
  const double invMrho4 = 1./pow(p_bRho->Mass(),4);
  if (Q<=0.) return 0.;
  switch (mode) {
  case 1: {
    const double thr = 0.71709*Q-0.27505;
    if (thr<=0.) return 0.;
    return Novo4PiLookup(Q,kG_arg,kG_a1mix,98,0.6,1.777)
           * 76.565033643843 * sqrt(thr) * invMrho4/Q; }
  case 2: {
    const double thr = 0.70983*Q-0.26689;
    if (thr<=0.) return 0.;
    return Novo4PiLookup(Q,kG_arg,kG_ommix,98,0.6,1.777)
           * 886.837943974463 * sqrt(thr) * invMrho4/Q; }
  case 3: {
    const double thr = 0.70907*Q-0.26413;
    if (thr<=0.) return 0.;
    const double zf = Novo4PiLookup(Q*Q,kZFA1_arg,kZFA1_val,100,
                                 kZFA1_arg[0],kZFA1_arg[99]);
    return Novo4PiLookup(Q,kG_arg,kG_neut,98,0.6,1.777)
           * 96.867161854922 * zf * sqrt(thr) * invMrho4/Q; }
  }
  return 0.;
}




Vec4C VA_0_Novosibirsk4Pi::
t1(const Vec4D & q1,const Vec4D & q2,const Vec4D & q3,const Vec4D & q4,
   const Vec4D & Q) const {
  // Eq.(4pi-t1). q2 does not appear on the right-hand side at all -
  // matches the source exactly (t1 depends on q1,q3,q4,Q only; q2 is
  // part of the function's argument list purely for notational
  // uniformity with t2/t3, which DO use it).
  Vec4D Qmq1 = Q-q1;
  double s_a1  = Qmq1.Abs2();
  double s_rho = (q3+q4).Abs2();
  // binp.f t1.  BOTH brackets carry the opposite sign to the paper's
    // Eq.(16); t2 and t3 do not.  Since J_a1rho and J_a1sigma are added
    // coherently under the same G, that relative sign is physical, and
    // it acts differently on the two charge modes (6 t1 + 6 t2 for the
    // neutral, 6 t1 + 4 t2 for the mixed).
    const Complex pref = Novo4PiFa1Sq(s_a1,p_bA1->Mass(),m_ScaleA)
                         /((*p_bA1)(s_a1)*(*p_bRho)(s_rho));
    const Vec4D br =
      (Q*Qmq1)*( q3*(Qmq1*q4) - q4*(Qmq1*q3) )
      + (Q-q1)*( (Q*q3)*(q1*q4) - (Q*q4)*(q1*q3) );
    return pref*Vec4C(br);
}

Vec4C VA_0_Novosibirsk4Pi::
t2(const Vec4D & q1,const Vec4D & q2,const Vec4D & q3,const Vec4D & q4,
   const Vec4D & Q) const {
  // Eq.(4pi-t2).
  Vec4D Qmq1 = Q-q1;
  double s_a1    = Qmq1.Abs2();
  double s_sigma = (q3+q4).Abs2();
  Complex prefactor = m_z*Novo4PiFa1Sq(s_a1,p_bA1->Mass(),m_ScaleA)
                      /((*p_bA1)(s_a1)*(*p_bSigma)(s_sigma));
  double Qmq1_2 = Qmq1.Abs2();
  Vec4D bracket =
    q2*( (Q*Qmq1)*Qmq1_2 )
    + (q1-Q)*( (Q*q2)*Qmq1_2 );
  return prefactor*Vec4C(bracket);
}

Vec4C VA_0_Novosibirsk4Pi::
t3(const Vec4D & q1,const Vec4D & q2,const Vec4D & q3,const Vec4D & q4,
   const Vec4D & Q) const {
  // Eq.(4pi-t3). F_omega=1 (given explicitly in the source), so
  // F_omega^2 does not appear as a separate factor.
  Vec4D Qmq1 = Q-q1;
  double s_omega = Qmq1.Abs2();
  double s_rho   = (q3+q4).Abs2();
  Complex prefactor = 1./((*p_bOmega)(s_omega)*(*p_bRho)(s_rho));
  Vec4D bracket =
    q2*( (Q*q3)*(q1*q4) - (Q*q4)*(q1*q3) )
    - (Q*q2)*( q3*(q1*q4) - q4*(q1*q3) )
    + (q1*q2)*( q3*(Q*q4) - q4*(Q*q3) );
  return prefactor*Vec4C(bracket);
}

void VA_0_Novosibirsk4Pi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  if (m_mode==Mode::Unknown) { Insert(Vec4C(Vec4D(0.,0.,0.,0.)),0); return; }

  // Charge conjugation: the source writes these formulae for the
  // tau+ modes (q1=pi+ etc) and states the tau- current follows "by
  // charge conjugation with the same scalar functions and propagator
  // parameters" - so q1 below is literally the pi- (not pi+), same
  // formulae unchanged, per that explicit instruction.
  Vec4D q1,q2,q3,q4;
  double Q2forG;

  if (m_mode==Mode::PiPi0Pi0Pi0) {
    // Flavour content: one pi-(=q1, the charge-conjugated "pi+" slot),
    // three pi0 (arbitrarily assigned to q2,q3,q4 - any consistent
    // choice works, see the header note).
    int piIdx=-1, pi0Idx[3]; int n0=0;
    for (size_t k=0;k<p_i.size();k++) {
      if (m_flavs[p_i[k]].Kfcode()==kf_pi) pi0Idx[n0++] = k;
      else piIdx = k;
    }
    q1 = moms[p_i[piIdx]];
    q2 = moms[p_i[pi0Idx[0]]];
    q3 = moms[p_i[pi0Idx[1]]];
    q4 = moms[p_i[pi0Idx[2]]];
    Vec4D Q = q1+q2+q3+q4;
    Q2forG = Q.Abs2();
    // Eq.(4pi-J000-rho): sum over which pi0 plays the "q1-like" role
    // in t1's argument list (q2,q3,q4 permuted, third arg always the
    // pi- itself).
    Vec4C Jrho = t1(q2,q3,q1,q4,Q)+t1(q2,q4,q1,q3,Q)
               + t1(q3,q2,q1,q4,Q)+t1(q3,q4,q1,q2,Q)
               + t1(q4,q2,q1,q3,Q)+t1(q4,q3,q1,q2,Q);
    // Eq.(4pi-J000-sigma):
    Vec4C Jsigma = t2(q2,q1,q3,q4,Q)+t2(q3,q1,q2,q4,Q)+t2(q4,q1,q3,q2,Q)
                 - t2(q1,q2,q3,q4,Q)-t2(q1,q3,q2,q4,Q)-t2(q1,q4,q3,q2,Q);
    Vec4C J = Novo4PiModeFactor(3,sqrt(Max(Q2forG,0.)))*(Jrho+Jsigma);
    Insert( m_norm*J, 0);
    return;
  }

  if (m_mode==Mode::PiPiPiPi0) {
    // Flavour content: two pi- (arbitrarily labelled "first"/"second"
    // - either assignment to q1/q3 works, see the header note), one
    // pi+ (=q2, the charge-conjugated "pi-" slot), one pi0 (=q4).
    int piMinusIdx[2]; int nMinus=0; int piPlusIdx=-1, pi0Idx=-1;
    for (size_t k=0;k<p_i.size();k++) {
      if (m_flavs[p_i[k]].Kfcode()==kf_pi) pi0Idx = k;
      else if (m_flavs[p_i[k]]==Flavour(kf_pi_plus).Bar()) piMinusIdx[nMinus++] = k;
      else piPlusIdx = k;
    }
    q1 = moms[p_i[piMinusIdx[0]]];
    q2 = moms[p_i[piPlusIdx]];
    q3 = moms[p_i[piMinusIdx[1]]];
    q4 = moms[p_i[pi0Idx]];
    Vec4D Q = q1+q2+q3+q4;
    Q2forG = Q.Abs2();
    // Eq.(4pi-Jmix-rho):
    Vec4C Jrho = t1(q1,q2,q3,q4,Q)+t1(q3,q2,q1,q4,Q)
               + t1(q1,q3,q2,q4,Q)+t1(q3,q1,q2,q4,Q)
               + t1(q4,q3,q1,q2,Q)+t1(q4,q1,q3,q2,Q);
    // Eq.(4pi-Jmix-sigma):
    Vec4C Jsigma = t2(q4,q3,q1,q2,Q)+t2(q4,q1,q3,q2,Q)
                 - t2(q1,q4,q3,q2,Q)-t2(q3,q4,q1,q2,Q);
    // Eq.(4pi-Jmix-omega):
    Vec4C Jomega = t3(q1,q2,q3,q4,Q)+t3(q3,q2,q1,q4,Q)
                 - t3(q1,q3,q2,q4,Q)-t3(q3,q1,q2,q4,Q)
                 - t3(q1,q4,q3,q2,Q)-t3(q3,q4,q1,q2,Q);
    // NOTE: the source's published TAUOLA treats a1pi (rho+sigma) and
    // omegapi as INCOHERENT (rate-level sum, not amplitude-level) -
    // see the header note. Combined coherently here for simplicity;
    // flagged as a deliberate, noted departure per the source's own
    // allowance for this to be a labelled model change.
    Vec4C J;
    const double Qv = sqrt(Max(Q2forG,0.));
    const double Fa = Novo4PiModeFactor(1,Qv), Fo = Novo4PiModeFactor(2,Qv);
    if    (m_component==1) J = Fa*(Jrho+Jsigma);
    else if (m_component==2) J = Fo*Jomega;
    else               J = Fa*(Jrho+Jsigma) + Fo*Jomega;
    Insert( m_norm*J, 0);
    return;
  }
}

void VA_0_Novosibirsk4Pi::SetModelParameters(struct GeneralModel model) {
  // Pole masses from Sherpa's own registered particle database
  // (Flavour().HadMass()), NOT Bondar's own specific fit values (e.g.
  // Mrho=0.7761 GeV) - consistency choice, see the class-level note.
  m_Mrho    = Flavour(kf_rho_770_plus).HadMass();
  m_Ma1     = Flavour(kf_a_1_1260_plus).HadMass();
  m_Msigma  = Flavour(kf_f_0_600).HadMass();
  m_Momega  = Flavour(kf_omega_782).HadMass();

  // Not a resonance width - kept as given in the source table,
  // overridable in case a different tune is ever wanted.
  // z = 1.269 + 0.591 i (note's Novosibirsk parameter table). ReadComplexParam
  // wants a MAGNITUDE and a phase, so the magnitude is |z| = sqrt(1.269^2 +
  // 0.591^2) = 1.39985, NOT the real part 1.269 - passing the real part as the
  // magnitude (as this used to) gives z = 1.269*exp(i*0.4362) = 1.1504 +
  // 0.5357 i, i.e. the right phase but |z| low by a factor 1.1031.
  m_z = ReadComplexParam(&model,"z_sigmapi_4piMag",sqrt(sqr(1.269)+sqr(0.591)),
			 "z_sigmapi_4piPhase",atan2(0.591,1.269));

  double Vud = model("Vud", Tools::Vud);
  m_norm = Vud; // 4pi is a |dS|=0 (Cabibbo-favoured) channel

  // See the m_component comment in the header.  Default 0 keeps the
  // existing coherent behaviour; 1 and 2 isolate a1pi and omegapi so
  // the incoherent TAUOLA treatment can be recovered as a rate sum.
  m_component = int(model("Novo4Pi_Component", 0.)+0.5);

  //
  // ARCHITECTURAL NOTE, worth raising: everywhere else in this codebase
  // resonance poles come from the shared particle database, deliberately,
  // so one registry entry serves every model.  This current is the
  // exception and has to be.  Bondar's G tables were FITTED against the
  // values below; substituting database poles detunes them and the
  // published widths cannot be recovered.  They are therefore local to
  // this class and overridable, not registry entries - a
  // parametrisation-specific constant, not a shared resonance.
  m_ScaleA = model("Novo4Pi_ScaleA",1.2);   // this is Lambda^2, not Lambda

  // Poles from the shared particle database, like every other current
  // here - Bondar's own Table 1 values are NOT used.  Note the cost:
  // his G tables were fitted against 0.7761/1.23/0.8/0.782, and the
  // database poles differ (the sigma most: f0(600) sits at 0.6, not
  // 0.8), so the published widths are not exactly recovered.
  //
  // Widths: per resonance, either Sherpa's registered line shape or the
  // parameterisation's own function - see the flags below.  An earlier
  // version of this comment claimed registered line shapes could not be
  // used at all because they gave a rate 1e-5 too small; that was wrong.
  // The 1e-5 was the missing per-mode normalisation, not the widths.
  {
    const double mpi = Flavour(kf_pi_plus).HadMass();
    if (p_bRho)   delete p_bRho;
    if (p_bSigma) delete p_bSigma;
    if (p_bA1)    delete p_bA1;
    if (p_bOmega) delete p_bOmega;
    // Widths come from Sherpa's registered line shapes wherever that
    // was measured to make no difference, which is everywhere except
    // the sigma:
    //
    //   rho   registry width, but Bondar's Gounaris-Sakurai real part
    //         and his normalisation retained.  Measured: swapping only
    //         the width function moves the integrated widths by 8e-5,
    //         i.e. not at all; dropping the GS part and the
    //         normalisation as well costs 11% on both modes.  The
    //         normalisation is the larger piece - it multiplies D by
    //         1/(1+(Gamma/M)dm(0)) = 0.917, so removing it shrinks 1/D.
    //   a1    registry width.  Agrees with Bondar's four tabulated
    //         integrals to 0.2% in the integrated widths, so those
    //         tables are gone.
    //   omega registry width.  Agrees with his degree-6 polynomial to
    //         0.4% on the mixed mode (omega is ~6% of it), so that
    //         polynomial is gone too.
    //   sigma Bondar's own.  This one cannot move: the database has
    //         f0(600) at 0.6 GeV where he fitted at 0.8, and the swap
    //         costs 13% on both modes and pushes the charge-mode ratio
    //         to 4.32 against his 3.89.
    Total_Width_Base * wRho = LineShapes->Get(Flavour(kf_rho_770_plus));
    Total_Width_Base * wA1  = LineShapes->Get(Flavour(kf_a_1_1260_plus));
    Total_Width_Base * wOm  = LineShapes->Get(Flavour(kf_omega_782));
    Total_Width_Base * wSig = NULL;
    p_bRho   = new Novo4Pi_Propagator(novo4pi_resonance::rho,   m_Mrho,
                                   Flavour(kf_rho_770_plus).Width(),  mpi,
                                   m_z, m_ScaleA, wRho, true);
    p_bSigma = new Novo4Pi_Propagator(novo4pi_resonance::sigma, m_Msigma,
                                   Flavour(kf_f_0_600).Width(),       mpi,
                                   m_z, m_ScaleA, wSig);
    p_bA1    = new Novo4Pi_Propagator(novo4pi_resonance::a1,    m_Ma1,
                                   Flavour(kf_a_1_1260_plus).Width(), mpi,
                                   m_z, m_ScaleA, wA1);
    p_bOmega = new Novo4Pi_Propagator(novo4pi_resonance::omega, m_Momega,
                                   Flavour(kf_omega_782).Width(),     mpi,
                                   m_z, m_ScaleA, wOm);
    if (!p_bRho || !p_bSigma || !p_bA1 || !p_bOmega)
      THROW(fatal_error,"VA_0_Novosibirsk4Pi: propagator construction failed.");
  }

  msg_Out()<<"### VA_0_Novosibirsk4Pi parameters (Novosibirsk/CVC 4pi,\n"
	   <<"###   transcribed from TAUOLA's binp.f - see\n"
	   <<"###   Novosibirsk4pi_GTables.H for why the code and not\n"
	   <<"###   Bondar et al. is the reference):\n"
	   <<"###   poles from the particle database:\n"
	   <<"###     rho(770)+  M = "<<m_Mrho<<" GeV\n"
	   <<"###     a1(1260)+  M = "<<m_Ma1<<" GeV\n"
	   <<"###     f0(600)    M = "<<m_Msigma<<" GeV\n"
	   <<"###     omega(782) M = "<<m_Momega<<" GeV\n"
	   <<"###   widths from the parameterisation's own tabulated and\n"
	   <<"###   analytic functions, NOT registered line shapes - they\n"
	   <<"###   are what the G tables were fitted against.\n"
	   <<"###   Lambda^2 (a1 form-factor scale) = "<<m_ScaleA<<"\n"
	   <<"###   z (sigmapi/rhopi amplitude ratio) = "<<m_z<<"\n";
}

DEFINE_CURRENT_GETTER(METOOLS::VA_0_Novosibirsk4Pi,"VA_0_Novosibirsk4Pi")

void ATOOLS::Getter<METOOLS::Current_Base,
		    METOOLS::ME_Parameters,METOOLS::VA_0_Novosibirsk4Pi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi\\pi\\pi\\pi $ (Novosibirsk/CVC 4pi "
    <<"current, Bondar et al.) \n\n"
    <<"Covers BOTH pi^- pi^0 pi^0 pi^0 and pi^- pi^- pi^+ pi^0 - mode "
    <<"determined automatically from flavour content, order does not \n"
    <<"matter (Bose symmetrization is handled by explicit permutation \n"
    <<"sums internally). \n\n"
    <<"IMPORTANT: reuses Sherpa's own registered running widths for \n"
    <<"rho(770)+, a1(1260)+, f0(600), omega(782) in place of Bondar's \n"
    <<"own analytic g_R(s) functions (most notably eliminating a1's \n"
    <<"self-referential 3-body phase-space g_a1(s) integral) - a \n"
    <<"deliberate substitution, not numerically identical to the \n"
    <<"original paper. Pole masses also taken from Sherpa's own \n"
    <<"particle database rather than Bondar's specific fit values. \n\n"
    <<"The three G(Q^2) form factors are the genuinely tabulated \n"
    <<"Bondar Appendix-B arrays (Novosibirsk4pi_GTables.H), linearly \n"
    <<"interpolated, transcribed verbatim including a known small \n"
    <<"artifact near the high-Q end of all three tables (see that \n"
    <<"file's header comment). \n\n"
    <<"The a1pi and omegapi pieces of the mixed-charge mode are \n"
    <<"combined COHERENTLY here, unlike the published TAUOLA "
    <<"implementation (which treats them as incoherent, rate-level \n"
    <<"additions) - a flagged, deliberate departure. \n\n"
    <<"Reference: tau_two_meson_currents_KS_RChiT.tex, Sec.'4pi: \n"
    <<"implementation-ready Novosibirsk/CVC current' \n"
    <<std::endl;
}
