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
  m_Lambda(1.2), m_z(1.269,0.591), m_kappaF0(0.,0.), m_norm(1.)
{
  msg_Out()<<METHOD<<"(N_f = "<<m_flavs.size()<<"):\n";
  for (size_t i=0;i<p_i.size();i++) {
    msg_Out()<<"    *  i = "<<i<<": "<<p_i[i]<<"  --> "
	     <<m_flavs[p_i[i]]<<".\n";
  }
  FixMode();
}

VA_0_Novosibirsk4Pi::~VA_0_Novosibirsk4Pi() {}

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

Complex VA_0_Novosibirsk4Pi::
DR(const double & s,Total_Width_Base * width,const double & M) const {
  // Bare propagator 1/(s-M^2+iM*Gamma(s)), Gamma(s) from Sherpa's own
  // registered running width - see the class-level "reuse lineshape
  // machinery" note for why this replaces Bondar's own g_R(s)/g_R(M^2)
  // construction (Eq.4pi-DR).
  if (width==NULL) return Complex(0.,0.);
  double Gamma_s = (*width)(s);
  return 1./Complex(s-sqr(M), M*Gamma_s);
}

Complex VA_0_Novosibirsk4Pi::DRrho(const double & s) const {
  // rho(770) + kappa*f0(500)/sigma admixture - NEW addition (not in
  // Bondar's original current), requested to fill in the same
  // missing low-mass ππ scalar strength identified in the 3pi
  // channels. Reuses Sherpa's own registered f0(600) running width.
  Complex Drho = DR(s,p_rho_width,m_Mrho);
  if (p_sigma_width==NULL || m_kappaF0==Complex(0.,0.)) return Drho;
  return Drho + m_kappaF0*DR(s,p_sigma_width,m_Msigma);
}

double VA_0_Novosibirsk4Pi::Fa1sq(const double & s) const {
  // Eq.(4pi-Fa1): F_a1(q) = (1+Ma1^2/Lambda^2)/(1+q^2/Lambda^2).
  // Kept exactly as given - not a resonance width, no lineshape
  // substitute applies here.
  double num = 1.+sqr(m_Ma1)/sqr(m_Lambda);
  double den = 1.+s/sqr(m_Lambda);
  double Fa1 = num/den;
  return Fa1*Fa1;
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
  Complex prefactor = Fa1sq(s_a1)*DR(s_a1,p_a1_width,m_Ma1)*DRrho(s_rho);
  Vec4D bracket =
    (Q*Qmq1)*( q4*(Qmq1*q3) - q3*(Qmq1*q4) )
    + (Q-q1)*( (Q*q4)*(q1*q3) - (Q*q3)*(q4*q1) );
  return prefactor*Vec4C(bracket);
}

Vec4C VA_0_Novosibirsk4Pi::
t2(const Vec4D & q1,const Vec4D & q2,const Vec4D & q3,const Vec4D & q4,
   const Vec4D & Q) const {
  // Eq.(4pi-t2).
  Vec4D Qmq1 = Q-q1;
  double s_a1    = Qmq1.Abs2();
  double s_sigma = (q3+q4).Abs2();
  Complex prefactor = m_z*Fa1sq(s_a1)*DR(s_a1,p_a1_width,m_Ma1)*DR(s_sigma,p_sigma_width,m_Msigma);
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
  Complex prefactor = DR(s_omega,p_omega_width,m_Momega)*DRrho(s_rho);
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
  Complex Gval;

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
    Gval = Complex(NovosibirskGLookup(sqrt(Q2forG), kG_pi000), 0.);
    // Eq.(4pi-J000-rho): sum over which pi0 plays the "q1-like" role
    // in t1's argument list (q2,q3,q4 permuted, third arg always the
    // pi- itself).
    Vec4C Jrho = t1(q2,q3,q1,q4,Q)+t1(q2,q4,q1,q3,Q)
               + t1(q3,q2,q1,q4,Q)+t1(q3,q4,q1,q2,Q)
               + t1(q4,q2,q1,q3,Q)+t1(q4,q3,q1,q2,Q);
    // Eq.(4pi-J000-sigma):
    Vec4C Jsigma = t2(q2,q1,q3,q4,Q)+t2(q3,q1,q2,q4,Q)+t2(q4,q1,q3,q2,Q)
                 - t2(q1,q2,q3,q4,Q)-t2(q1,q3,q2,q4,Q)-t2(q1,q4,q3,q2,Q);
    Vec4C J = Gval*(Jrho+Jsigma);
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
    double Qval = sqrt(Q2forG);
    Complex Grho   = Complex(NovosibirskGLookup(Qval, kG_pimix), 0.);
    Complex Gomega = Complex(NovosibirskGLookup(Qval, kGomega_pimix), 0.);
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
    Vec4C J = Grho*(Jrho+Jsigma) + Gomega*Jomega;
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
  p_rho_width   = LineShapes->Get(Flavour(kf_rho_770_plus));
  p_a1_width    = LineShapes->Get(Flavour(kf_a_1_1260_plus));
  p_sigma_width = LineShapes->Get(Flavour(kf_f_0_600));
  p_omega_width = LineShapes->Get(Flavour(kf_omega_782));

  // Not resonance widths - kept exactly as given in the source table,
  // overridable in case a different tune is ever wanted.
  m_Lambda = model("Lambda_a1_4pi", 1.2);
  m_z = ReadComplexParam(&model,"z_sigmapi_4piMag",1.269,
			 "z_sigmapi_4piPhase",atan2(0.591,1.269));
  // f0(500)/sigma admixture in the rho tower used by t1/t3 - NEW
  // addition (request: "add a set of new models... plus the f0"),
  // not part of Bondar's original current. Default nonzero (unlike
  // the K-omega channel's h_R*t_R couplings, which default to 0
  // because there is NO known-reasonable value to guess at all) -
  // here a modest, same-order-of-magnitude guess relative to the
  // rho'/f0 admixture weights used elsewhere in this codebase
  // (deltaMag_pipi_f0 etc.) is at least plausible, but is still a
  // guess - retune against data.
  m_kappaF0 = ReadComplexParam(&model,"kappaF0_4piMag",-0.2,"kappaF0_4piPhase");

  double Vud = model("Vud", Tools::Vud);
  m_norm = Vud; // 4pi is a |dS|=0 (Cabibbo-favoured) channel

  msg_Out()<<"### VA_0_Novosibirsk4Pi parameters (lineshape-sourced "
	   <<"pole masses/widths, not Bondar's own specific values):\n"
	   <<"###   rho(770)+: M = "<<m_Mrho<<" GeV (running width via LineShapes)\n"
	   <<"###   a1(1260)+: M = "<<m_Ma1<<" GeV (running width via LineShapes,\n"
	   <<"###      replacing Bondar's own self-referential g_a1(s))\n"
	   <<"###   f0(600)/sigma: M = "<<m_Msigma<<" GeV (running width via LineShapes)\n"
	   <<"###   omega(782): M = "<<m_Momega<<" GeV (running width via LineShapes)\n"
	   <<"###   Lambda (a1 form factor scale) = "<<m_Lambda<<" GeV\n"
	   <<"###   z (sigmapi/rhopi amplitude ratio) = "<<m_z<<"\n"
	   <<"###   kappa_F0 (NEW rho-tower f0(500) admixture, not in "
	   <<"Bondar's original) = "<<m_kappaF0<<"\n";
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
