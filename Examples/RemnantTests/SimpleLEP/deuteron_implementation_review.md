# Deuteron formation in Sherpa: implementation review

This document puts every formula in the deuteron-formation code next to the publication it comes from, so that it can be checked by hand. Each item gives the code (copied from the file, with a link to the line), the source (paper, equation or table, link), and a verdict.

Line numbers refer to the working tree on branch `nucleon_recombination` as of 13 October 2026 (uncommitted changes); they move when the files are edited.

**Verdicts:** ✅ matches the source · ⚠️ deliberate deviation or interpretation, check · ❌ disagrees with the source · ⓘ no published source (our own choice)

## Points that need a reviewer's decision

| # | Item | Verdict | Section |
|---|---|---|---|
| 1 | Parameter $e$ of the $pp\to d\pi^+\pi^0$ fit is `2.866e16` in the code, **2.868e16** in Table VI of Dal & Raklev | ❌ (0.07 % effect) | [1.6](#16-pp--dpipi0-and-nn--dpi-pi0-eqs-13-12-table-vi) |
| 2 | Eq. (15) of Dal & Raklev is printed with $m^2_{d\pi}$; the code uses $m^2_{\pi\pi}$ (standard three-body kinematics) | ⚠️ (we read the paper as a typo) | [2.5](#25-three-body-final-state-dpipi-sec-iii-d-eq-15) |
| 3 | The `Exponential` model is a toy model; its relation to the Wigner function of arXiv:1905.01192, Eq. (11), is given in 3.3 | ⓘ toy model | [3.3](#33-exponential-model) |
| 4 | The impact-parameter condition of `HardDisk`/`Exponential` ($b<b_{\rm coal}$), defined with Sherpa's hadron positions | ⓘ our choice; literature uses a spatial separation instead | [3.4](#34-impact-parameter-b) |
| 5 | "Expected" spectra: every pair weighted with $P_i$ — the analytic limit of the paper's $N_{\rm samp}$ method, not the method itself | ⚠️ | [4.4](#44-expected-yields-weighted-pairs-sec-iii-g) |
| 6 | $f_{\rm ND}$ from ALICE's diffraction measurement, to compare with data normalised to *inelastic* events (the paper used V0AND-triggered data and Table IX) | ⚠️ | [5.3](#53-non-diffractive-fraction-f_nd) |
| 7 | $p_0$ from a single normalisation via $p_0\propto s^{1/3}$ (the paper's "roughly $p_0^3$" remark), not a scan | ⚠️ approximation | [5.4](#54-coalescence-momentum-from-the-normalisation) |
| 8 | Multiplicity classes: Sherpa's percentiles from its own non-diffractive INEL>0 events | ⚠️ | [4.8](#48-alice-multiplicity-classes-v0m-and-inel0) |

---

## 1. Cross-section fits — [`HADRON_RESCATTERING/XSecs/BaryonBaryon.C`](../../../HADRON_RESCATTERING/XSecs/BaryonBaryon.C)

Source: L. A. Dal and A. R. Raklev, *Alternative formation model for antideuterons from dark matter*, Phys. Rev. D 91 (2015) 123536, [arXiv:1504.07242](https://arxiv.org/abs/1504.07242), [doi:10.1103/PhysRevD.91.123536](https://doi.org/10.1103/PhysRevD.91.123536) — referred to as **D&R** below. All cross sections are in µb.

### 1.1 The variable $k$ (Eq. (1))

[`BaryonBaryon.C:75-77`](../../../HADRON_RESCATTERING/XSecs/BaryonBaryon.C#L75) (same in L138–140, L159–161, L176–178)
```cpp
    // k = |p1 - p2| in the c.m. frame, i.e. twice the c.m. momentum of each
    // nucleon (arXiv:1504.07242, after Eq. (1)); the fits are in kappa = k/GeV.
    double p_rel = 2.*sqrt(p1cm[1]*p1cm[1]+p1cm[2]*p1cm[2]+p1cm[3]*p1cm[3]);
```
D&R, text after Eq. (1): "where $k=|\vec p_{\bar p}-\vec p_{\bar n}|_{\rm COM}$"; Eqs. (7), (13), (14) use $\kappa=k/(1\,{\rm GeV})$. — ✅ (the factor 2 was missing before; fixed in this work)

### 1.2 $pn\to d\gamma$ (Eq. (7), Table II)

[`BaryonBaryon.C:64-92`](../../../HADRON_RESCATTERING/XSecs/BaryonBaryon.C#L64)
```cpp
    std::vector<double> aVec ={2.30346, -9.366346*10, 2.565390*std::pow(10,3),
                               -2.5594101*std::pow(10,4), 1.43513109*std::pow(10,5),
                               -5.0357289*std::pow(10,5), 1.14924802*std::pow(10,6),
                               -1.72368391*std::pow(10,6), 1.67934876*std::pow(10,6),
                               -1.01988855*std::pow(10,6), 3.4984035*std::pow(10,5),
                               -5.1662760*std::pow(10,4)};
    std::vector<double> bVec = {-5.1885, 2.9196};
...
if (p_rel < 1.28) {
    double kappa = p_rel; // kappa = k/(1 GeV), already dimensionless
    sigma = aVec[0] / kappa;          // a_{-1} * kappa^{-1}
    double kp = 1.0;
    for (size_t n = 1; n < aVec.size(); ++n) {
        sigma += aVec[n] * kp;        // a_{n-1} * kappa^{n-1}
        kp *= kappa;
    }
}
else {
    sigma = std::exp(-bVec[0]*p_rel - bVec[1]*p_rel*p_rel);
}
```
D&R Eq. (7): $\sigma(\kappa)/(1\,\mu{\rm b})=\sum_{n=-1}^{10}a_n\kappa^n$ for $\kappa<1.28$, $\exp(-b_1\kappa-b_2\kappa^2)$ for $\kappa\ge1.28$; parameters in Table II ($a_{-1}=2.30346,\ a_0=-9.366346\times10^{1},\dots,a_{10}=-5.1662760\times10^{4}$, $b_1=-5.1885$, $b_2=2.9196$). — ✅ all 14 parameters and the branch at $\kappa=1.28$ checked.

The $\kappa^{-1}$ term diverges at $k\to0$; D&R (text after Eq. (7)): "we therefore take care to restrict … $P\le1$". Done in [`Deuteron_Production.C:147`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L147), see 2.2. — ✅

### 1.3 $NN\to d\pi$ (Eq. (10), Table III; isospin Eqs. (8), (9))

[`BaryonBaryon.C:99-130`](../../../HADRON_RESCATTERING/XSecs/BaryonBaryon.C#L99)
```cpp
  double a = 170., b = 1.34, c= 1.77, d = 0.38, e = 0.096;// a is in micro barns
...
  double kallen = (s - std::pow(mDeuteron + mPion, 2))
              * (s - std::pow(mDeuteron - mPion, 2));
  double q = (kallen > 0.) ? std::sqrt(kallen) / (2.0 * sqrtS) : 0.0;
  double etaDefined = (q)/mPion;
...
  double denom = (c-std::exp(d*etaDefined));
  double sigma = a * std::pow(etaDefined,b)/(denom*denom + e);
...
double BaryonBaryon::nnToDPiMinus(...) { return NNToDPi(A,B); }
double BaryonBaryon::pnToDPiZero(...)  { return 0.5*NNToDPi(A,B); }
double BaryonBaryon::ppToDPiPlus(...)  { return NNToDPi(A,B); }
```
D&R Eq. (10): $\sigma(\eta)=a\eta^b/[(c-e^{d\eta})^2+e]$ with $\eta=q/m_{\pi^+}$, $q$ the pion momentum in the c.m. frame; Table III: $a=170\,\mu$b, $b=1.34$, $c=1.77$, $d=0.38$, $e=0.096$ (fit of H. Machner and J. Niskanen, Nucl. Phys. A 776 (2006) 172, [nucl-ex/0511027](https://arxiv.org/abs/nucl-ex/0511027)). Eq. (8): $\sigma_{pn\to d\pi^0}=\frac12\sigma_{pp\to d\pi^+}$; Eq. (9): $\sigma_{nn\to d\pi^-}=\sigma_{pp\to d\pi^+}$. Footnote 3: use $m_{\pi^+}$ in $\eta$ also for $pn\to d\pi^0$ and $nn\to d\pi^-$. — ✅ (`mPion` = 0.139570 GeV = $m_{\pi^+}$; $q$ from the Källén function with $m_d$ and $m_{\pi^+}$)

### 1.4 $pn\to d\pi^0\pi^0$ (Eq. (13), Table IV)

[`BaryonBaryon.C:154,164`](../../../HADRON_RESCATTERING/XSecs/BaryonBaryon.C#L154)
```cpp
  std::vector<double> piZeroPiZero = {2.855*std::pow(10,6), 1.311*10, 2.961*1000, 5.572, 1.461*std::pow(10,6)};
...
  sigma = piZeroPiZero[0]*std::pow(p_rel,piZeroPiZero[1])/(std::pow( (piZeroPiZero[2] - std::exp(piZeroPiZero[3]*p_rel)),2 ) + piZeroPiZero[4]);
```
D&R Eq. (13): $\sigma(\kappa)=a\kappa^b/[(c-e^{d\kappa})^2+e]$; Table IV: $a=2.855\times10^6\,\mu$b, $b=1.311\times10^1$, $c=2.961\times10^3$, $d=5.572$, $e=1.461\times10^6$. — ✅

### 1.5 $pn\to d\pi^+\pi^-$ (Eq. (14), Table V)

[`BaryonBaryon.C:170-171,181-182`](../../../HADRON_RESCATTERING/XSecs/BaryonBaryon.C#L170)
```cpp
  std::vector<double> piPlusPiMinusOne = {6.465*std::pow(10,6), 1.051*10, 1.979*1000, 5.363, 6.045*std::pow(10,5)};
  std::vector<double> piPlusPiMinusTwo = {2.549*std::pow(10,15), 1.657*10, 2.330*std::pow(10,7), 1.119*10, 2.868*std::pow(10,16)};
...
  sigma = One[0]*pow(k,One[1])/(pow(One[2]-exp(One[3]*k),2)+One[4]) + Two[0]*pow(k,Two[1])/(pow(Two[2]-exp(Two[3]*k),2)+Two[4]);
```
(last line abbreviated; see the file). D&R Eq. (14): sum of two terms of the form (13); Table V: $a_1=6.465\times10^6$, $b_1=10.51$, $c_1=1.979\times10^3$, $d_1=5.363$, $e_1=6.045\times10^5$; $a_2=2.549\times10^{15}$, $b_2=16.57$, $c_2=2.330\times10^7$, $d_2=11.19$, $e_2=2.868\times10^{16}$. — ✅

### 1.6 $pp\to d\pi^+\pi^0$ and $nn\to d\pi^-\pi^0$ (Eqs. (13), (12), Table VI)

[`BaryonBaryon.C:133,143,147-150`](../../../HADRON_RESCATTERING/XSecs/BaryonBaryon.C#L133)
```cpp
  std::vector<double> piPlusPiZero = {5.099*std::pow(10,15), 1.656*10, 2.333*std::pow(10,7), 1.133*10, 2.866*std::pow(10,16)};
...
  double BaryonBaryon::nnToPiMinusPiZero(...) { return ppToPiPlusPiZero(A,B); }
```
D&R Table VI: $a=5.099\times10^{15}\,\mu$b, $b=1.656\times10^1$, $c=2.333\times10^7$, $d=1.133\times10^1$, **$e=2.868\times10^{16}$**; Eq. (12): $\sigma_{nn\to d\pi^-\pi^0}=\sigma_{pp\to d\pi^+\pi^0}$.
— ❌ **$e$ is `2.866e16` in the code, `2.868e16` in Table VI.** Changes $\sigma$ by at most 0.07 % near the peak. Suggested fix: `2.868*std::pow(10,16)`.

---

## 2. Formation per nucleon pair — [`HADRON_RESCATTERING/Main/Deuteron_Production.C`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C)

### 2.1 Channels and kinematic thresholds (Table I; Secs. III C, III D)

[`Deuteron_Production.C:118-122`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L118)
```cpp
double Deuteron_Coalescence::Sigma(const Deuteron_Channel & channel,
                                   Particle * A,Particle * B) {
  // The fits are not all zero below threshold (e.g. the d pi pi ones),
  // so enforce the kinematic thresholds here, as in the paper.
  if ((A->Momentum()+B->Momentum()).Mass()<=Threshold(channel)) return 0.;
```
D&R Table I lists the 8 channels ($\bar p\bar n\to\bar d\gamma,\ \bar d\pi^0,\ \bar d\pi^+\pi^-,\ \bar d\pi^0\pi^0$; $\bar p\bar p\to\bar d\pi^-,\ \bar d\pi^-\pi^0$; $\bar n\bar n\to\bar d\pi^+,\ \bar d\pi^+\pi^0$); Sec. III C: "We set the cross sections to zero below the kinematic thresholds for the processes, as this is not ensured by the fit"; Sec. III D: "As in the $\bar N_1\bar N_2\to\bar d\pi$ case, we set the cross sections to zero below the kinematic thresholds". Channel list and charge conjugation: [`Deuteron_Production.C:26-67`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L26) (`Channels`, `Products`). — ✅

### 2.2 Formation probability (Eqs. (3), (4))

[`Deuteron_Production.C:141-148`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L141)
```cpp
  if (m_model==CoalescenceModel::CrossSectionModel) {
    for (const auto & ch : Channels()) {
      if (!Matches(ch,A,B)) continue;
      double sigma = Sigma(ch,A,B);
      if (sigma<=0.) continue;
      // sigma in microbarn, 1/sigma_0 in 1/barn.
      probs.push_back({ch, sigma, std::min(1.,sigma*1.e-6*m_inv_sigma0)});
```
D&R Eq. (4): $P(\bar N_1\bar N_2\to\bar dX_i\,|\,k)=\sigma_{\bar N_1\bar N_2\to\bar dX_i}(k)/\sigma_0$, with $\sigma_0$ "the same for all processes"; Sec. III A: $\sigma_{\bar N\bar N\to\bar dX}=\sigma_{NN\to d\bar X}$ is assumed. $1/\sigma_0$ is quoted in barn$^{-1}$ (Tables VII, VIII), hence the factor $10^{-6}$. — ✅

### 2.3 Choice of channel (Sec. III F)

[`Deuteron_Production.C:173-182`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L173)
```cpp
  for (const auto & cp : Probabilities(A,B,b2_fm2)) {
    if (ran->Get()<cp.prob) { success.push_back(cp); sigsum += cp.sigma; }
  }
  if (success.empty())   return Deuteron_Channel::none;
  if (success.size()==1 || sigsum<=0.) return success[0].channel;
  double r = ran->Get()*sigsum;
  for (const auto & cp : success) {
    if ((r -= cp.sigma)<=0.) return cp.channel;
  }
```
D&R Sec. III F, 3rd bullet: "Draw a random number $r_i$ uniformly on the unit interval for each possible formation process. If $r_i<P$ for one of the processes, the antinucleon pair forms an antideuteron through that process. If $r_i<P$ for more than one process, pick one of the processes randomly using probabilities equal to the relative cross sections." — ✅

### 2.4 Two-body final states $d\gamma$, $d\pi$ (Secs. III B, III C, III F)

[`Deuteron_Production.C:191-216`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L191)
```cpp
  // isotropic two-body decay in the rest frame of P, boosted back with P
  double M2 = P.Abs2(), M = sqrt(M2);
  double lambda = (M2-sqr(m1+m2))*(M2-sqr(m1-m2));
  double p = sqrt(std::max(0.,lambda))/(2.*M);
  Vec3D dir;
  Isotropic(p,dir);
  p1 = Vec4D(sqrt(p*p+m1*m1), dir);
  p2 = Vec4D(sqrt(p*p+m2*m2),-dir);
  Poincare boost(P);
  boost.BoostBack(p1);
  boost.BoostBack(p2);
```
D&R Sec. III B: "we let the antideuteron and photon be emitted back-to-back in the COM system, in a random direction drawn from an isotropic distribution"; Sec. III C: the same for $d\pi$; Sec. III F: "Boost the antideuteron to the lab frame." Momentum from four-momentum conservation (standard two-body kinematics). — ✅

### 2.5 Three-body final state $d\pi\pi$ (Sec. III D, Eq. (15))

[`Deuteron_Production.C:229-245`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L229)
```cpp
  for (size_t trial=0;trial<10000;trial++) {
    double m23 = m23min+ran->Get()*(m23max-m23min);
    double m12 = m12min+ran->Get()*(m12max-m12min);
    double M23 = sqrt(m23);
    double E2  = (m23-m3*m3+m2*m2)/(2.*M23), E1 = (s-m23-m1*m1)/(2.*M23);
    double q2  = sqrt(std::max(0.,E2*E2-m2*m2)), q1 = sqrt(std::max(0.,E1*E1-m1*m1));
    double lo  = sqr(E1+E2)-sqr(q1+q2), hi = sqr(E1+E2)-sqr(q1-q2);
    if (m12>=lo && m12<=hi) { m23sq = m23; break; }
  }
  ...
  double E_d = (s+m1*m1-m23sq)/(2.*sqrtS);
  double p_d = sqrt(std::max(0.,E_d*E_d-m1*m1));
  ...
  TwoBody(ppipi_cm,m2,m3,p2_cm,p3_cm);
```
(particle 1 = deuteron, 2 and 3 = pions; `m23` = $m^2_{\pi\pi}$, `m12` = $m^2_{d\pi}$)

D&R Sec. III D: "We determine the deuteron momentum by first drawing random invariant masses $m^2_{\pi\pi}$ and $m^2_{d\pi}$ uniformly within the kinematically allowed region. The momentum of the deuteron in the COM frame is then given by
$p_d=\sqrt{\left(\frac{s+m_d^2-m^2_{d\pi}}{2\sqrt s}\right)^2-m_d^2}$ (15), and we draw its direction from an isotropic distribution."

— ⚠️ **The code uses $m^2_{\pi\pi}$ where Eq. (15) is printed with $m^2_{d\pi}$.** The energy of particle 1 in a three-body state is $E_1=(s+m_1^2-m_{23}^2)/(2\sqrt s)$, with the invariant mass of the *other two* particles (Particle Data Group, *Review of Particle Physics*, "Kinematics" review, three-body decays / Dalitz plot); with $m^2_{d\pi}$ the expression is not the deuteron energy. We therefore read Eq. (15) as a typo. The Dalitz-region bounds for $m^2_{d\pi}$ at fixed $m^2_{\pi\pi}$ (variables `lo`, `hi`) are the standard ones from the same review. For a fixed $m_{\pi\pi}$ a uniform Dalitz density is equivalent to an isotropic $\pi\pi$ decay in the $\pi\pi$ rest frame, which is how the pion momenta are generated (`TwoBody`). The sampled $|p_d^*|$ is checked against the analytic phase-space distribution in `Tests()` (notebook Sec. 1). **Reviewer: please confirm the typo reading.**

### 2.6 The deuteron as a particle

[`SHERPA/SoftPhysics/Hadron_Init.C:53`](../../../SHERPA/SoftPhysics/Hadron_Init.C#L53)
```cpp
  InitHadron(kf_deuterium,   true, 1.875613, 2.1 ,0.,  3,2,true,true,1,"deuterium", "d");
```
Mass 1.875613 GeV (nuclear mass), charge $+1$ (`3` in units of $e/3$), spin 1 (`2` in units of 1/2), has an antiparticle (`true`), stable. Value of the mass: PDG/CODATA deuteron mass $m_d=1875.61\,$MeV. The same value is used in [`BaryonBaryon.H`](../../../HADRON_RESCATTERING/XSecs/BaryonBaryon.H) (`mDeuteron`). — ✅ (before this work the deuteron was registered as its own antiparticle, with spin 0 and the atomic mass)

---

## 3. Coalescence models — [`Deuteron_Production.C:151-166`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L151)

### 3.1 Momentum difference $k$

[`Deuteron_Production.C:99-104`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L99)
```cpp
double Deuteron_Coalescence::K(const Vec4D & pA,const Vec4D & pB) {
  Poincare boost(pA+pB);
  Vec4D pAcm = pA;
  boost.Boost(pAcm);
  return 2.*Vec3D(pAcm).Abs();
}
```
D&R Eq. (1) and the text after it ($k=|\vec p_p-\vec p_n|_{\rm COM}$). — ✅

### 3.2 `HardDisk` (classic coalescence)

```cpp
  // Coalescence-type models: only p n -> d gamma conserves charge.
  if (!Matches(Deuteron_Channel::pn_dgamma,A,B)) return probs;
  double k = K(A->Momentum(),B->Momentum()), P = 0.;
  bool   space = m_b_coal>0.;
  if (m_model==CoalescenceModel::HardDisk) {
    P = ((!space || b2_fm2<m_b_coal*m_b_coal) && k<m_p_coal) ? 1. : 0.;
  }
```
D&R Eq. (1): $P(\bar p\bar n\to\bar d\,|\,k)=\theta(p_0-k)$ — the coalescence model of S. T. Butler and C. A. Pearson, Phys. Rev. 129 (1963) 836, [doi:10.1103/PhysRev.129.836](https://doi.org/10.1103/PhysRev.129.836), and A. Schwarzschild and Č. Zupančič, Phys. Rev. 129 (1963) 854, [doi:10.1103/PhysRev.129.854](https://doi.org/10.1103/PhysRev.129.854). Final state: D&R Sec. II describes the radiative-capture view $\bar p\bar n\to\bar d\gamma$ "for a full kinematical description"; the code uses the $d\gamma$ kinematics of 2.4.
— ✅ for `B_COAL: 0` (momentum cut only). With `B_COAL > 0` the extra cut $b<b_{\rm coal}$ is ⓘ, see 3.4.

### 3.3 `Exponential` model

```cpp
  else {
    double C_b = (m_b_coal*m_b_coal)/(0.75*M_PI);
    double C_p = (m_p_coal*m_p_coal)/(0.75*M_PI);
    P = 0.75*(space ? std::exp(-b2_fm2/C_b) : 1.)*std::exp(-(k*k)/C_p);
  }
```
— ⓘ **Toy model** (our own, no published source). Closest published form: the Wigner function of a Gaussian deuteron wave function, M. Kachelrieß, S. Ostapchenko, J. Tjemsland, Eur. Phys. J. A 56 (2020) 4, [arXiv:1905.01192](https://arxiv.org/abs/1905.01192), **Eq. (11)**: $D(r,q)=8\,e^{-r^2/d^2}e^{-q^2d^2}$, with $d=3.2$ fm (reproduces the deuteron rms radius 1.96 fm). $D$ is not yet a probability: with the spin–isospin factor $S=3/8$ of Eq. (8), a $pn$ pair at separation $r$ and momentum difference $2q=k$ forms a deuteron with
$P=S\,D=3\,e^{-r^2/d^2}e^{-q^2d^2}$ (the same per-pair weight is Eq. (14) of the follow-up paper, Eur. Phys. J. A 57 (2021) 167, [arXiv:2012.04352](https://arxiv.org/abs/2012.04352)).

| | toy `Exponential` | from Eq. (11) |
|---|---|---|
| prefactor | 0.75 | 3 |
| distance variable | impact parameter $b$ | separation $r$ |
| spatial factor | $e^{-b^2/C_b}$, $C_b=b_{\rm coal}^2/(0.75\pi)$ | $e^{-r^2/d^2}$ |
| momentum factor | $e^{-k^2/C_p}$, $C_p=p_{\rm coal}^2/(0.75\pi)$ | $e^{-q^2d^2}=e^{-k^2d^2/(4(\hbar c)^2)}$ |

The toy model has the shape of Eq. (11) and reduces to it for prefactor 3, $b\to r$, $b_{\rm coal}=d\sqrt{0.75\pi}\simeq4.9$ fm and $p_{\rm coal}=2\hbar c\sqrt{0.75\pi}/d\simeq0.19$ GeV. The original default $p_{\rm coal}=\hbar c/3\,$fm was motivated in a code comment by the uncertainty principle.

### 3.4 Impact parameter $b$

[`Deuteron_Production.C:106-116`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L106)
```cpp
  // positions are in mm; in the pair frame A and B move along -+n, so the
  // separation transverse to n does not depend on the production times
  Poincare boost(A->Momentum()+B->Momentum());
  Vec4D xA = A->Position(), xB = B->Position(), pA = A->Momentum();
  boost.Boost(xA); boost.Boost(xB); boost.Boost(pA);
  Vec3D dx = Vec3D(xA)-Vec3D(xB), n = Vec3D(pA);
  if (n.Abs()>0.) n = n/n.Abs();
  double dl = dx*n;
  return std::max(0.,dx.Sqr()-dl*dl)*1.e24;   // mm^2 -> fm^2
```
— ⓘ Our definition: the distance of closest approach of the two straight trajectories in the pair frame, from Sherpa's hadron production positions (the same geometry as `Hadron_Collision::B2()` used by the rescattering, which is in mm²; the conversion to fm² is new). The literature's spatial coalescence conditions use the **separation** $\Delta r$ instead, e.g. $\Delta r<2$ fm in N. Fornengo, L. Maccione, A. Vittino, JCAP 09 (2013) 031, [arXiv:1306.4171](https://arxiv.org/abs/1306.4171) (D&R ref. [6]). **Caveat:** Sherpa's hadron positions in pp spread nucleons over tens to hundreds of fm (notebook Sec. 8), so any spatial condition is currently dominated by that.

---

## 4. The event phase — [`SHERPA/Single_Events/Deuteron_Formation.C`](../../../SHERPA/Single_Events/Deuteron_Formation.C)

Registered after the hadron decays in [`SHERPA/Main/Sherpa.C:204`](../../../SHERPA/Main/Sherpa.C#L204).

### 4.1 Nucleons from weak decays are excluded (Sec. IV)

[`Deuteron_Formation.C:60`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L60), [`:117-125`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L117)
```cpp
  m_widthcut = 0.1973269804/100.;
...
  while (blob && blob->Type()==btp::Hadron_Decay && blob->NInP()>0) {
    Particle * mother = blob->InParticle(0);
    if (mother->Flav().Width()<m_widthcut) return true;
    blob = mother->ProductionBlob();
  }
```
D&R Sec. IV: "antinucleons from weak decays should be excluded …, and we thus set all particles with mean lifetimes above 100 fm/c to be stable." A lifetime above 100 fm/c is a width below $\hbar c/(100\,{\rm fm})$; nucleons whose decay chain contains such a particle are dropped. — ✅ (implemented as exclusion instead of switching the decays off; the result for the nucleon sample is the same)

### 4.2 All unique pairs (Sec. III F, 1st bullet)

[`Deuteron_Formation.C:202-213`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L202)
```cpp
  for (size_t i=0;i<nucleons.size();i++) {
    for (size_t j=i+1;j<nucleons.size();j++) {
      ...
      // k = |pA - pB| in the pair c.m. frame, as in the paper
      double k = 2.*Vec3D(pAcm).Abs();
```
D&R Sec. III F: "For each event, iterate over all possible unique antinucleon–antinucleon pairs, avoiding double counting $\bar p\bar p$ and $\bar n\bar n$-pairs." — ✅ ($j>i$)

### 4.3 Used nucleons are excluded (Sec. III F, 3rd bullet)

[`Deuteron_Formation.C:251`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L251)
```cpp
      used[i] = used[j] = formed = true;
```
D&R Sec. III F: "If an antideuteron is formed, exclude the involved antinucleons from being used in the formation of other antideuterons." — ✅ (applies to the deuterons actually put into the event)

### 4.4 Expected yields: weighted pairs (Sec. III G)

[`Deuteron_Formation.C:224-230`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L224)
```cpp
      for (const auto & cp : p_deuterons->Probabilities(A,B,b2)) {
        ...
        for (size_t n=0;n<m_nsamp;n++) {
          std::vector<Vec4D> moms;
          if (!p_deuterons->FinalState(cp.channel,pA,pB,moms)) continue;
          FillDeuteron("expected",anti,moms[0],cp.prob/double(m_nsamp));
```
D&R Sec. III G, Eqs. (16)–(18): most Monte Carlo events form no antideuteron, so each event $i$ is evaluated $N_{\rm samp}$ times with the procedure of Sec. III F (new random numbers for the formation decision and the emission direction each time). $N^{\bar d}_{b,i}$ is the number of antideuterons from these $N_{\rm samp}$ evaluations that fall into histogram bin $b$, and
$w_{b,i}=N^{\bar d}_{b,i}/N_{\rm samp}$ (Eq. (17)) is the contribution of event $i$ to bin $b$, i.e. the average number per evaluation; the bin content is $N^{\bar d}_b=\sum_i w_{b,i}$ (Eq. (16)), its error $\sigma_b=\sqrt{\sum_i w_{b,i}^2}$ (Eq. (18)).

What the code does: instead of deciding at random whether a pair forms a deuteron (probability $P_i$), it adds the average of that decision, $P_i$, for every pair and open channel. `FinalState` draws one random emission direction and returns the deuteron momentum `moms[0]`; `FillDeuteron` enters it with weight $P_i/N_{\rm samp}$; this is repeated `m_nsamp` $=N_{\rm samp}$ times (setting `NSAMP`, default 10), so each pair contributes $P_i$ in total, spread over $N_{\rm samp}$ directions. Errors are from the sum of squared weights (histogram type 1), as in Eq. (18).
— ⚠️ The code fills **every** pair with its probability $P_i$ and samples only the kinematics $N_{\rm samp}$ times. This is the $N_{\rm samp}\to\infty$ expectation of Eqs. (16)–(17) for the formation decision, without the exclusion of used nucleons, which D&R (footnote 9) estimate to matter in <0.1 % of events. The deuterons actually formed (4.3) are kept as separate `formed_*` histograms; the validation run with $1/\sigma_0=18$/b (notebook Sec. 3) shows that both agree.

### 4.5 ALEPH and OPAL windows

[`Deuteron_Formation.C:24,27`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L24), filled in `FillDeuteron` ([`:127`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L127))
```cpp
  const double aleph_pmin = 0.62, aleph_pmax = 1.03, aleph_cosmax = 0.95;
  const double opal_pmin = 0.35, opal_pmax = 1.1;
```
ALEPH, S. Schael et al., Phys. Lett. B 639 (2006) 192, [hep-ex/0604023](https://arxiv.org/abs/hep-ex/0604023), Eq. (2): $(5.9\pm1.8\pm0.5)\times10^{-6}$ antideuterons per hadronic Z decay for $0.62<p<1.03$ GeV/c and $|\cos\theta|<0.95$. OPAL range $0.35<p<1.1$ GeV: D&R Sec. IV E (OPAL, R. Akers et al., Z. Phys. C 67 (1995) 203, [doi:10.1007/BF01571281](https://doi.org/10.1007/BF01571281)). — ✅

### 4.6 ALEPH $B_2$ inputs (ALEPH Eqs. (4), (5))

[`Deuteron_Formation.C:33,192-196`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L33)
```cpp
  const double aleph_half_pmin = 0.5*aleph_pmin, aleph_half_pmax = 0.5*aleph_pmax;
...
      if (p>aleph_half_pmin && p<aleph_half_pmax && std::abs(mom[3])/p<aleph_cosmax)
        m_evcounts["nucleon_"+sp+"_ALEPH_halfwindow"] += 1.;
```
ALEPH Eqs. (4), (5) are ALEPH's **measured values**: the invariant $\bar d$ cross section averaged over the sensitive angular range, $(3.7\pm1.1\pm0.3)\times10^{-6}$ GeV$^{-2}c^3$ (Eq. (4)), and $B_2=(3.3\pm1.0\pm0.8)\times10^{-3}$ GeV$^2$ (Eq. (5)), the former divided by the square of the antiproton invariant cross section at half the momentum, using antiprotons "from direct production, excluding those from weak decays". The code does not use these numbers; it counts what Sherpa needs to compute its **own** values the same way: the $\bar d$ in the ALEPH window (4.5) and, here, the (anti)protons in the window at half the momentum with the same angular cut. The antiprotons are the nucleons of 4.1 (no weak-decay products). — ✅ (the invariant yields and $B_2$ are formed in the notebook, 5.5)

**`m_evcounts`** (used here and in 4.5): a per-event counter, a map from a name (e.g. `"nucleon_pbar_ALEPH_halfwindow"`) to a number, cleared at the start of each event ([`:156`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L156)). At the end of the event every entry is added to `m_counts` (sum) and `m_counts2` (sum of squares); `Finish` writes the per-event mean and its error, $\sqrt{(\langle x^2\rangle-\langle x\rangle^2)/N_{\rm ev}}$, to `summary.dat`.

### 4.7 ALICE invariant yields

[`Deuteron_Formation.C:136`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L136)
```cpp
    m_histos[pre+"invyield_midy"]->Insert(pT,weight/(2.*pi*pT));
```
The observable of ALICE, Phys. Rev. C 97 (2018) 024615, [arXiv:1709.08522](https://arxiv.org/abs/1709.08522) (Fig. 5; HEPData [ins1625294](https://www.hepdata.net/record/ins1625294) Tables 1–3, column `E*D3N/DP**3`) and of D&R Fig. 7: the invariant yield
$E\frac{d^3N}{dp^3}=\frac{1}{2\pi p_T}\frac{d^2N}{dy\,dp_T}$ per event, at $|y|<0.5$ (the identity follows from $d^3p=p_T\,dp_T\,d\varphi\,dp_z$ and $dp_z=E\,dy$, averaged over $\varphi$).
How the code builds it: every deuteron with $|y|<0.5$ is entered at its $p_T$ with weight $w/(2\pi p_T)$, using its own $p_T$; in `Finish` the histogram is divided by $N_{\rm ev}$ and by the bin width $\Delta p_T$; the rapidity window has width $\Delta y=1$. Each bin therefore holds $\frac{1}{N_{\rm ev}}\sum w/(2\pi p_T)/(\Delta p_T\,\Delta y)$. — ✅

### 4.8 ALICE multiplicity classes: V0M and INEL>0

[`Deuteron_Formation.C:163-173`](../../../SHERPA/Single_Events/Deuteron_Formation.C#L163)
```cpp
  for (Particle * part : bloblist->ExtractParticles(part_status::active,1)) {
    if (part->DecayBlob() || part->Flav().IntCharge()==0) continue;
    ...
    double eta = mom.Eta();
    if ((eta>=2.8 && eta<=5.1) || (eta>=-3.7 && eta<=-1.7)) nV0++;
    if (std::abs(eta)<0.5) nch_mid++;
    if (std::abs(eta)<1.)  inel0 = true;
  }
```
What the variables are:
* `nV0`: the number of charged final-state particles in the acceptance of ALICE's two V0 scintillators — Sherpa's stand-in for the V0M amplitude, by which ALICE sorts events into multiplicity classes. The notebook orders all events by `nV0` and splits them at ALICE's percentiles.
* `nch_mid`: the number of charged particles at $|\eta|<0.5$; summed over a class and divided by its number of events it is $\langle dN_{ch}/d\eta\rangle$ (the $\eta$ interval has width 1), which ALICE quotes per class (Table 1) and uses as the abscissa of $d/p$ and $B_2$.
* `inel0`: true if the event has at least one charged particle at $|\eta|<1$ — ALICE's INEL>0 selection; only these events are classified.
* The first line skips neutral and non-final-state particles; the elided lines read the momentum and skip particles without transverse momentum (no $\eta$).

ALICE, Phys. Lett. B 794 (2019) 50, [arXiv:1902.09290](https://arxiv.org/abs/1902.09290), Sec. 2.1: V0-A $2.8\le\eta\le5.1$, V0-C $-3.7\le\eta\le-1.7$; Sec. 2.2: "The V0M amplitude is linearly proportional to the total number of charged particles produced in the V0 detectors acceptance", INEL>0 = "at least one charged particle being produced in the pseudorapidity interval $|\eta|<1$"; Table 1: class boundaries 0–4.7–9.5–19–38–100 % and $\langle dN_{ch}/d\eta\rangle$.
— ⚠️ (a) Charged final-state particles are used as V0M (no detector response). (b) Sherpa's percentiles are computed in the notebook from its own INEL>0 events, which are non-diffractive only (Amisic); ALICE's percentiles are of the full INEL>0 cross section. (c) Sherpa keeps particles with $c\tau>10$ mm stable by default, close to ALICE's primary-particle definition ($c\tau>1$ cm).

---

## 5. Analysis in the notebook — [`XSecs/Plotting.ipynb`](XSecs/Plotting.ipynb)

### 5.1 Fit of $1/\sigma_0$ (Sec. IV)

Notebook Sec. 6 (ALICE cell)
```python
best = (val*mc1/err**2).sum()/(mc1**2/err**2).sum()
```
Least-squares fit of one normalisation; the yields are linear in $1/\sigma_0$ (D&R Sec. IV: "the antideuteron spectrum in the cross section based model scales as $\propto1/\sigma_0$"), with only experimental uncertainties ("we only use the experimental uncertainty in calculating $\chi^2$"). Data errors: statistical ⊕ systematic in quadrature. — ✅ (the quadrature treatment of the correlated systematic error is our choice)

### 5.2 LEP: ALEPH + OPAL (Eqs. (21), (22))

Notebook Sec. 4
```python
chi2 = lambda x: (A - a*x)**2/sA**2 + eps*Nev*b*x
x_comb  = max(0., (A - eps*Nev*b*sA**2/(2*a))/a)      # d chi2/dx = 0
```
D&R Eq. (21): $s=\epsilon N_{ev}n_{\bar d,{\rm MC}}$ with $\epsilon=0.234$, $N_{ev}=1.64\times10^6$; Eq. (22): $\chi^2_{\rm OPAL}=s$. ALEPH term with $\sigma_A=\sqrt{1.8^2+0.5^2}\times10^{-6}$. — ✅

### 5.3 Non-diffractive fraction $f_{\rm ND}$

Notebook Sec. 6, `f_nd`: $f_{\rm ND}=1-\sigma_{\rm SD}/\sigma_{\rm inel}-\sigma_{\rm DD}/\sigma_{\rm inel}$ from ALICE, Eur. Phys. J. C 73 (2013) 2456, [arXiv:1208.4968](https://arxiv.org/abs/1208.4968), HEPData [ins1181770](https://www.hepdata.net/record/ins1181770) Tables 1, 2 (0.68 at all three energies).
— ⚠️ D&R Eqs. (19), (20) and Table IX convert non-diffractive Monte Carlo to V0AND-*triggered* events, because their (preliminary) ALICE data were normalised that way. The final ALICE data used here are normalised to *inelastic* events, so the analogous factor is $f_{\rm ND}=\sigma_{\rm ND}/\sigma_{\rm inel}$; like D&R we assume diffractive events produce no deuterons. ALICE's SD definition is $M_X<200$ GeV and DD $\Delta\eta>3$, so $f_{\rm ND}$ is approximate.

### 5.4 Coalescence momentum from the normalisation

Notebook Sec. 8: `p0 = P_COAL*s**(1/3)`.
D&R Sec. IV A: "In the coalescence model, the spectrum scales roughly as $p_0^3$". — ⚠️ Approximation; an exact $p_0$ needs runs at several $p_0$ (the yield is not exactly $\propto p_0^3$, and `Exponential` scales differently in detail).

### 5.5 Coalescence parameter $B_2$

Notebook Sec. 9 (ALICE) and Sec. 10 (multiplicity classes)
```python
v = scale*dv/(fnd*pv**2)                       # ALICE 2018: d/(p^2), yields per inelastic event
B2 = inv_d/inv_p**2                             # ALEPH
def inv_window(n, m, pmin, pmax):               # ALEPH: E/(p^2 dOmega dp) * N, |cos(theta)|<0.95
    p0 = 0.5*(pmin+pmax); E = np.sqrt(p0**2+m**2)
    return n*E/(p0**2*dOmega*(pmax-pmin))
```
ALICE 2018, Eq. (1): $E_A\frac{d^3N_A}{dp_A^3}=B_A\left(E_p\frac{d^3N_p}{dp_p^3}\right)^Z\left(E_n\frac{d^3N_n}{dp_n^3}\right)^N$ at $\vec p_p=\vec p_n=\vec p_A/A$, with neutrons = protons; ALICE 2019 per class with $(d+\bar d)/2$ and $(p+\bar p)/2$ (HEPData ins1721729 Tables 6–10). ALEPH Eqs. (4), (5): invariant cross sections at the centre of the momentum range, averaged over the sensitive solid angle. — ✅ (the per-event normalisation of $B_2\propto N_d/N_p^2$ is handled through $f_{\rm ND}$ for ALICE 2018; within a multiplicity class it cancels)

### 5.6 Ratio $d/p$

Notebook Sec. 10: `2*s_mult*d_dndy/(p_dndy+pbar_dndy)`. ALICE 2019 Table 12 / Fig. 8: "$2d/(p+\bar p)$", $p_T$-integrated yields. — ✅ (Sherpa integrates over all $p_T$; ALICE extrapolates outside its measured range)

### 5.7 Continuum limits

Notebook Sec. 9. ARGUS, H. Albrecht et al., Phys. Lett. B 236 (1990) 102, [doi:10.1016/0370-2693(90)90602-3](https://doi.org/10.1016/0370-2693(90)90602-3), abstract: $N_{\bar d}/N_{q\bar q}<1.7\times10^{-5}$ (90 % C.L.). CLEO, D. M. Asner et al., Phys. Rev. D 75 (2007) 012009, [hep-ex/0612019](https://arxiv.org/abs/hep-ex/0612019), Eq. (10): $\sigma(e^+e^-\to\bar dX)<0.031$ pb at $\sqrt s=10.5$ GeV. The Sherpa cross section uses $\sigma_{\rm had}=\sigma(e^+e^-\to\bar dX)/[\sigma(\bar dX)/\sigma_{\rm had}]=9.63\,{\rm fb}/3.01\times10^{-6}$ from BaBar, Phys. Rev. D 89 (2014) 111102, [arXiv:1403.4409](https://arxiv.org/abs/1403.4409), HEPData [ins1286317](https://www.hepdata.net/record/ins1286317) Tables 4, 5. — ⚠️ Sherpa is run at 10.58 GeV, not at ARGUS's ~10 GeV or CLEO's 10.5 GeV.

---

## 6. Code changes without new physics

| File | Change | Why |
|---|---|---|
| [`HADRON_RESCATTERING/Main/Hadron_Rescatterings.C:71`](../../../HADRON_RESCATTERING/Main/Hadron_Rescatterings.C#L71) | incoming hadrons of a rescattering are set to `decayed` | before, the original and the rescattered hadrons both stayed in the final state (double counting, energy not conserved) |
| [`SHERPA/SoftPhysics/Hadron_Rescattering_Handler.C`](../../../SHERPA/SoftPhysics/Hadron_Rescattering_Handler.C) | `delete` of collision candidates before `erase` (3 places) | memory leak of ~0.6 GB/min per 7 TeV process |
| [`HADRON_RESCATTERING/XSecs/BaryonBaryon.C`](../../../HADRON_RESCATTERING/XSecs/BaryonBaryon.C#L194) `Tests()` | particles built with the flavours of each channel | the output `pnToDGamma.dat` had one non-zero bin out of 10 000 |
| [`SHERPA/Main/Sherpa.C:204`](../../../SHERPA/Main/Sherpa.C#L204), [`SHERPA/Single_Events/CMakeLists.txt`](../../../SHERPA/Single_Events/CMakeLists.txt) | the `Deuteron_Formation` phase is added and built | — |

## 7. How the code is checked

`DEUTERON_FORMATION: {TEST: true}` runs [`Deuteron_Coalescence::Tests`](../../../HADRON_RESCATTERING/Main/Deuteron_Production.C#L321) at start-up and prints:
- four-momentum conservation (≲10⁻¹⁵ relative), mass shells (≲10⁻¹³ GeV²) and charge for all 8 channels, nucleons and antinucleons, at rest and boosted;
- the channel selection of 2.3 against the exact probability of the paper's rule (pulls);
- impact parameter and coalescence probabilities against hand-calculated values;
- `test_sigma_*.dat`, `test_P_*.dat`, `test_pd_cm_*.dat` (vs. phase space), `test_costheta_d_*.dat` (isotropy), plotted in notebook Sec. 1.

## References

1. L. A. Dal, A. R. Raklev, Phys. Rev. D 91 (2015) 123536, [arXiv:1504.07242](https://arxiv.org/abs/1504.07242)
2. H. Machner, J. Niskanen, Nucl. Phys. A 776 (2006) 172, [nucl-ex/0511027](https://arxiv.org/abs/nucl-ex/0511027) (source of Eq. (10), via D&R)
3. S. T. Butler, C. A. Pearson, Phys. Rev. 129 (1963) 836, [doi:10.1103/PhysRev.129.836](https://doi.org/10.1103/PhysRev.129.836)
4. A. Schwarzschild, Č. Zupančič, Phys. Rev. 129 (1963) 854, [doi:10.1103/PhysRev.129.854](https://doi.org/10.1103/PhysRev.129.854)
5. N. Fornengo, L. Maccione, A. Vittino, JCAP 09 (2013) 031, [arXiv:1306.4171](https://arxiv.org/abs/1306.4171)
6. M. Kachelrieß, S. Ostapchenko, J. Tjemsland, Eur. Phys. J. A 56 (2020) 4, [arXiv:1905.01192](https://arxiv.org/abs/1905.01192); Eur. Phys. J. A 57 (2021) 167, [arXiv:2012.04352](https://arxiv.org/abs/2012.04352)
7. ALEPH, S. Schael et al., Phys. Lett. B 639 (2006) 192, [hep-ex/0604023](https://arxiv.org/abs/hep-ex/0604023)
8. OPAL, R. Akers et al., Z. Phys. C 67 (1995) 203, [doi:10.1007/BF01571281](https://doi.org/10.1007/BF01571281)
9. BaBar, J. P. Lees et al., Phys. Rev. D 89 (2014) 111102, [arXiv:1403.4409](https://arxiv.org/abs/1403.4409), HEPData [ins1286317](https://www.hepdata.net/record/ins1286317)
10. ARGUS, H. Albrecht et al., Phys. Lett. B 236 (1990) 102, [doi:10.1016/0370-2693(90)90602-3](https://doi.org/10.1016/0370-2693(90)90602-3)
11. CLEO, D. M. Asner et al., Phys. Rev. D 75 (2007) 012009, [hep-ex/0612019](https://arxiv.org/abs/hep-ex/0612019)
12. ALICE, Phys. Rev. C 97 (2018) 024615, [arXiv:1709.08522](https://arxiv.org/abs/1709.08522), HEPData [ins1625294](https://www.hepdata.net/record/ins1625294)
13. ALICE, Phys. Lett. B 794 (2019) 50, [arXiv:1902.09290](https://arxiv.org/abs/1902.09290), HEPData [ins1721729](https://www.hepdata.net/record/ins1721729)
14. ALICE, Eur. Phys. J. C 73 (2013) 2456, [arXiv:1208.4968](https://arxiv.org/abs/1208.4968), HEPData [ins1181770](https://www.hepdata.net/record/ins1181770)
15. ALICE, Eur. Phys. J. C 75 (2015) 226, [arXiv:1504.00024](https://arxiv.org/abs/1504.00024), HEPData [ins1357424](https://www.hepdata.net/record/ins1357424)
16. Particle Data Group, *Review of Particle Physics*, "Kinematics" review (three-body decays, Dalitz plot), [pdg.lbl.gov](https://pdg.lbl.gov)
