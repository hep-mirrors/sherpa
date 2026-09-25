/*!
  \file Ceex_Real.C

  The contractions that fold the emission matrices U and V into a Born.

  This file used to hold beta_1^0 as well, transcribed from KKMC's HiniPlus
  and HfinPlus. Those are gone: they are a formula for e+e- -> f fbar rather
  than a method - the photon's spinor index goes into a beam slot, which only
  means anything for a single s-channel 2 -> 2 - and they gave a CEEX cross
  section 3.1x the NLO for e+e- -> nu nu. beta_1 now comes from Comix's own
  one-photon amplitude with the eikonal subtracted, in
  Ceex_Base::ComixInfraredSubtracted_1_0. The contractions below are still
  used by beta_2 and by the Comix alignment.
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


void Ceex_Base::AddU(Complex &sum, const Amplitude &Born, const Amplitude &U, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex CSum(0, 0);
          for (int  j = 0; j <= 1; j++) {
            CSum += fac * Born.m_A[Idx(j,h1,h2,h3)] * U.m_U[j][h1];
          }
          sum += CSum;
        }
      }
    }
  }
}


void Ceex_Base::AddU(Amplitude &out, const Amplitude &Born,
                     const Amplitude &U, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0)
    for (int h1 = 0; h1 <= 1; ++h1)
      for (int h2 = 0; h2 <= 1; ++h2)
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex c(0., 0.);
          for (int j = 0; j <= 1; ++j)
            c += fac * Born.m_A[Idx(j,h1,h2,h3)] * U.m_U[j][h0];
          out.m_A[Idx(h0,h1,h2,h3)] += c;
        }
}



void Ceex_Base::AddV(Complex &sum, const Amplitude &Born, const Amplitude &V, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex CSum(0, 0);
          for (int  j = 0; j <= 1; j++) {
            CSum += fac * (V.m_V[h2][j]) * Born.m_A[Idx(h0,j,h2,h3)];
          }
          sum += CSum;
        }
      }
    }
  }
}


void Ceex_Base::AddV(Amplitude &out, const Amplitude &Born,
                     const Amplitude &V, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0)
    for (int h1 = 0; h1 <= 1; ++h1)
      for (int h2 = 0; h2 <= 1; ++h2)
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex c(0., 0.);
          for (int j = 0; j <= 1; ++j)
            c += fac * V.m_V[h1][j] * Born.m_A[Idx(h0,j,h2,h3)];
          out.m_A[Idx(h0,h1,h2,h3)] += c;
        }
}

void Ceex_Base::AddUF(Amplitude &out, const Amplitude &U,
                      const Amplitude &Born, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0)
    for (int h1 = 0; h1 <= 1; ++h1)
      for (int h2 = 0; h2 <= 1; ++h2)
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex c(0., 0.);
          for (int j = 0; j <= 1; ++j)
            c += fac * U.m_U[h2][j] * Born.m_A[Idx(h0,h1,j,h3)];
          out.m_A[Idx(h0,h1,h2,h3)] += c;
        }
}


void Ceex_Base::AddVF(Amplitude &out, const Amplitude &Born,
                      const Amplitude &V, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0)
    for (int h1 = 0; h1 <= 1; ++h1)
      for (int h2 = 0; h2 <= 1; ++h2)
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex c(0., 0.);
          for (int j = 0; j <= 1; ++j)
            c += fac * Born.m_A[Idx(h0,h1,h2,j)] * V.m_V[j][h3];
          out.m_A[Idx(h0,h1,h2,h3)] += c;
        }
}




void Ceex_Base::SumAmplitude(Complex &sum, const Amplitude &Amp, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          sum += fac * Amp.m_A[Idx(h0,h1,h2,h3)];
        }
      }
    }
  }
}



void Ceex_Base::SumAmplitude(Complex &sum, const Amplitude &Amp1, const Amplitude &Amp2, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          for (int  j = 0; j <= 1; j++) {
            sum += fac * Amp2.m_A[Idx(j,h1,h2,h3)] * Amp1.m_U[j][h1];
          }
        }
      }
    }
  }
}
