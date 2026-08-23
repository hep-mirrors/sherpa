#ifndef HADRONS_PS_Library_Five_Body_PSs_H
#define HADRONS_PS_Library_Five_Body_PSs_H

#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Math/Vector.H"
#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Vegas.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"
#include "PHASIC++/Channels/Rambo.H"

namespace HADRONS {
  // ------------------------------------------------------------------
  // ThreeResonances: 1->5 phase-space channel with a THREE-resonance
  // sequential cascade, direct extension of Four_Body_PSs.H's
  // TwoResonances (1->4, two resonances) by one more nesting level.
  // Needed for 4-pseudoscalar tau decay final states (4 hadrons + 1
  // neutrino = 5 outgoing particles), analogous to how TwoResonances
  // already serves the 3-pseudoscalar (a1->rho+pi -style) case.
  //
  // Decay topology (M = fl[0], the tau's hadronic system momentum):
  //   M -> dir + res1;  res1 -> k + res2;  res2 -> l + res3;
  //   res3 -> i + j.
  // "dir" (the outermost spectator) is deduced automatically, exactly
  // as TwoResonances deduces its own m_dir - it is whichever of the 5
  // outgoing indices is none of i,j,k,l.
  // ------------------------------------------------------------------
  class ThreeResonances : public PHASIC::Single_Channel {
    ATOOLS::Vec4D     m_P;
    int           	  m_i, m_j, m_k, m_l, m_dir;
    int           	  m_chnumber;
    SimpleResonanceFlavour  m_prop1, m_prop2, m_prop3;
    ATOOLS::Flavour	* p_fl;

    ATOOLS::Info_Key m_kI_1234_5,m_kI_123_4,m_kI_12_3,m_kI_1_2;
    PHASIC::Vegas* p_vegas;
    ATOOLS::Integration_Info * p_info;

  public :
    ThreeResonances(const ATOOLS::Flavour * fl,SimpleResonanceFlavour prop1,
		    const int _k,SimpleResonanceFlavour prop2,
		    const int _l,SimpleResonanceFlavour prop3,
		    const int _i,const int _j );
    ~ThreeResonances();
    void GeneratePoint(ATOOLS::Vec4D * p,PHASIC::Cut_Data * =NULL,
		       double * _ran=NULL);
    void GenerateWeight(ATOOLS::Vec4D * p,PHASIC::Cut_Data * =NULL);
    int  ChNumber()                 { return m_chnumber;       }
    void SetChNumber(int _chnumber) { m_chnumber = _chnumber;  }
	void   AddPoint(double value )	{
	  Single_Channel::AddPoint( value );
	  p_vegas->AddPoint(value, p_rans);
	}
    void   Optimize()                { p_vegas->Optimize(); }
    void   EndOptimize()             { p_vegas->EndOptimize(); }
    void   WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void   ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }
    void MPISync() {};
  }; // end of class

  /*!
    \file Five_Body_PSs.H
    \brief Declares the class HADRONS::ThreeResonances

    This file can be found in the directory <tt>PS_Library</tt>.
  */
  /*!
    \class ThreeResonances
    \brief Tool to handle 5-particle PS with 3 sequential resonances

    This class is a subclass of PHASIC::Single_Channel. It is the
    default integrator for decays of the form with three resonances
    (internal propagators)
    \f[ M \to l + res_1; \: res_1 \to k+res_2; \: res_2 \to m+res_3;
        \: res_3 \to i+j. \f]
    Direct extension of TwoResonances (Four_Body_PSs.H) by one more
    sequential resonance, for 1->5 (e.g. 4 pseudoscalars + neutrino)
    final states.
  */

} // end of namespace

#endif
