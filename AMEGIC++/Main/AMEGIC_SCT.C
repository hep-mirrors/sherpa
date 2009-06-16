#include "AMEGIC++/Main/AMEGIC_SCT.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MyComplex.H"

namespace AMEGIC {

  AMEGIC_SCT::AMEGIC_SCT(std::vector<int> A1, std::vector<int> A2, 
			 Helicity* Hel, std::vector<int>* particles, size_t pPos)
  {
    if (pPos == particles->size()) {
      m_particle = -1;
      m_A1 = A1;
      m_A2 = A2;
    }
    else if (pPos < particles->size()) {
      m_particle = (*particles)[pPos];
      size_t max(Hel->MaxHel(m_particle));

      std::vector<std::vector<int> > B1(max), B2(max);
      for (size_t i=0; i<A1.size(); ++i) {
	int pol(Hel->GetPol(m_particle, A1[i]));
	pol = AMEGIC_to_HADRONS(pol, max);
	  B1[pol].push_back(A1[i]);
      }
      for (size_t i=0; i<A2.size(); ++i) {
	int pol(Hel->GetPol(m_particle, A2[i]));
	pol = AMEGIC_to_HADRONS(pol, max);
	B2[pol].push_back(A2[i]);
      }

      if (m_mode == scmode::Full) {
	p_next = new std::vector<Spin_Correlation_Tensor*> (max*max);	
	for (size_t i=0; i<max; ++i) 
	  for (size_t j=0; j<max; ++j) 
	    (*p_next)[max*i+j] = new AMEGIC_SCT(B1[i], B2[j], Hel, particles, pPos+1);
      }
      if (m_mode == scmode::Diagonal) {
	p_next = new std::vector<Spin_Correlation_Tensor*> (max);
	for (size_t i=0; i<max; ++i)
	  (*p_next)[i] = new AMEGIC_SCT(B1[i], B2[i], Hel, particles, pPos+1);
      }
    }
  }



  AMEGIC_SCT::~AMEGIC_SCT() {}
    

  int AMEGIC_SCT::HADRONS_to_AMEGIC(size_t pol, size_t maxPol) 
  {
    switch (maxPol) {
      case 1: return 0;
      case 2: switch(pol)
        { case 0: return 1;
          case 1: return 0; }
      case 3: switch(pol) 
	{ case 0: return 1;
	  case 1: return 2;
	  case 2: return 0; }
      default: msg_Error()<<"Warning in AMEGIC_SCT::HADRONS_to_AMEGIC:"<<std::endl
	                  <<"Particles with a maximum number of "<<maxPol
			  <<" are not supported, yet! Abort the run."<<std::endl;
 	  return 0;
    }
    return 0;
  }

  size_t AMEGIC_SCT::AMEGIC_to_HADRONS(int pol, size_t maxPol) 
  {
    switch (maxPol) {
    case 1: return 0;
    case 2: switch(pol) {
    case -1 : return 1;
    case 1: return 0;
    }
    case 3: switch(pol) {
    case -1 : return 2;
    case 0: return 1;
    case 1: return 0;
    }
    default: msg_Error()<<"Warning in AMEGIC_SCT::HADRONS_to_AMEGIC:"<<std::endl
      <<"Particles with a maximum number of "<<maxPol
      <<" are not supported, yet! Abort the run."<<std::endl;
      return 0;
    }
    return 0;
  }

  Spin_Correlation_Tensor* AMEGIC_SCT::CreateSCT(std::vector<Amplitude_Base* > *graphs,
					      CFColor *col, Helicity* hel)
  {
    if (m_mode == scmode::None) return NULL;

    SCT_DUMMY* SCT = new SCT_DUMMY();
    if (m_particle==-1) {
      SCT->m_particle=-1;
      Complex val(0., 0.);
      for (size_t ampl=0; ampl<m_A1.size(); ++ampl)
	for (size_t c1=0; c1<graphs->size(); ++c1)
	  for (size_t c2=0; c2<graphs->size(); ++c2) {
	    val +=  (*graphs)[c1]->Zvalue(m_A1[ampl])
	           *conj((*graphs)[c2]->Zvalue(m_A2[ampl]))
	           *col->Mij(c1, c2)*hel->PolarizationFactor(m_A1[ampl]);
	    if (hel->Multiplicity(ampl) != 1) {
	      msg_Error()<<"Encoutered process for which the current version of AMEGIC cannot "
			 <<"create spin informations."<<std::endl
			 <<"Please restart the run with option SPIN_CORRELATIONS=0 or "
			 <<"contact the authors at support@sherpa-mc.de for a solution."<<std::endl;
	      abort();
	    }
	  }
      SCT->m_value=val;
    }
    else {
      SCT->m_particle = m_particle;
      SCT->p_next = new std::vector<Spin_Correlation_Tensor*>(p_next->size());
      for (size_t i=0; i<p_next->size(); ++i)
	(*SCT->p_next)[i] = ((AMEGIC_SCT*)(*p_next)[i])->CreateSCT(graphs, col, hel);
    }
    return SCT;
  }
} // end of namespace AMEGIC
