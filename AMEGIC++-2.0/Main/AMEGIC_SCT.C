#include "AMEGIC_SCT.H"
#include "Message.H"
#include "MyComplex.H"

#include "prof.hh"

#define DEBUG_AMEGIC_SCT

namespace AMEGIC {

  AMEGIC_SCT::AMEGIC_SCT(std::vector<int> A1, std::vector<int> A2, 
			 Helicity* Hel, std::vector<int>* particles, size_t pPos)
  {
    PROFILE_HERE;
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
      default: msg.Error()<<"Warning in AMEGIC_SCT::HADRONS_to_AMEGIC:"<<std::endl
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
    default: msg.Error()<<"Warning in AMEGIC_SCT::HADRONS_to_AMEGIC:"<<std::endl
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

    PROFILE_HERE;
    SCT_DUMMY* SCT = new SCT_DUMMY();
    if (m_particle==-1) {
      SCT->m_particle=-1;
      Complex val(0., 0.);
#ifdef DEBUG_AMEGIC_SCT
      if (m_A1.size() != m_A2.size()) {
	PRINT_INFO("Now you really have something to think about!");
	abort();
      }
#endif
      for (size_t ampl=0; ampl<m_A1.size(); ++ampl)
	for (size_t c1=0; c1<graphs->size(); ++c1)
	  for (size_t c2=0; c2<graphs->size(); ++c2)
	    val +=  (*graphs)[c1]->Zvalue(m_A1[ampl])
	           *conj((*graphs)[c2]->Zvalue(m_A2[ampl]))
	           *col->Mij(c1, c2)*hel->PolarizationFactor(m_A1[ampl]);
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










  /*   AMEGIC_SCT::AMEGIC_SCT(std::vector< std::vector<Complex> > *Ampls, 
			  CFColor *ColorMatrix, Helicity *Hel,
			  std::vector<int> comb1, std::vector<int> comb2)
   {
#ifdef DEBUG_AMEGIC_SCT
     if (comb1.size() != comb2.size()) {
       PRINT_INFO("Error: comb1 and comb2 do not have the same length!");
       abort();
     }
#endif
     if (comb1.size() == Hel->Nflavs()) {
       // If an end-node is reached, compute the entry
       m_particle = -1;

       // get the helicity numbers
       size_t i1=Hel->GetAmplitudeNumber(&comb1);
       size_t i2=Hel->GetAmplitudeNumber(&comb2);

       // Get the vectors for the amplitudes of the different graphs.
       std::vector<Complex> *v1 = &(*Ampls)[i1];
       std::vector<Complex> *v2 = &(*Ampls)[i2];

#ifdef DEBUG_AMEGIC_SCT
       if (v1->size() != v2->size()) {
	 PRINT_INFO("Error in Constructor of AMEGIC_SCT:"<<std::endl<<
		    "The amplitudes belonging to the "<<"helicities "<<i1<<" and "<<i2
		    <<" have a different number of associated graphs:");
	 PRINT_INFO(v1->size()<<" vs. "<<v2->size());
       }
#endif

       // Create and store the value M M*
       for (size_t i=0; i<v1->size(); ++i)
	 for (size_t j=0; j<v2->size(); ++j)
	   m_value += (*v1)[i] * conj((*v2)[j]) * ColorMatrix->Mij(i,j);
     
     } else {

       // Not an end-node, yet. Add all the possible helicity combinations for the current
       // particle to the end of the comb1/2 list and go on with recursive creation.

       size_t pos = m_particle = comb1.size();
       size_t numHel = Hel->MaxHel(pos);
    
       p_next = new std::vector<Spin_Correlation_Tensor*>(numHel*numHel);
       
       // extend length of comb1 and comb2
       comb1.push_back(0);
       comb2.push_back(0);
       
       for (size_t i=0; i<numHel; ++i) {
	 comb1.back() = HADRONS_to_AMEGIC(i, numHel);
	 for (size_t j=0; j<numHel; ++j) {
	   comb2.back() = HADRONS_to_AMEGIC(j, numHel);
	   // *** either numHel*i+j or numHel*j+i *** check!
	   (*p_next)[numHel*i+j] 
	     = new AMEGIC_SCT(Ampls, ColorMatrix, Hel, comb1, comb2);
	 }
       }
     }
     } */


} // end of namespace AMEGIC
