#include "AMEGIC_SCT.H"
#include "Message.H"

#define DEBUG_AMEGIC_SCT

namespace AMEGIC {

   AMEGIC_SCT::AMEGIC_SCT(std::vector< std::vector<Complex> > *Ampls, 
			  CFColor *ColorMatrix, Helicity *Hel,
			  std::vector<int> comb1, std::vector<int> comb2)
   {
     //     PRINT_INFO("Creating...");
     //     PRINT_INFO("sizeof(Ampls) = "<<Ampls->size());
     //PRINT_INFO("comb1.size() = "<<comb1.size()); 
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
   };

  AMEGIC_SCT::~AMEGIC_SCT() {};


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
}
