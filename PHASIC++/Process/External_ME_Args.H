#ifndef PHASIC_Process_External_ME_Args_H
#define PHASIC_Process_External_ME_Args_H

#include "ATOOLS/Phys/Flavour.H"

namespace PHASIC {

  /* This struct is to be passed to the getters for external MEs and
     correlated MEs. Specifies flavours, couplings and source (i.e
     OpenLoops, MadGraph, etc) */

  struct External_ME_Args {
    
    External_ME_Args(const ATOOLS::Flavour_Vector& inflavs,
		     const ATOOLS::Flavour_Vector& outflavs,
		     const std::vector<double>& orders,
		     const std::string& source = "")
      : m_inflavs(inflavs), m_outflavs(outflavs), m_orders(orders), m_source(source) {}
      
    const ATOOLS::Flavour_Vector m_inflavs;
    const ATOOLS::Flavour_Vector m_outflavs;
    const std::vector<double>    m_orders;
    const std::string            m_source;
  };

}

#endif
