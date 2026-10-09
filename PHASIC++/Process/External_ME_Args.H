#ifndef PHASIC_Process_External_ME_Args_H
#define PHASIC_Process_External_ME_Args_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Scoped_Settings.H"

namespace PHASIC {

  /* This struct is to be passed to the getters for external MEs and
     correlated MEs. Specifies flavours, couplings and source (i.e
     OpenLoops, MadGraph, etc) */

  struct External_ME_Args {
    
    External_ME_Args(const ATOOLS::Flavour_Vector& inflavs,
		     const ATOOLS::Flavour_Vector& outflavs,
		     const std::vector<double>& orders,
		     const std::string& source = "",
		     const ATOOLS::Scoped_Settings& currents
		       = ATOOLS::Scoped_Settings(),
		     bool hascurrents = false)
      : m_inflavs(inflavs), m_outflavs(outflavs), m_orders(orders),
        m_source(source), m_currents(currents),
        m_hascurrents(hascurrents) {}
      
    const ATOOLS::Flavour_Vector m_inflavs;
    const ATOOLS::Flavour_Vector m_outflavs;
    const std::vector<double>    m_orders;
    const std::string            m_source;
    // The process' Currents: block. Only meaningful when m_hascurrents is
    // set: a default-constructed Scoped_Settings resolves to the settings
    // ROOT, so "is it empty" is not a usable test here.
    ATOOLS::Scoped_Settings      m_currents;
    bool                         m_hascurrents;

    ATOOLS::Flavour_Vector Flavours() const
    {
      ATOOLS::Flavour_Vector flavs = m_inflavs;
      flavs.insert(flavs.end(),
		   m_outflavs.begin(),
		   m_outflavs.end());
      return flavs;
    }
  };

}

#endif
