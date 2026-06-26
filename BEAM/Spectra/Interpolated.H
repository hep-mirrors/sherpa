#ifndef BEAM_Main_Interpolated_Neutrinos_H
#define BEAM_Main_Interpolated_Neutrinos_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Lookup_Tables.H"


namespace BEAM {
  class Interpolated_Neutrinos : public Beam_Base {
  private:
    ATOOLS::OneDim_Flexible_Table * p_table;
    void Initialise();
  public:
    Interpolated_Neutrinos(const ATOOLS::Flavour & flav,
			   const double & energy,const int dir);
    ~Interpolated_Neutrinos() {
      if (p_table) delete p_table;
    }
    Beam_Base     * Copy();
    ATOOLS::Flavour Remnant();
    bool CalculateWeight(const double,const double);
  };
}

#endif
