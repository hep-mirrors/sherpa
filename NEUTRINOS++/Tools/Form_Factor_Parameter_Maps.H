#ifndef NEUTRINOS_Tools_Form_Factor_Parameter_Maps_H
#define NEUTRINOS_Tools_Form_Factor_Parameter_Maps_H

#include "NEUTRINOS++/Tools/Form_Factor_Library.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <set>

namespace ATOOLS {
  class Scoped_Settings;
}

namespace NEUTRINOS {
  struct cpl_info {
    enum code {
      unknown      = 0,
      scalar       = 1,
      pseudoscalar = 2,
      vector       = 3,
      axialvector  = 4,
      tensor       = 5
    };
  };
  std::ostream & operator<<(std::ostream & s,const cpl_info::code & cpl);
  
  struct ff_info {
    cpl_info::code      m_cpl;
    ff_type::code       m_type;
    std::vector<double> m_params;
    ff_info(const cpl_info::code & cpl, const ff_type::code & type, 
	    const size_t & nparams=2) :
      m_cpl(cpl), m_type(type) { m_params.resize(nparams); }
  };
  std::ostream & operator<<(std::ostream & s,const ff_info & info);
  
  
  struct ModelParameters: public std::map<std::string, double> {
  };
  
  class Form_Factor_Entry {
  private:
    std::pair<kf_code, kf_code > m_pair;
    std::map<kf_code, std::list<ff_info * > > m_entries;
  public:
    Form_Factor_Entry(std::pair<kf_code, kf_code > pair) : m_pair(pair) {}
    ~Form_Factor_Entry();
    void Add(const kf_code & prop, std::list<ff_info * > & ffinfos);
    std::pair<kf_code, kf_code > GetNucleonPair() const { return m_pair; }
    std::map<kf_code, std::list<ff_info * > > * GetAllEntries() {
      return &m_entries;
    }
    std::list<ff_info * > * GetEntries(kf_code & kfc) {
      if (m_entries.find(kfc)==m_entries.end()) return NULL;
      return &m_entries[kfc];
    }
  };
  std::ostream & operator<<(std::ostream & s,Form_Factor_Entry & entry);

  class Form_Factor_Parameter_Maps :
    public std::map<std::pair<kf_code, kf_code > , Form_Factor_Entry * > {
  private:
    ATOOLS::Scoped_Settings m_alltransitions;
    ModelParameters         m_parameters;
    
    void InitialiseMaps();
    void InitialiseConstants();
    void InitialiseFFMaps();
    
    bool InitialiseFormFactor(std::string & name);
    bool ExtractFlavours(std::string & buffer,std::vector<kf_code> & kfcs);
    bool ExtractParameters(std::string & buffer,std::vector<double> & parameters);
    inline cpl_info::code CplInfo(const std::string & cpl) {
      if (cpl=="scalar")       return cpl_info::scalar;
      if (cpl=="pseudoscalar") return cpl_info::pseudoscalar;
      if (cpl=="vector")       return cpl_info::vector;
      if (cpl=="axialvector")  return cpl_info::axialvector;
      if (cpl=="tensor")       return cpl_info::tensor;
      return cpl_info::unknown;
    }
    inline ff_type::code FFType(const std::string & form) {
      if (form=="none")        return ff_type::none;
      if (form=="dipole")      return ff_type::dipole;
      if (form=="exponential") return ff_type::exponential;
      if (form=="Gaussian")    return ff_type::Gaussian;
      return ff_type::unknown;
    } 
    inline Form_Factor_Entry * Entry(kf_code & n1,kf_code & n2) {
      std::pair<kf_code,kf_code> kfpair = std::make_pair(n1,n2);
      if (find(kfpair)==end()) (*this)[kfpair] = new Form_Factor_Entry(kfpair);
      return (*this)[kfpair];
    }
  public:
    Form_Factor_Parameter_Maps();
    ~Form_Factor_Parameter_Maps();
    void Output();
  };
}
#endif
