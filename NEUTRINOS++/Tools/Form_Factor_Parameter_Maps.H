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
  class Form_Factor_Entry {
  private:
    std::pair<kf_code, kf_code > m_pair;
    std::map<kf_code, std::map<cpl_info::code, ff_info * > > m_entries;
  public:
    Form_Factor_Entry(std::pair<kf_code, kf_code > pair) : m_pair(pair) {}
    ~Form_Factor_Entry();
    void Add(const kf_code & prop, std::map<cpl_info::code, ff_info * >  & ffinfos);
    std::pair<kf_code, kf_code > GetNucleonPair() const { return m_pair; }
    std::map<kf_code, std::map<cpl_info::code, ff_info * > > * GetAll() { return &m_entries; }
    std::map<cpl_info::code, ff_info * >  * Get(kf_code & kfc) {
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
    
    void InitializeConstants();
    void InitializeFFMaps();
    
    bool InitializeFormFactor(std::string & name);
    bool ExtractFlavours(std::string & buffer,std::vector<kf_code> & kfcs);
    inline cpl_info::code CplInfo(const std::string & cpl) {
      if (cpl=="scalar")       return cpl_info::scalar;
      if (cpl=="pseudoscalar") return cpl_info::pseudoscalar;
      if (cpl=="vector")       return cpl_info::vector;
      if (cpl=="axialvector")  return cpl_info::axialvector;
      if (cpl=="tensor")       return cpl_info::tensor;
      if (cpl=="G_E")          return cpl_info::GE;
      if (cpl=="G_M")          return cpl_info::GM;
      if (cpl=="f1")           return cpl_info::f1;
      if (cpl=="f2")           return cpl_info::f2;
      if (cpl=="f3")           return cpl_info::f3;
      if (cpl=="g1")           return cpl_info::g1;
      if (cpl=="g2")           return cpl_info::g2;
      if (cpl=="g3")           return cpl_info::g3;
      return cpl_info::unknown;
    }
    inline ff_type::code FFType(const std::string & form) {
      if (form=="none")             return ff_type::none;
      if (form=="constant")         return ff_type::constant;
      if (form=="dipole")           return ff_type::dipole;
      if (form=="lambda_dipole")    return ff_type::lambda_dipole;
      if (form=="neutron_electric") return ff_type::neutron_electric;
      if (form=="exponential")      return ff_type::exponential;
      if (form=="Gaussian")         return ff_type::Gaussian;
      if (form=="Kelly")            return ff_type::Kelly;
      if (form=="BBBA")             return ff_type::BBBA;
      if (form=="ArringtonHill")    return ff_type::ArringtonHill;
      if (form=="Helm")             return ff_type::Helm;
      if (form=="Lovato")           return ff_type::Lovato;
      return ff_type::unknown;
    } 
    inline Form_Factor_Entry * Entry(kf_code & n1,kf_code & n2) {
      std::pair<kf_code,kf_code> kfpair = std::make_pair(n1,n2);
      if (find(kfpair)==end()) (*this)[kfpair] = new Form_Factor_Entry(kfpair);
      return (*this)[kfpair];
    }
    inline Form_Factor_Entry * FindEntry(kf_code & n1,kf_code & n2) {
      std::pair<kf_code,kf_code> kfpair = std::make_pair(n1,n2);
      if (find(kfpair)==end()) return NULL;
      return (*this)[kfpair];
    }
    inline ff_info * FindInfo(kf_code & n1,kf_code & n2,kf_code & prop,cpl_info::code & cpl) {
      Form_Factor_Entry * entry = FindEntry(n1,n2);
      if (entry!=NULL) {
        map<cpl_info::code, ff_info * > * infos = entry->Get(prop);
        if (infos!=NULL) {
          if (infos->find(cpl)!=infos->end()) return (*infos)[cpl];
        }
      }
      return NULL;
    }
  public:
    Form_Factor_Parameter_Maps();
    ~Form_Factor_Parameter_Maps();

    void Initialize();
    void Output();

    Form_Factor_Base * GetFF(kf_code & n1,kf_code & n2,kf_code & prop,cpl_info::code & cpl);
    double GetModelParms(std::string constantName);
    bool ContainsFormFactorType(cpl_info::code & cpl);
  };

  extern Form_Factor_Parameter_Maps * ffs;

}
#endif
