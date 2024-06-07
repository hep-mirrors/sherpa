#ifndef ATOOLS_Phys_NLO_Types_H
#define ATOOLS_Phys_NLO_Types_H

#include <string>
#include <iostream>

namespace ATOOLS {

  struct nlo_type {
    enum code {
      lo   = 0,
      born = 1,
      loop = 2,
      vsub = 4,
      real = 8,
      rsub = 16
    };
  };// end of struct nlo_type

  inline nlo_type::code operator~(const nlo_type::code c1)
  { return (nlo_type::code)(~(int)c1); }
  inline nlo_type::code operator|(const nlo_type::code c1,
                                  const nlo_type::code c2)
  { return (nlo_type::code)((int)c1|(int)c2); }
  inline const nlo_type::code &operator|=(nlo_type::code &c1,
                                          const nlo_type::code c2)
  { return c1=(nlo_type::code)((int)c1|(int)c2); }
  inline nlo_type::code operator&(const nlo_type::code c1,
                                  const nlo_type::code c2)
  { return (nlo_type::code)((int)c1&(int)c2); }
  inline const nlo_type::code &operator&=(nlo_type::code &c1,
                                          const nlo_type::code c2)
  { return c1=(nlo_type::code)((int)c1&(int)c2); }

  std::ostream &operator<<(std::ostream &str,const nlo_type::code &c);
  std::istream &operator>>(std::istream &str,nlo_type::code &c);

  struct asscontrib {
    enum type {
      none = 0,
      EW   = 1,
      LO1  = 2,
      LO2  = 4,
      LO3  = 8
    };
  };

  inline asscontrib::type operator|(const asscontrib::type at1,
                                 const asscontrib::type at2)
  { return (asscontrib::type)((int)at1|(int)at2); }
  inline const asscontrib::type &operator|=(asscontrib::type &at1,
                                         const asscontrib::type at2)
  { return at1=(asscontrib::type)((int)at1|(int)at2); }
  inline asscontrib::type operator&(const asscontrib::type at1,
                                 const asscontrib::type at2)
  { return (asscontrib::type)((int)at1&(int)at2); }
  inline const asscontrib::type &operator&=(asscontrib::type &at1,
                                         const asscontrib::type at2)
  { return at1=(asscontrib::type)((int)at1&(int)at2); }
  inline asscontrib::type operator^(const asscontrib::type at1,
                                 const asscontrib::type at2)
  { return (asscontrib::type)((int)at1^(int)at2); }
  inline const asscontrib::type &operator^=(asscontrib::type &at1,
                                         const asscontrib::type at2)
  { return at1=(asscontrib::type)((int)at1^(int)at2); }

  std::ostream & operator<<(std::ostream & s,const asscontrib::type & at);
  std::istream & operator>>(std::istream & s,asscontrib::type &at);

  struct nlo_mode {
    enum code {
      none       = 0,
      fixedorder = 1,
      powheg     = 2,
      mcatnlo    = 3,
      unknown    = 99
    };
  };// end of struct nlo_type

  std::ostream &operator<<(std::ostream &str,const nlo_mode::code &c);
  std::istream &operator>>(std::istream &str,nlo_mode::code &c);

  struct vtype {
    enum type {
      none          = 0,
      devidedByBorn = 1,
      includesI     = 2,
      needsCoupling = 4
    };
  };// end of struct nlo_type

  inline vtype::type operator|(const vtype::type vt1,
                               const vtype::type vt2)
  { return (vtype::type)((int)vt1|(int)vt2); }
  inline const vtype::type &operator|=(vtype::type &vt1,
                                       const vtype::type vt2)
  { return vt1=(vtype::type)((int)vt1|(int)vt2); }
  inline vtype::type operator&(const vtype::type vt1,
                               const vtype::type vt2)
  { return (vtype::type)((int)vt1&(int)vt2); }
  inline const vtype::type &operator&=(vtype::type &vt1,
                                       const vtype::type vt2)
  { return vt1=(vtype::type)((int)vt1&(int)vt2); }

  std::ostream &operator<<(std::ostream &str,const vtype::type &vt);
  std::istream &operator>>(std::istream &str,vtype::type &vt);

  struct cs_itype {
    enum type {
      none = 0,
      I    = 1,
      K    = 2,
      P    = 4
    };
  };// end of struct cs_itype

  inline cs_itype::type operator|(const cs_itype::type csit1,
                                  const cs_itype::type csit2)
  { return (cs_itype::type)((int)csit1|(int)csit2); }
  inline const cs_itype::type &operator|=(cs_itype::type &csit1,
                                          const cs_itype::type csit2)
  { return csit1=(cs_itype::type)((int)csit1|(int)csit2); }
  inline cs_itype::type operator&(const cs_itype::type csit1,
                                  const cs_itype::type csit2)
  { return (cs_itype::type)((int)csit1&(int)csit2); }
  inline const cs_itype::type &operator&=(cs_itype::type &csit1,
                                          const cs_itype::type csit2)
  { return csit1=(cs_itype::type)((int)csit1&(int)csit2); }

  std::ostream &operator<<(std::ostream &str,const cs_itype::type &csit);
  std::istream &operator>>(std::istream &str,cs_itype::type &csit);

  struct cs_kcontrib {
    enum type {
      none = 0,
      Kb   = 1,
      KFS  = 2,
      t    = 4,
      Kt   = 8
    };
  };// end of struct cs_kcontrib

  inline cs_kcontrib::type operator|(const cs_kcontrib::type csit1,
                                     const cs_kcontrib::type csit2)
  { return (cs_kcontrib::type)((int)csit1|(int)csit2); }
  inline const cs_kcontrib::type &operator|=(cs_kcontrib::type &csit1,
                                             const cs_kcontrib::type csit2)
  { return csit1=(cs_kcontrib::type)((int)csit1|(int)csit2); }
  inline cs_kcontrib::type operator&(const cs_kcontrib::type csit1,
                                     const cs_kcontrib::type csit2)
  { return (cs_kcontrib::type)((int)csit1&(int)csit2); }
  inline const cs_kcontrib::type &operator&=(cs_kcontrib::type &csit1,
                                             const cs_kcontrib::type csit2)
  { return csit1=(cs_kcontrib::type)((int)csit1&(int)csit2); }

  std::ostream &operator<<(std::ostream &str,const cs_kcontrib::type &csit);
  std::istream &operator>>(std::istream &str,cs_kcontrib::type &csit);

  struct sbt {
    enum subtype {
      none = 0,
      qcd  = 1,
      qed  = 2
    };
  };// end of struct sbt

  inline sbt::subtype operator|(const sbt::subtype s1,const sbt::subtype s2)
  { return (sbt::subtype)((int)s1|(int)s2); }
  inline const sbt::subtype &operator|=(sbt::subtype &s1,const sbt::subtype s2)
  { return s1=(sbt::subtype)((int)s1|(int)s2); }
  inline sbt::subtype operator&(const sbt::subtype s1,const sbt::subtype s2)
  { return (sbt::subtype)((int)s1&(int)s2); }
  inline const sbt::subtype &operator&=(sbt::subtype &s1,const sbt::subtype s2)
  { return s1=(sbt::subtype)((int)s1&(int)s2); }
  inline sbt::subtype operator^(const sbt::subtype s1,const sbt::subtype s2)
  { return (sbt::subtype)((int)s1^(int)s2); }
  inline const sbt::subtype &operator^=(sbt::subtype &s1,const sbt::subtype s2)
  { return s1=(sbt::subtype)((int)s1^(int)s2); }
  std::ostream &operator<<(std::ostream &ostr,const sbt::subtype &st);
  std::istream &operator>>(std::istream &str,sbt::subtype &st);

  struct subscheme {
    enum code {
      CS        = 0,
      Dire      = 1,
      CSS       = 2,
      Alaric    = 3
    };
  };// end of struct subtraction_variant

  std::ostream &operator<<(std::ostream &ostr,const subscheme::code &ss);
  std::istream &operator>>(std::istream &str,subscheme::code &ss);

  struct dpt {
    enum dipoletype {
      none = 0,
      f_f  = 1,
      f_i  = 2,
      i_f  = 3,
      i_i  = 4,
      f_fm = 5,
      f_im = 6,
      i_fm = 7
    };
  };// end of struct dpt

  std::ostream &operator<<(std::ostream &ostr,const dpt::dipoletype &dt);

  struct spt {
    enum splittingtype {
      // g - gluon, q - quark, s - squark, G - gluino, V - massive vector
      none  = 0,
      q2qg  = 1,
      q2gq  = 2,
      g2qq  = 3,
      g2gg  = 4,
      s2sg  = 5,
      s2gs  = 6,
      G2Gg  = 7,
      G2gG  = 8,
      V2Vg  = 9,
      V2gV  = 10
    };
  };// end of struct spt

  std::ostream &operator<<(std::ostream &ostr,const spt::splittingtype &st);

  struct ist {
    enum itype {
      // g - gluon, q - quark
      none  = 0,
      q     = 1,
      g     = 2,
      Q     = 3,
      V     = 4,
      sQ    = 5,
      sG    = 6
    };
  };// end of struct ist

  std::ostream &operator<<(std::ostream &ostr,const ist::itype &it);

}

#endif
