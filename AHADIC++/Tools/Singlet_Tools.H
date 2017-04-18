#ifndef AHADIC_Tools_Singlet_Tools_H
#define AHADIC_Tools_Singlet_Tools_H

#include "AHADIC++/Tools/Proto_Particle.H"
#include "AHADIC++/Tools/Constituents.H"
#include <list>

namespace AHADIC {
  class Singlet : public std::list<Proto_Particle *> {
    bool ReorderCriterion(Proto_Particle * first=NULL);
  private:
  public:
    Singlet();
    ~Singlet();
    void Reorder(Proto_Particle * first=NULL);
    void Erase(Proto_Particle * ref);
    bool Combine(Proto_Particle * part1,Proto_Particle * part2);
    double        Mass2() const;
    ATOOLS::Vec4D Momentum() const;
  };

  std::ostream & operator<<(std::ostream &,const Singlet &);

  class Singlet_Tools {
  protected:
    Constituents * p_constituents;
    Singlet      * p_singlet;

    double m_minQmass, m_mass;
  public:
    Singlet_Tools();
    ~Singlet_Tools() {}
    virtual void Init();
    
    virtual bool CheckMass(Proto_Particle * part1,Proto_Particle * part2);
 };
}
    
#endif
