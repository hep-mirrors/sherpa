#ifndef Tree_H
#define Tree_H

#include "Run_Parameter.H"

#include "Knot.H"
#include "Poincare.H"
#include "Particle_List.H"
#include <deque>

//using namespace ATOOLS;

namespace APACIC {

  typedef std::deque<Knot*> Knot_List;

  class Tree {
  private:
    Knot   * p_root, * p_save_root; 
    double   m_vetoscale;

    static Knot_List *s_knots;

    bool        CheckMomentumConservation(Knot *const knot) const;
    static void BoRoDaughters(ATOOLS::Poincare &, Knot *);

    Knot      * CopyKnot(Knot *,Knot *);
    void        CopyBackKnot(Knot *, Knot *);
    void        DeleteKnot(Knot *);

    friend std::ostream &operator<<(std::ostream &,const Tree &tree);

  public:

    // constructors
    Tree();
    Tree(Tree *tree);
    ~Tree();

    Knot * NewKnot(ATOOLS::Flavour, ATOOLS::Vec4D, double, double);
    Knot * NewKnot(ATOOLS::Particle * _inpart=NULL);
    Knot * NewKnot(Knot *);

    void   SetVetoScale(const double);
    double VetoScale() const;

    void Links(Knot *,Knot *);
    void Reset(); 
    void ResetKnots();
    void ClearStore();
    void Store();
    void Restore();

    Knot * GetRoot() const;
    Knot * GetInitiator() const;

    bool CheckStructure(bool fixit);
    bool SingleCheckStructure(Knot *mo, Knot*gr, bool fixit);
    void ExtractPartons(ATOOLS::Particle_List *);

    void        BoRo(ATOOLS::Poincare &);
    static void BoRo(ATOOLS::Poincare &, Knot * );
    static void UpdateDaughters(Knot *);
    bool        CheckMomentumConservation() const;
  private:
  };

  //-----------------------------------------------------------------------
  //--------------------------- Inline functions --------------------------
  //----------------------------------------------------------------------- 

  inline void   Tree::Reset()           { ResetKnots();  }
  inline Knot * Tree::GetRoot()   const { return p_root; }
  inline void   Tree::SetVetoScale(const double scale) { m_vetoscale = scale; }
  inline double Tree::VetoScale() const { return m_vetoscale; }

  /*!
    \file
    \brief contains the class APACIC::Tree
  */

  /*!
    \class Tree
    \brief  This class is used as a representation of a parton shower history
       in terms of binary splittings.

    This class provides all neccessary routines to handle a tree with
    binary splittings, esp. :
     - BoRo() ,
     - Store() ,
     - Restore() and
     - Reset()
  */


}
#endif






