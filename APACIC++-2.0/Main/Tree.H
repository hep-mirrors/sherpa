#ifndef Tree_H
#define Tree_H

#include "Run_Parameter.H"

#include "Knot.H"
#include "Poincare.H"
#include "Knot_List.H"
#include "Particle_List.H"

//using namespace ATOOLS;

namespace APACIC {
  class Tree {
    Knot      * p_root; 
    static Knot_List * s_knots;

    Knot      * p_save_root;
  public:
    friend std::ostream & operator<<(std::ostream &, Tree *);
    //-----------------------------------------------------------------------
    //--------------------------- Constructors ------------------------------
    //----------------------------------------------------------------------- 
    //! constructors
    Tree();
    Tree(Tree *);
    ~Tree();
    //-----------------------------------------------------------------------
    //--------------------------- Initialize New knots ----------------------
    //----------------------------------------------------------------------- 
    /*! 
      Returns a pointer to a free knot with a new parton.
      No links are set so far, in case no root of the knotlist was specified so far, 
      the new knot plays this role.
    */
    //    Knot * NewKnot();
    /*! 
      Returns a pointer to a free knot with a halfway specified parton,
      i.e. momentum and flavour will be set.
    */
    Knot * NewKnot(ATOOLS::Flavour, ATOOLS::Vec4D, double, double);
    Knot * NewKnot(ATOOLS::Particle * _inpart=NULL);
    /*! 
      Returns a pointer to a free knot in the list while copying another knot with
      zero links into it.
    */
    Knot * NewKnot(Knot *);
    /*! 
      Copies, if existing the links, of the incoming knot into the actual one.
      In doing that, new knots will be initialized.
    */
    void Links(Knot *,Knot *);
    //-----------------------------------------------------------------------
    //--------------------------- Deleting ----------------------------------
    //----------------------------------------------------------------------- 
    //! Resets the complete tree, i.e. knots etc..
    void Reset(); 
    //! Resets all knots in a tree.
    void ResetKnots();
    //! deletes a saved parton status (cf. Store())
    void ClearStore();
    //! saves a parton status (can be reactivated by Restore())
    void Store();
    //! Restores the tree down to partons with info = 'H', 'G', or 'M'
    void Restore();
    //    void ResetDaughters(Knot *);
    //-----------------------------------------------------------------------
    //--------------------------- Specific knots ----------------------------
    //----------------------------------------------------------------------- 
    //! Returns the root knot of the tree (perhaps iterator)
    Knot * GetRoot();
    //! Returns the initiator knot of the tree (i.e. the most previous one)
    Knot * GetInitiator();
    //-----------------------------------------------------------------------
    //--------------------------- Extracting partons ------------------------
    //----------------------------------------------------------------------- 
    bool CheckStructure(bool fixit);
    bool SingleCheckStructure(Knot *mo, Knot*gr, bool fixit);
    void ExtractPartons(ATOOLS::Particle_List *);
    //-----------------------------------------------------------------------
    //--------------------------- New frames for the tree -------------------
    //----------------------------------------------------------------------- 
    //! Boost the entire tree
    void BoRo(ATOOLS::Poincare &);
    //! Boost the branch only
    static void BoRo(ATOOLS::Poincare &, Knot * );
    static void UpdateDaughters(Knot *);
  private:
    static void BoRoDaughters(ATOOLS::Poincare &, Knot *);

    // some recursive helpers
    Knot * CopyKnot(Knot *,Knot *);
    void CopyBackKnot(Knot *, Knot *);
    //! deletes a Knot();
    void DeleteKnot(Knot *);
  };

  //-----------------------------------------------------------------------
  //--------------------------- Inline functions --------------------------
  //----------------------------------------------------------------------- 

  inline void Tree::ResetKnots() {
    if (!s_knots) return;
    for (Knot_Iterator kit=s_knots->begin(); kit!=s_knots->end(); ++kit) {
      delete (*kit); 
    }
    s_knots->erase(s_knots->begin(),s_knots->end());
    p_root = 0;
  }

  inline void Tree::Reset() {
    ResetKnots();    
  }

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






