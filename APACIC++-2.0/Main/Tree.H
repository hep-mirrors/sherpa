#ifndef Tree_H
#define Tree_H

#include "Run_Parameter.H"

#include "Knot.H"
#include "Poincare.H"
#include "Knot_List.H"
#include "Parton_List.H"

//using namespace APHYTOOLS;

namespace APACIC {
  class Tree {
    Knot      * root; 
    Knot_List * knots;

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
    Knot * NewKnot();
    /*! 
      Returns a pointer to a free knot with a halfway specified parton,
      i.e. momentum and flavour will be set.
    */
    Knot * NewKnot(APHYTOOLS::Flavour, AMATOOLS::Vec4D, double, double);
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
    //! deletes a Knot();
    void DeleteKnot(Knot *);
    //! Resets the complete tree, i.e. knots etc..
    void Reset(); 
    //! Resets all knots in a tree.
    void ResetKnots();
    //! Restores the tree down to partons with info = 'H', 'G', or 'M'
    void Restore();
    void Restore(Knot *);
    void ResetDaughters(Knot *);
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
    void ExtractPartons(APHYTOOLS::Parton_List *);
    //-----------------------------------------------------------------------
    //--------------------------- New frames for the tree -------------------
    //----------------------------------------------------------------------- 
    //! Boost the entire tree
    void BoRo(AMATOOLS::Poincare &);
    //! Boost the branch only
    static void BoRo(AMATOOLS::Poincare &, Knot * );
  private:
    static void BoRoDaughters(AMATOOLS::Poincare &, Knot *);
  };

  //-----------------------------------------------------------------------
  //--------------------------- Inline functions --------------------------
  //----------------------------------------------------------------------- 

  inline void Tree::Reset() {
    ResetKnots();    
    if (AORGTOOLS::rpa.gen.Debugging()) 
      std::cout<<"Passed resetting the tree."<<std::endl;
  }

  inline void Tree::ResetKnots() {
    for (Knot_Iterator kit=knots->begin(); kit!=knots->end(); ++kit) {
      delete (*kit); 
    }
    knots->erase(knots->begin(),knots->end());
    root = 0;
  }

}
#endif






