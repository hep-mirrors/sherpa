#ifndef Combine_Table_H
#define Combine_Table_H

#include <vector>
#include <list>
#include "Jet_Finder.H"  
#include "Point.H" 
#include "Particle.H" 



namespace SHERPA {
  // Look for more comments at the end of this file!

  /*---------------------------------------------------------------------------

                             Leg                       

    ---------------------------------------------------------------------------*/
  //! Stores all information of a point and an extra "anti"-flag


  class Leg {  // can be used as synonym for "Point *"
    AMEGIC::Point     * point;
    int                 anti;    // extra anti for propagator
  public:
    // constructor
    Leg(AMEGIC::Point * _point = 0, int _anti = 1) :  
      point(_point), anti(_anti) { }
    Leg(const Leg & l) : point(l.point), anti(l.anti) { }

    // access function for point data
    const AMEGIC::Point * operator&()  const { return point; }
    const AMEGIC::Point * operator->() const { return point; }

    ATOOLS::Flavour Flav()              
    { 
      if (anti==-1) return point->fl.Bar(); 
      return point->fl;
    }
    // access extra anti
    int   ExtraAnti()                const { return anti; }
    void  SetAnti(int a)                   { anti   = a; }
    void  SetPoint(AMEGIC::Point * p)      { point  = p; }
  };


  class Combine_Data;
  class Combine_Table;

  std::ostream & operator<<(std::ostream &,Combine_Data &);

  /*---------------------------------------------------------------------------

                             Combine_Data                       

    ---------------------------------------------------------------------------*/
  //! Stores the data of one row of the Combine_Table.
  class Combine_Data {
    int                i,j;       // combine leg "i" and "j"
    double             pt2ij;
    int                strong;    // there is a strong graph for this
    std::vector<int>   graphs;    // all graphs where this combination is possible
    Combine_Table    * down;      // new Combine table with n-1 legs (i and j combined)
  public:
    Combine_Data();
    Combine_Data(const Combine_Data &);
    Combine_Data(int, int, double, int);
    ~Combine_Data();

    // only friends can access the internal data
    friend class Combine_Table;
    friend std::ostream& operator<< (std::ostream&,Combine_Data &);
    friend std::ostream& operator<< (std::ostream&,Combine_Table *);
  };

  //! standard list and iterator for Combine_Data
  typedef  std::list<Combine_Data>  CD_List;
  typedef  CD_List::iterator        CD_Iterator;




  std::ostream & operator<<(std::ostream &,Combine_Table *);

  //! This is the main structure for the combination scheme.
  class Combine_Table {
    static int        all;
    int               no;

    Combine_Table  *  up;
    Leg            ** legs;
    int               nlegs, nampl;
    CD_List           combinations;

    CD_Iterator       cwin;
    int               gwin;

    int               m_isr1on, m_isr2on, m_isrshoweron;

    ATOOLS::Jet_Finder * jf;
    ATOOLS::Vec4D       * moms;

    double x1,x2;

    bool  Combinable(const Leg &, const Leg &, int & strong) const ;
    Leg * CombineLegs(Leg *, int, int, int );

    void CombineMoms(ATOOLS::Vec4D *, int, int, int);
    void CombineMoms(ATOOLS::Vec4D *, int, int, int, ATOOLS::Vec4D * &);

    void AddPossibility(int, int, int, int);
    Leg  CombinedLeg(Leg *, int, int);
    double Sprime();
  public:
    Combine_Table(ATOOLS::Jet_Finder *, ATOOLS::Vec4D *, Combine_Table *, int, int);
    ~Combine_Table();


    void            FillTable(Leg **,int, int); 
    Combine_Table * CalcJet(int, double, double, ATOOLS::Vec4D * _moms=0);

    //! general access routines (return value is determined by winner combination 
    //! and winner graph number)
    //@{
    const int   NLegs()                   { return nlegs; }
    const Leg & GetLeg(int ngraph, int l) { return legs[ngraph][l]; }
    const Leg & GetLeg(int l)             { return legs[gwin][l]; }
    double GetWinner(int & i, int & j)    { 
      i = cwin->i;j = cwin->j;
      return sqrt(cwin->pt2ij);
    }
    void GetX1X2(double & _x1, double & _x2) {
      _x1=x1;
      _x2=x2;
    }

    ATOOLS::Flavour   Flav(int i)      { return legs[gwin][i].Flav(); }
    ATOOLS::Vec4D    & Momentum(int i)  { return moms[i]; }
    Combine_Table      * Up()             { return up; }
    //@}

    // only friends can access internal data
    friend std::ostream & operator<< (std::ostream&,Combine_Table *);
    friend std::ostream & operator<< (std::ostream& s ,Combine_Data& cd);
  };
}





/*!
  \file Combine_Table.H
  \brief Defines the basic structures for the combination of matrix elements and 
         parton shower. Therefore a class APACIC::Combine_Table, APACIC::Combine_Data, 
	 and APACIC::Leg is defined, in order to determine and store all possiblities 
	 to combine two legs of a given matrix element. The information of all graphs 
	 (used in AMEGIC) is employed to seperate allowed combination. The aim is to 
	 have all neccessary information to establish a parton histories for a n-jet 
	 matrix element configuration easily.
*/

/*! 
  \class APACIC::Leg

  The idea is to access, all information of the Point with same syntax as an pointer 
  to a Point. For this reason the operator-> is overloaded. In addition a anti-flag 
  is stored. This is due to keep track of the information in case a propagator has 
  the "wrong".

  Note: The overloaded operator-> can sometimes lead to confusion, especially the 
        anti flag can not be accessed via this operator in case a pointer to a leg is 
	used. In this case the operator* together with the dot has to be used.
        So always think "Leg" as a synonym for "Point*".
*/

/*! 
  \class APACIC::Combine_Data

  This includes the number of combined leg (i and j), the jet measure of these legs 
  (yij), a list of numbers graphs where this combination is possible (graphs) and a 
  link to the new table where those legs have been combined (down).
*/

/*! 
  \class APACIC::Combine_Table

  it stores: 
    - a table with (e.g. 6) combinations: (i and j are leg numbers)
      \verbatim
       i&j   y_ij         graphs         pointer to next table
       -------------------------------------------------------
       0&2   0.0810366    2,3,4,5,8      #
       0&4   0.0691623    6,7            8
       1&3   0.0844243    0,1,6,7,8      #
       1&4   0.293399     4,5            #
       2&4   0.111385     0,1            #
       3&4   0.215127     2,3            #
      \endverbatim
    - a list of corresponding momenta
      \verbatim
       0    (45.5942,0,0,45.5942)
       1    (27.4694,-8.00725,4.40136,-29.941)
       2    (21.7338,-13.9122,15.4912,6.23155)
       3    (21.9827,18.2065,7.01464,-10.127)
       4    (29.347,-12.3016,-18.1044,19.5487)
      \endverbatim
    - a matrix of all dangling legs (one row for each graph)
    - a link to a jet_measure calculator (Jet_Finder)
    - a reference to the winner combination and graph number
    - a link "up" to the table with one combination less performed

*/

#endif



