#ifndef Combine_Table_H
#define Combine_Table_H

#include <vector>
#include <map>
#include "Point.H" 
#include "Vector.H"

namespace ATOOLS {
  class Jet_Finder;
}


namespace SHERPA {
  // Look for more comments at the end of this file!

  /*---------------------------------------------------------------------------

                             Leg                       

    ---------------------------------------------------------------------------*/


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

  std::ostream & operator<<(std::ostream &,const Combine_Data &);

  /*---------------------------------------------------------------------------

                             Combine_Data                       

    ---------------------------------------------------------------------------*/
  //! Stores the data of one row of the Combine_Table.
  class Combine_Key {
    int      i,j;             // combine leg "i" and "j"
    ATOOLS::Flavour  flav;    // results in propagator flavor, only for ew merging
  public:
    Combine_Key();
    Combine_Key(int _i, int _j, const ATOOLS::Flavour & _flav=ATOOLS::Flavour(ATOOLS::kf::none));
    friend class Combine_Table;
    friend std::ostream& operator<< (std::ostream&,const Combine_Key &);
    friend std::ostream& operator<< (std::ostream&,const Combine_Table &);
    friend bool operator<(const Combine_Key & a, const Combine_Key & b);
  };

  class Combine_Data {
    double  pt2ij;
    double  sij;
    double  prop;
    double  coupling;
    double  weight;
    int     strong;    // there is a strong graph for this

    std::vector<int>   graphs;    // all graphs where this combination is possible
    Combine_Table    * down;      // new Combine table with n-1 legs (i and j combined)
  public:
    Combine_Data();
    Combine_Data(const Combine_Data &);
    Combine_Data(double, int);
    ~Combine_Data();
    friend class Combine_Table;
    friend std::ostream& operator<< (std::ostream&,const Combine_Data &);
    friend std::ostream& operator<< (std::ostream&,const Combine_Table &);
  };




  //! standard list and iterator for Combine_Data
  typedef  std::map<Combine_Key,Combine_Data> CD_List;
  typedef  CD_List::iterator                  CD_Iterator;
  typedef  CD_List::const_iterator            CD_Const_Iterator;

  struct Vertex_Info {
    ATOOLS::Flavour fl;
    Complex * cpl;
    AMEGIC::Color_Function * color;
    unsigned int mode;
    inline Vertex_Info();
  };

  std::ostream & operator<<(std::ostream &,const Combine_Table &);

  class Combine_Table {
    int        m_mode;
    static int        s_all;
    int               m_no;

    Combine_Table  *  p_up;
    Leg            ** p_legs;
    int               m_nlegs, m_nampl;
    CD_List           m_combinations;

    CD_Iterator       m_cwin;
    int               m_gwin;

    int               m_isr1on, m_isr2on, m_isrshoweron;

    ATOOLS::Jet_Finder  * p_jf;
    ATOOLS::Vec4D       * p_moms;

    double m_x1, m_x2;

    double m_kt2min;

    bool  Combinable(const Leg &, const Leg &, 
		     Vertex_Info &, Vertex_Info &) const ;
    static Leg * CombineLegs(Leg *, int, int, int );
    static Leg   CombinedLeg(Leg *, int, int);

    void CombineMoms(ATOOLS::Vec4D *, int, int, int);
    static void CombineMoms(ATOOLS::Vec4D *, int, int, int, ATOOLS::Vec4D * &);

    void AddPossibility(int i, int j, int ngraph, 
			const Vertex_Info & vinfo1, const Vertex_Info & vinfo2);
    double ColorFactor(int i, int j, AMEGIC::Color_Function * const color,
		       unsigned int mode);
    double Sprime() const;
  public:
    Combine_Table(ATOOLS::Jet_Finder * jf, ATOOLS::Vec4D * moms, 
		  Combine_Table * up, int isrmode, int isrshoweron, int mode);
    ~Combine_Table();


    void            FillTable(Leg **,int, int); 
    Combine_Table * CalcJet(int, double, double, ATOOLS::Vec4D * moms=0);
    CD_Iterator     CalcPropagator(CD_Iterator & cit);

    //! general access routines (return value is determined by winner combination 
    //! and winner graph number)
    //@{
    const int   NLegs()                   { return m_nlegs; }
    const Leg & GetLeg(int ngraph, int l) { return p_legs[ngraph][l]; }
    const Leg & GetLeg(int l)             { return p_legs[m_gwin][l]; }
    double GetWinner(int & i, int & j);
    int AddCouplings(int & nqed, int & nqcd);
    double MinKt();
    void GetX1X2(double & x1, double & x2) {
      x1=m_x1;
      x2=m_x2;
    }

    ATOOLS::Flavour  Flav(int i)      { return p_legs[m_gwin][i].Flav(); }
    ATOOLS::Vec4D  & Momentum(int i)  { return p_moms[i]; }
    Combine_Table  * Up()             { return p_up; }
    //@}

    // only friends can access internal data
    friend std::ostream & operator<< (std::ostream&,const Combine_Table &);
    friend std::ostream & operator<< (std::ostream&,const Combine_Data &);
  };

  inline Vertex_Info::Vertex_Info() :
    cpl(0), color(0), mode(0) {}

  /*!
    \file Combine_Table.H
    \brief Defines the basic structures for the combination of matrix elements and 
           parton shower. Therefore a class SHERPA::Combine_Table, SHERPA::Combine_Data, 
	   and SHERPA::Leg is defined, in order to determine and store all possiblities 
 	   to combine two legs of a given matrix element. The information of all graphs 
	   (used in AMEGIC) is employed to seperate allowed combination. The aim is to 
	   have all neccessary information to establish a parton histories for a n-jet 
	   matrix element configuration easily.
  */

  /*! 
    \class Leg
    \brief Stores all information of a point and an extra "anti"-flag
    
    A Leg represents a particle dangling from a Feynman diagram, either before or
    while the clustering takes place. It stores all information of an AMEGIC::Point and an 
    extra "anti"-flag.

    The idea is to access, all information of the Point with same syntax as an pointer 
    to an AMEGIC:: Point. For this reason the operator-> is overloaded. In addition a anti-flag 
    is stored. This is due to keep track of the information in case a propagator has 
    the "wrong".
    
    Note: The overloaded operator-> can sometimes lead to confusion, especially the 
          anti flag can not be accessed via this operator in case a pointer to a leg is 
          used. In this case the operator* together with the dot has to be used.
          So always think "Leg" as being synonymous to "Point*".
  */
  /*! 
    \class Combine_Key
    \brief is one of the basic elements for the creation of a Combine_Table.
  
    It includes the numbers of combinable legs (i and j) and the flavour of the resulting propagator. 
    It is used as a key in a fast access map in a Combine_-Table in order to access the
    information placed in a Combine_Data object.
  */
  /*!
    \struct Vertex_Info 
    \brief A wrapper containing relevant information of AMEGIC::Vertex
  */
  /*! 
    \class Combine_Data
    \brief is the basic element for the determination of clustering when using a Combine_Table.
 
    This includes the number of combined leg (i and j), the jet measure of these legs 
    (yij), a list of (numbered) graphs where this combination is possible (graphs) and a 
    link to the new table where those legs have been combined (down).

    In case the extended merging algorithm is active, additional information is
    included, namely: the virtuality of the resulting propagator sij, the estimate of the propagator prop, 
    and the coupling of the corresponding vertices.
  */
  /*!
    \fn Combine_Data::Combine_Data(double _pt2ij, int _ngraph)
    \brief Constructor with a value for the kperp and one graph (the first one), where the corresponding
    combination of particles appears.
  */


  /*! 
    \class Combine_Table
    \brief This is the main structure for storing a combination scheme.
    
    This class provides the structure for storing histories of successive
    clusterings. The structure fills itself recursively, with the only
    input being the Feynman diagrams of the process under consideration
    and the four-momenta of the current event. A Combine_Table consists of 
    a list of possible clusterings (particles i, j and the flavour of the resulting 
    intermediate particle) and the kperp values and Feynman diagrams 
    associated with them. This information is realised through the classes Combine_Key 
    and Combine_Data. In addition, a number of methods allows a Combine_Table to create these 
    data and to construct the sequence of Combine_Tables representing the clustering
    history.

    A Combine_Table stores: 
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
  /*!
    \var Combine_Table::m_mode
    
    m_mode =
    -# 0 : Standard CKKW merging
    -# 1 : electroweak merging a la hep-ph/0503281
  */
  /*! 
    \fn Combine_Table::FillTable()  
    has two tasks to fulfil. First of all, a set of given Legs,
    i.e. particles, are filled into the table. Then, all pairs of them are checked whether 
    they can be clustered or not. A clustering is possible only, if it occurs in 
    a corresponding Feynman diagram, which disables unphysical parton histories. This 
    check is performed through the method Combinable().
  */
  /*! 
    \fn Combine_Table::CalcJet()
    evaluates the kperp distance of all allowed parton pairs (i,j) 
    created by FillTable() with the Jet_Finder. After that,
    a pair to be clustered is selected according to the merging
    prescription, and a new Combine_Table is constructed, where
    the number of Legs is reduced by one. Consequently, after each
    clustering step, the set of Feynman diagrams is pruned, to include only
    those where the selected combination is possible. The four-momentum of
    the new (joined) Leg is given by the corresponding combination
    of the two individual four momenta. The algorithm continues
    recursively with corresponding calls to FillTable() and 
    CalcJet() until only a 2\to 2 process remains.
  */
  /*! 
    \fn Combine_Table::CalcPropagator()
    performs all basic calculations for the determination of cluster
    probability for a given pair (i, j). This usually includes the
    evaluation of the kperp measure, and the invariant mass
    sij. In case of the extended clustering algorithm, an estimate
    for the branching probability is also computed, which includes the
    couplings of that branching process, as well as
    the corresponding propagator. The couplings are available in the
    Feynman diagrams provided by AMEGIC.
  */
  /*! \fn Combine_Table::Combinable()
    determines whether two particles, i.e. Legs can be
    clustered. To this end, the two Points related to the 
    Legs are checked whether they have a common third Point,
    i.e. vertex, linked to them.    
  */
  /*! \fn Combine_Table::AddPossibility(int i, int j, int ngraph,
                                        const Vertex_Info & vinfo1,const Vertex_Info & vinfo2)
      adds possible clusterings ij with the corresponding graph number ngraph into the CD_List 
      (map of Combine_Data). If m_mode=0, it is always assumed that the resulting intermediate particle is 
      unambigous - this is the original CKKW approach which assumes a fixed minimal number of electroweak 
      interactions. If, in contrast m_mode=1, i.e. the electroweak merging is the mehtod of choice, then 
      couplings need to be added to the Combine_Data, and the Combine_Key also carries information about the 
      resulting propagator flavour, accessible through the Vertex_Info. Then, for each different propagator 
      flavour for the same clustering ij a new Combine_Key entry is created.
   */
}


#endif



