#ifndef CFColor_h
#define CFColor_h

#include "String_Tree.H"
#include "Single_Amplitude.H"

namespace AMEGIC {
  class CFColor {
    static std::string noname;

    String_Tree    st;
    Complex     ** CFC;

    int *     id;        // points to graph with same Color structure (mod. sign)
    int *     map;       // points to entry of color matrix for each graph
    int mcount, ncount;  // number of graphs, number of graphs with generic Color!

    void ReplaceT(sknot*);
    void ReplaceD(sknot*);
    void ReplaceG(sknot* m,sknot* m0=0);
    void ReplaceF(sknot*,char&);

    void Single_ReplaceF(sknot*,int&,char&);
    void Single_ReplaceFT(sknot*,int&,char&);

    int  Compare(Color_Function*,Color_Function*);
    int  CompareArg(int,int,int,Color_Function*,Color_Function*);
  public:
    /*!
      The constructor of the colour class that basically organizes the
      full evaluation of the colour matrix :
      - After initializing colour functions for every amplitude it is checked
      whether there are colour factors at all. If not, the colour matrix has
      only unity elements. 
      - Then the colour structures are analyzed in order to find identical ones
      to further alleviate the calculation.
      - Finally, pairs of surviving colour functions are evaluated.
      In doing so first structure constants f^abc are replaced by the usual
      relations :  
      Then, products of structure constants T^a_ij T^a_kl are replaced via :
      Furthermore, delta functions are eliminated.

    */
    CFColor(int,Single_Amplitude*,bool,std::string& pID=noname);
    ~CFColor();
    //! Outputs the colour matrix.
    void Output(std::string &);
    /*!
      Simple access to individual elements of the colour matrix.
    */
    inline const Complex& Mij(const int i,const int j) const {return CFC[i][j];}
    inline const int CFMap(int i) const {return map[i]; }
    inline const int CFSign(int i) const {return (id[i]>=0) ? 1 : -1; }
    inline const int MatrixSize() const {return ncount;}
  };
}
#endif
