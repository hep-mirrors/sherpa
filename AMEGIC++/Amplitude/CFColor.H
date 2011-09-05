#ifndef AMEGIC_Amplitude_CFColor_H
#define AMEGIC_Amplitude_CFColor_H

#include "ATOOLS/Math/Matrix.H"
#include "AMEGIC++/String/String_Tree.H"
#include "AMEGIC++/Amplitude/Single_Amplitude.H"
#include <map>

namespace AMEGIC {

  //!the container for the T and F products
  typedef std::map<std::string, Complex>  TFprod_Table;
  typedef TFprod_Table::iterator          TF_Iterator;
  //!the map of indices
  typedef std::map<char, char>            Char_Map;
  typedef std::map<char, int>             CharNum_Map;
  
  class CFColor {
    static std::string noname;
    
    String_Tree    st;
    ATOOLS::CMatrix * CFC;

    int *     id;        // points to graph with same Color structure (mod. sign)
    int *     map;       // points to entry of color matrix for each graph
    int mcount, ncount;  // number of graphs, number of graphs with generic Color!

    void ReplaceT(sknot*);
    void ReplaceD(sknot*,sknot*);
    void ReplaceG(sknot* m,sknot* m0=0);
    void ReplaceF(sknot*,CharNum_Map&,char&);


    int  SingleReplaceD(sknot *,sknot *,char, char);
    void SingleReplaceF(sknot*,int&,CharNum_Map&,char&);
    void SingleReplaceFT(sknot*,int&,CharNum_Map&,char&);

    int  Compare(MODEL::Color_Function*,MODEL::Color_Function*);
    int  CompareArg(int,int,int,MODEL::Color_Function*,MODEL::Color_Function*);
    
    void ExtractIndices(sknot*,CharNum_Map&);
    char DeliverIndex(CharNum_Map&,char&);
    std::string BuildTChain(std::vector<std::string>);
    std::string MapFChain(std::vector<std::string>);
    TFprod_Table t_table;
    TFprod_Table f_table;
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
    CFColor(int,Single_Amplitude*,ATOOLS::Flavour *,char,char,std::string pID=noname,bool=false);
    ~CFColor();
    //! Outputs the colour matrix.
    void Output(std::string &);
    /*!
      Simple access to individual elements of the colour matrix.
    */
    inline const Complex& Mij(const int i,const int j) const {return (*CFC)[i][j];}
    inline const ATOOLS::CMatrix* GetCMatrix() const {return CFC;}
    inline int CFMap(int i) const {return map[i]; }
    inline int CFSign(int i) const {return (id[i]>=0) ? 1 : -1; }
    inline int MatrixSize() const {return ncount;}
  };
  
}
#endif
