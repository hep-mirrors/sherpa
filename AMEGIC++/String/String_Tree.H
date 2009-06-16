#ifndef AMEGIC_String_String_Tree_H
#define AMEGIC_String_String_Tree_H

#include "AMEGIC++/String/MyString.H"
#include "ATOOLS/Math/Kabbala.H"
#include "ATOOLS/Math/MyComplex.H"
#include <list>
#include <vector>

namespace AMEGIC {
  class sknot {
    std::string*   strp;
    static std::string emptystring;
  public:
    sknot*   left;
    sknot*   right;
    ATOOLS::Kabbala* value;
    char     op;

    sknot() {strp = 0; value =0;}
    sknot(const sknot& cp) {
      left  = cp.left;
      right = cp.right;
      value = cp.value;
      op    = cp.op;
      if (cp.strp!=0) strp = new std::string(cp.Str());
                 else strp = 0;
    }
    ~sknot() {if (strp!=0) delete strp;}
    void KillString() {if (strp!=0) delete strp; strp=0;}
    void KillValue()  {if (value!=0) delete value; value=0;}
    const std::string& Str() const { if (strp!=0) return *strp;return emptystring;}
    void SetString(const std::string& str) {
      if (strp!=0) delete strp;
      strp = new std::string(str);
    }
    sknot& operator=(const sknot& cp) {
      if (this!=&cp) {
      left  = cp.left;
      right = cp.right;
      value = cp.value;
      op    = cp.op;
      if (strp!=0) delete strp;
      if (cp.strp!=0) strp = new std::string(cp.Str());
                 else strp = 0;
      }
      return *this;
    }
  };

  class String_Tree {
    static sknot zero;
    static const int block_size;
    int    skpos;
    std::vector<sknot*> sblocks;
    std::list<sknot*> leaflist;
    sknot* Leaf(std::string&,sknot*,int);
   
    void Single_Expand(sknot*,int&);
    void Single_Delete(sknot*&,sknot*,const std::string&); 

    void SingleDeleteMinus(sknot* &,int&);
    int  CountMinus(sknot* &); 


    void    LinearPM(sknot*);
    void    OrderPM(sknot*,sknot*);
    void    DetermineLeafAndSign(sknot*,std::vector<sknot*>&,std::vector<int>&,int&);
    void    SetLeafAndSign(sknot*,std::vector<sknot*>&,std::vector<int>&,int&);
    void    DeleteEquals(std::vector<sknot*>&,std::vector<int>&);

    int     CountFactorNumber(sknot*,std::vector<sknot*>*&,sknot*,std::vector<sknot*>*&,int);
    void    CollectLeafs(sknot*,std::vector<sknot*>&,int);
    
  public:
    String_Tree();
    ~String_Tree();
    void Reset();
    sknot*  newsk(); 
    void    popsk(); 
    sknot*  String2Tree(std::string,int fixed = 0); 
    Complex eval(sknot*);
    Complex evalcolor(sknot*);
    Complex Evaluate(sknot*);
    std::string  Tree2String(sknot*,sknot*);
    std::string  Tree2Tex(sknot*,sknot*);
    void    Expand(sknot*);
    void    DetermineDepth(sknot*,char,int&);
    void    ExpandToDepth(sknot*,int,std::list<sknot*>&);
    void    Linear(sknot*);
    void    LinearOrderPM(sknot*);
    sknot*  Copy(sknot*,int fixed = 0);
    void    Sort(sknot*);
    void    Delete(sknot*&,const std::string&); 
    void    GetEnd(sknot*,std::list<sknot*> &);
    void    Find(sknot*,const std::string&,int&);
    void    Simplify(sknot*&);
    void    DeleteMinus(sknot* &);
    int     SknotListSize() {return skpos+1;}
    void    Cluster(sknot*,sknot*,int full = 0);
    void    Addends(sknot*,std::list<sknot*> &);
    void    Factors(sknot*,std::list<sknot*> &);
    void    CleanValues();
  };
}
#endif
