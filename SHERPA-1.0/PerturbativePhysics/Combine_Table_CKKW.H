#ifndef Combine_Table_CKKW_H
#define Combine_Table_CKKW_H

#include "Combine_Table_Base.H"

namespace SHERPA {

  class Combine_Table_CKKW: public Combine_Table_Base {
  private:

    ATOOLS::Vec4D *p_save_moms;
    size_t         m_nl;

    bool Combinable(const Leg &a,const Leg &b) const ;
    void AddPossibility(const int i,const int j,const int ngraph);
    bool InitStep(ATOOLS::Vec4D *moms,const int nl);
    bool SelectWinner(const bool did_boost);
    bool TestMomenta(const int i,const int j);

    Combine_Table_CKKW *CreateNext(const bool did_boost);
    Combine_Table_CKKW *NextTable(Combine_Table_CKKW *tab,
				  const double x1,const double x2);

  public:

    //constructor
    Combine_Table_CKKW(ATOOLS::Jet_Finder * jf, ATOOLS::Vec4D * moms, 
		       Combine_Table_CKKW * up, int isrmode, int isrshoweron);
    // destructor
    ~Combine_Table_CKKW();

    // member functions
    CD_List::iterator   CalcPropagator(CD_List::iterator &cit);
    Combine_Table_Base *CalcJet(int nl,const double x1,const double x2,
				ATOOLS::Vec4D *moms=NULL);

    void   FillTable(Leg **legs,const int nlegs,const int nampl); 
    double GetWinner(int &i, int &j);

  };// end of class Combine_Table_CKKW

}// end of namespace SHERPA

#endif



