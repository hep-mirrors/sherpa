#ifndef Cluster_Partons_CKKW_H
#define Cluster_Partons_CKKW_H

#include "Cluster_Partons_Base.H"

namespace SHERPA {

  class Cluster_Partons_CKKW: public Cluster_Partons_Base {
  private:

    std::vector<double> m_last_q;
    std::vector<int>    m_last_i;

    double m_as_jet[2];
    int    m_strong, m_nlegs;
    PDF::PDF_Base *p_pdf[2];
    
    void InitWeightCalculation();
    double CouplingWeight(const bool is,Leg &leg,const double &kt);
    double InternalWeight(const bool is,Leg &leg,const double upper,
			  const double actual,double min);
    double ExternalWeight(const bool is,Leg &leg,
			  const double actual,double min);
    void StoreOldValues(const int i,const int j,
			const int si,const double ptij);
    void WeightHardProcess();

    Combine_Table_Base *CreateTable(ATOOLS::Jet_Finder *jf,ATOOLS::Vec4D *amoms,Combine_Table_Base *ct,
				    const int isrmode,const int isrshoweron);

    void GenerateTables(const std::string &path);
    
  public :
    Cluster_Partons_CKKW(Matrix_Element_Handler *,ATOOLS::Jet_Finder * jf,
			 PDF::ISR_Handler * isr,const int,const int,const int);
    ~Cluster_Partons_CKKW();

    EXTRAXS::XS_Base * GetXS(EXTRAXS::XS_Group * group, ATOOLS::Flavour * fl);

    bool ClusterConfiguration(ATOOLS::Blob * blob,double x1=1.,double x2=1.);
    int  SetColours(EXTRAXS::XS_Base * xs, ATOOLS::Vec4D * p, ATOOLS::Flavour * fl);
    void CalculateWeight(const double &meweight);
  };
}

#endif
