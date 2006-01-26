#ifndef Cluster_Partons_CKKW_H
#define Cluster_Partons_CKKW_H

#include "Cluster_Partons_Base.H"

namespace SHERPA {

  class Cluster_Partons_CKKW: public Cluster_Partons_Base {
  private:
    double m_AcceptMisClusters;
    bool   m_LowestFromME;

    std::vector<double> m_last_q;
    std::vector<int>    m_last_i;

    double m_as_amegic, m_as_jet, m_as_hard, m_as_qcd;
    int    m_strong, m_nlegs, m_count_startscale;
    
    void InitWeightCalculation();
    bool ApplyCombinedInternalWeight(const bool is,const ATOOLS::Flavour & fl,
				     const double upper,const double actual,
				     const double asref,const int order);
    bool ApplyExternalWeight(const bool is,const ATOOLS::Flavour & fl,const double actual);
    void StoreOldValues(const int i,const int j,
			const int si,const double ptij);
    void OverwriteScales(const int);

    Combine_Table_Base *CreateTable(ATOOLS::Jet_Finder *jf,ATOOLS::Vec4D *amoms,Combine_Table_Base *ct,
				    const int isrmode,const int isrshoweron);

    void GenerateTables(const std::string &path);
    
  public :
    Cluster_Partons_CKKW(Matrix_Element_Handler *,ATOOLS::Jet_Finder * jf,const int,const int,const int,const int);
    ~Cluster_Partons_CKKW();

    EXTRAXS::XS_Base * GetXS(EXTRAXS::XS_Group * group, ATOOLS::Flavour * fl);

    bool ClusterConfiguration(ATOOLS::Blob * blob,double x1=1.,double x2=1.);
    int  SetColours(EXTRAXS::XS_Base * xs, ATOOLS::Vec4D * p, ATOOLS::Flavour * fl);
    void CalculateWeight();
  };
}

#endif
