#ifndef AMEGIC_Cluster_Cluster_Algorithm_H
#define AMEGIC_Cluster_Cluster_Algorithm_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "AMEGIC++/Cluster/Combine_Table.H"
#include "EXTRA_XS/Main/Single_Process.H"

namespace ATOOLS { class Cluster_Definitions_Base; }

namespace AMEGIC {

  class Process_Base;

  class Cluster_Algorithm {
  public:

    typedef std::map<ATOOLS::Flavour_Vector,EXTRAXS::ME2_Base*> Flav_ME_Map;

  private:

    ATOOLS::Mass_Selector *p_ms;

    Process_Base *p_proc;

    ATOOLS::Cluster_Amplitude        *p_ampl;
    ATOOLS::Cluster_Definitions_Base *p_clus;

    Combine_Table *p_ct, *p_combi;

    std::vector<double> m_qmin;

    Flav_ME_Map m_xsmap;

    int m_colors[4][2];

    EXTRAXS::ME2_Base *p_xs;

    bool SetColors();
    void Convert();

    bool FillLegs(Leg * alegs, Point * root, int & l, int maxl);
    Leg **CreateLegs(int &nampl,const int nlegs);
    void CreateTables(Leg ** legs,const int nampl,
		      const double x1,const double x2,const size_t mode) ;

    void SetNMax(ATOOLS::Cluster_Amplitude *const ampl,
		 const size_t &id,const size_t &nmax) const;

    void SetColours(ATOOLS::Cluster_Leg *const lij,
		    ATOOLS::Cluster_Leg *const li,
		    ATOOLS::Cluster_Leg *const lj);

    int SetDecayColours(const ATOOLS::Vec4D_Vector& p, ATOOLS::Flavour * fl,
			int col1,int col2);
    int SetColours(const ATOOLS::Vec4D_Vector& p,ATOOLS::Flavour * fl);
    int Set4Colours(const int nquark,const int ngluon,
		    const ATOOLS::Vec4D_Vector& p,ATOOLS::Flavour * fl);
    int Set3Colours(const int nquark,const int ngluon,
		    const ATOOLS::Vec4D_Vector& p,ATOOLS::Flavour * fl);
    int Set2Colours(const int nquark,const int ngluon,
		    const ATOOLS::Vec4D_Vector& p,ATOOLS::Flavour * fl);

    int SetColours(EXTRAXS::ME2_Base * xs, 
		   const ATOOLS::Vec4D_Vector& p, ATOOLS::Flavour * fl);
    EXTRAXS::ME2_Base *GetXS(const ATOOLS::Flavour_Vector &fl);

  public:

    Cluster_Algorithm(ATOOLS::Mass_Selector *const ms);

    ~Cluster_Algorithm();

    bool Cluster(Process_Base *const xs,const size_t mode);

    // inline functions
    inline ATOOLS::Cluster_Amplitude *Amplitude() { return p_ampl; }

    inline void SetClusterDefinitions
    (ATOOLS::Cluster_Definitions_Base *const cb) { p_clus=cb; }

  };// end of class Amegic_Cluster_Algorithm

}// end of namespace AMEGIC

#endif
