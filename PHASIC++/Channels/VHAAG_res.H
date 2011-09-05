#ifndef PHASIC_Channels_VHAAG_res_h
#define PHASIC_Channels_VHAAG_res_h

#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Vegas.H"
#include <map>

namespace PHASIC {
  class VHAAG_res : public Single_Channel {
    int      n_p1,m_type,n_d1,n_d2,n_b,m_rkf,n_ap;
    int      *p_perm;
    double   m_s0;
    ATOOLS::Vec4D* m_q;
    double*  p_s;
    Vegas* p_vegas;
    bool m_ownvegas;
    
    std::map<int,Vegas*>* p_sharedvegaslist;

    double PiFunc(double a1,double a2,
		  double s1b,double s2b,double c);
    void Split(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
	       ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int,int,double *ran);
    void Split0(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int,int,double *ran);
    void SingleSplit(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		     ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double,double *ran);
    void SingleSplitF(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		      ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double,double *ran);
    void SingleSplitF0(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		       ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double,double *ran);
    double SplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
		       ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int,int,double *ran);
    double Split0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int,int,double *ran);
    double SingleSplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D& Q,
			     ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double,double *ran);
    double SingleSplitFWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D& Q,
			      ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double,double *ran);
    double SingleSplitF0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			      ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double,double *ran);
    void GenerateBosonMass(ATOOLS::Vec4D *p,double *ran);
    double BosonWeight(ATOOLS::Vec4D *p,double *ran);
    void GenerateBranch(ATOOLS::Vec4D q1,ATOOLS::Vec4D Q,
			ATOOLS::Vec4D* q,double*,int n,double *ran);
    double BranchWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D &Q,
			ATOOLS::Vec4D* q,double*,int n,double *ran);
    void ConstructMomenta(double a1,double a2,
			  double s1,double s2,double s,
			  ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			  ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2);
    void ConstructMomenta(double a1,double phi,
			  double s1,double s2,double s,
			  ATOOLS::Vec4D q1,ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2);
    void CalculateS0(Cut_Data *);
     
    void Initialize(std::vector<int> perm, VHAAG_res* ovl);

  public:

    VHAAG_res(int _nin,int _nout,int pn,VHAAG_res* ovl);

    ~VHAAG_res();

    void AddPoint(double Value);
    void GenerateWeight(ATOOLS::Vec4D *,Cut_Data *);
    void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *,double *);
    std::string Name() {return name;}
    void   MPISync()                 { p_vegas->MPISync(); }
    void   Optimize()                { p_vegas->Optimize(); }
    void   EndOptimize()             { p_vegas->EndOptimize(); }
    void   WriteOut(std::string pId) { if (m_ownvegas) p_vegas->WriteOut(pId); }
    void   ReadIn(std::string pId)   { if (m_ownvegas) p_vegas->ReadIn(pId); }

    int    RKF()                     { return m_rkf; }
    int    RD1()                     { return n_d1; }
    int    RD2()                     { return n_d2; }
    int    Type()                    { return m_type; }
    int    OType(); 
    std::map<int,Vegas*>* GetSharedVegasList() { return p_sharedvegaslist; }

    bool   OptimizationFinished()  { return p_vegas->Finished(); }
  };
}
#endif
