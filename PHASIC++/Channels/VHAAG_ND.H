#ifndef PHASIC_Channels_VHAAG_ND_h
#define PHASIC_Channels_VHAAG_ND_h

#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Vegas.H"

namespace PHASIC {
  class VHAAG_ND : public Single_Channel {
    int      n_p1,m_type;
    int      *p_perm,*p_mrep;
    double   m_s0;
    ATOOLS::Vec4D* m_q;
    Vegas* p_vegas;
    bool m_ownvegas;
    
    Vegas** p_sharedvegaslist;

    double PiFunc(double a1,double a2,
		  double s1b,double s2b,double c);
    void Split(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
	       ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int,int,double *ran);
    void SplitF(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
		ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int,double *ran);
    void Split0(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
		ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int,double *ran);
    void Split1(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
		ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int,double *ran);
    void SingleSplit(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D Q,
		     ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,int,double *ran);
    void SingleSplitF(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D Q,
		      ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double *ran);
    void SingleSplitF0(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
		       ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2,double *ran);
    double SplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
		       ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int,int,double *ran);
    double SplitFWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int,double *ran);
    double Split0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int,double *ran);
    double Split1Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int,double *ran);
    double SingleSplitWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D& Q,
			     ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,int,double *ran);
    double SingleSplitFWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D& Q,
			      ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double *ran);
    double SingleSplitF0Weight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			      ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double *ran);
    void GenerateBranch(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D Q,
			ATOOLS::Vec4D* q,int n,double *ran);
    double BranchWeight(ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,ATOOLS::Vec4D &Q,
			ATOOLS::Vec4D* q,int n,double *ran);
    void ConstructMomenta(double a1,double a2,
			  double s1,double s2,double s,
			  ATOOLS::Vec4D q1,ATOOLS::Vec4D q2,
			  ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2);
    void ConstructMomenta(double a1,double phi,
			  double s1,double s2,double s,
			  ATOOLS::Vec4D q1,ATOOLS::Vec4D& p1,ATOOLS::Vec4D& p2);
    void CalculateS0(Cut_Data *);
     
    void Initialize(int _nin,int _nout,std::vector<int> perm, VHAAG_ND* ovl);

  public:

    VHAAG_ND(int _nin,int _nout,int pn, VHAAG_ND* ovl);
    VHAAG_ND(int _nin,int _nout,std::vector<size_t> tp, VHAAG_ND* ovl);

    ~VHAAG_ND();

    void AddPoint(double Value);
    void GenerateWeight(ATOOLS::Vec4D *,Cut_Data *);
    void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *,double *);
    std::string Name() {return name;}
    void   MPISync();
    void   Optimize();
    void   EndOptimize()             { p_vegas->EndOptimize(); }
    void   WriteOut(std::string pId) { if (m_ownvegas) p_vegas->WriteOut(pId); }
    void   ReadIn(std::string pId)   { if (m_ownvegas) p_vegas->ReadIn(pId); }

    int    OType(); 
    Vegas** GetSharedVegasList() { return p_sharedvegaslist; }

    bool   OptimizationFinished()  { return p_vegas->Finished(); }
  };
}
#endif
