#ifndef VHAAG_h
#define VHAAG_h

#include "Single_Channel.H"
#include "Vegas.H"

namespace PHASIC {
  class VHAAG : public Single_Channel {
    int      nin,nout,n_p1;
    int      *p_perm,*p_mrep;
    double   m_s0;
    ATOOLS::Vec4D* m_q;
    Vegas* p_vegas;
    bool m_ownvegas;

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
    void CalculateS0(ATOOLS::Cut_Data *);
     
  public:
    VHAAG(int,int,int);
    ~VHAAG();

    void AddPoint(double Value);
    void GenerateWeight(ATOOLS::Vec4D *,ATOOLS::Cut_Data *);
    void GeneratePoint(ATOOLS::Vec4D *,ATOOLS::Cut_Data *,double *);
    std::string Name() {return name;}
    void   Optimize()                { p_vegas->Optimize(); }
    void   EndOptimize()             { p_vegas->EndOptimize(); }
    void   WriteOut(std::string pId) { p_vegas->WriteOut(pId); }
    void   ReadIn(std::string pId)   { p_vegas->ReadIn(pId); }
  };
}
#endif
