#ifndef QCD_Process_C_H
#define QCD_Process_C_H

#include "XS_Group.H"
#include "Channel_Interface.H"

namespace EXTRAXS {

  class QCD_Processes_C: public XS_Group {
  public:
    enum Mode { gggg         = 100,
		qqbgg        = 200,
		ggqqb        = 300,
		qgqg         = 400,
		q1q1q1q1     = 510,
		q1q2q1q2     = 520,
		q1q2bq1q2b   = 530,
		q1q1bq1q1b   = 600,
		q1q1bq2q2b   = 700,
		All          =   1,
		None         =   0,
		Unknown      = 999 };
  private:
    int m_fsrmode;
    PHASIC::Channel_Interface *p_fsrinterface;
    double aS;
  public:
    // constructor
    QCD_Processes_C(PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *,
		  ATOOLS::Flavour *,ATOOLS::Selector_Data *,int,int,double,bool=true);
    // setters
    inline void SetFSRMode(int _m_fsrmode) { m_fsrmode=_m_fsrmode; }
    inline void SetFSRInterface(PHASIC::Channel_Interface *_p_fsrinterface)
    { p_fsrinterface=_p_fsrinterface; }
    // getters
    inline int                 FSRMode()      { return m_fsrmode;      }
    PHASIC::Channel_Interface *FSRInterface() { return p_fsrinterface; }
    // member functions
    void FillMode(Mode);
    void CreateFSRChannels();
  };
}

#endif
