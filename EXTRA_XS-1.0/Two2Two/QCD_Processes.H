#ifndef QCD_Process_H
#define QCD_Process_H

#include "XS_Group.H"

namespace EXTRAXS {

  class QCD_Processes: public XS_Group {
    double aS;
    double ycut;
    int    jetfinder;

    void Fill4qmodes();
    void Fill2q2gmodes();
    void Fill4gmodes();

    void CreateSelector();
  public:
    QCD_Processes();
    double Scale(AMATOOLS::vec4d *);
    double KFactor(double);
  };
}

#endif
