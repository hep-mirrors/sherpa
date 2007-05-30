#ifndef QCD_Vertices_H
#define QCD_Vertices_H

#include "QCD_Currents.H"

namespace EXTRAXS {

  class QCD_GGG: public Vertex {
  private:

    double m_cpl;

  public:

    QCD_GGG(Current_Base *const c);

    void Evaluate();

  };// end of class QCD_GGG

  class QCD_GGT: public Vertex {
  private:

    double m_cpl;

  public:

    QCD_GGT(Current_Base *const c);

    void Evaluate();

  };// end of class QCD_GGT

  class QCD_GTG: public Vertex {
  private:

    double m_cpl;

  public:

    QCD_GTG(Current_Base *const c);

    void Evaluate();

  };// end of class QCD_GTG

  class QCD_QQG: public Vertex {
  private:

    double m_cpl;

  public:

    QCD_QQG(Current_Base *const c);

    void Evaluate();

  };// end of class QCD_QQG

  class QCD_QGQ: public Vertex {
  private:

    double m_cpl;

  public:

    QCD_QGQ(Current_Base *const c);

    void Evaluate();

  };// end of class QCD_QGQ

}// end of namespace EXTRAXS

#endif
