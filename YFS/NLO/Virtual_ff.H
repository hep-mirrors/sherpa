#ifndef YFS_Virtual_Virtual_ff_H
#define YFS_Virtual_Virtual_ff_H


namespace YFS{

    class Virtual_ff{
    public:
        Virtual_ff(int order);
        ~Virtual_ff();
        bool CheckOrder();
        double Calculate(double beta, int order);
        double Del1();
        double Del2();
        double Del3();
        int m_order, m_orderMin, m_orderMax;
        double m_gamma;
        double m_beta;
        double m_virtual;
        long double m_del1, m_del2, m_del3;
        inline double GetVirtual() {return m_virtual ;}
    };
}

#endif
  /*!
    \file Virtual_ff.H
    \brief contains the class Virtual_ff
  */

  /*!
    \class Virtual_ff
    \brief This class calculates the pure virtual corrections to $e^+e^- \to \mu^+\mu^-$ according to 
    https://arxiv.org/pdf/hep-ph/0006359.pdf
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Virtual_ff
  ////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
    \var int Virtual_ff:m_order 
    \brief This is the pragmatic order to which the virtual corrections are caclulated

    \var double Virtual_ff::m_gamma
    \brief $\gamma = \frac{2\alpha}{\pi}*\left(L-1\right)$ which should be supplied by the YFS_handler

*/