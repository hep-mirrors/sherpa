#ifndef ALPACA_Tools_Split_Variables_H
#define ALPACA_Tools_Split_Variables_H

#include "ALPACA/Tools/Parton.H"

#include <memory>
#include <vector>
#include <iostream>

namespace ALPACA {

  class Split_Variables {
  //Contains all information needed for splitting processes

  private :
    std::shared_ptr<Parton>   p_part_in;            //Incoming parton that can split
    std::shared_ptr<Parton>   p_part_recoil;        //Recoil parton

    double  m_tausplit;                             //Tau of splitting
    double  m_x;                                    //energy fraction x of splitting
    double  m_kT2;                                  //Sampled transverse momentum squared k_T^2
    double  m_mg2, m_mq2, m_Tstar;                  //Dynamic quantities screening mass and Tstar;
    double  m_f_oe;                                 //Overestimate of instantaneous splitting probability used

    ATOOLS::Vec4D     m_P_original, m_R_original;   //Original momentum
    ATOOLS::Vec4D     m_Pp, m_Kp, m_Rp;             //Outgoing momenta
    ATOOLS::Flavour   m_flav_out_1, m_flav_out_2;   //Outgoing flavours
    
  public:
    Split_Variables();
    Split_Variables(std::shared_ptr<Parton> ptr_part_in);
    ~Split_Variables();

    void SetOEVariables(std::shared_ptr<Parton> p_part_recoil, double tausplit, double x, double kT2, 
                        double mg2, double mq2, double Tstar, double f_oe,
                        ATOOLS::Vec4D P_original, ATOOLS::Vec4D R_original,
                        ATOOLS::Vec4D Pp, ATOOLS::Vec4D Kp, ATOOLS::Vec4D Rp);

    inline std::shared_ptr<Parton>        GetPartIn() const { return p_part_in; }
    inline std::shared_ptr<Parton>        GetPartRecoil() const { return p_part_recoil; }
    inline std::vector<ATOOLS::Vec4D>     GetOriginalMomentum() const { 
                                            static std::vector<ATOOLS::Vec4D> momenta = {m_P_original, m_R_original};
                                            return momenta;
                                          }
    inline std::vector<ATOOLS::Vec4D>     GetOutMomentum() const { 
                                            static std::vector<ATOOLS::Vec4D> momenta = {m_Pp, m_Kp, m_Rp};
                                            return momenta;
                                          }
    inline std::vector<double>            Getm2Tstar() const { 
                                            static std::vector<double> m2Tstar = {m_mg2, m_mq2, m_Tstar};
                                            return m2Tstar;
                                          }
    inline std::vector<ATOOLS::Flavour>   GetOutFlavours() const { 
                                            static std::vector<ATOOLS::Flavour> flav_out = {m_flav_out_1, m_flav_out_2};
                                            return flav_out;
                                          }
    inline double                         Getx() const{ return m_x; }
    inline double                         GetkT2() const{ return m_kT2; }
    inline double                         GetfOE() const{ return m_f_oe; }
  };
}

#endif
