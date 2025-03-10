#ifndef ALPACA_Tools_Parton_H
#define ALPACA_Tools_Parton_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Flow.H"

#include <memory>
#include <iostream>
#include <tuple>
#include <array>

namespace ALPACA {

  typedef std::pair<double,double> distpair;
  typedef std::pair<ATOOLS::Vec4D,ATOOLS::Vec4D> momentum_save;

  class Parton {
  private : 
    static long unsigned int   s_totalnumber;
    int                        m_number;
    double                     m_lambda, m_inittau;
    ATOOLS::Flavour            m_flav;         
    ATOOLS::Vec4D              m_momentum, m_oldmom; 
    ATOOLS::Vec4D              m_initpos; 
    ATOOLS::Flow               m_flow;

    int                        m_timekeeper;      //Time evolution method (see 2306.12131). 
                                                  //Set the corresponding for each initialized parton:
                                                  //    0 - Original, no timekeeper. DEFAULT.
                                                  //    1 - m_ahat=(1,0,0,0) in the lab frame. 
                                                  //    2 - m_ahat=(1,0,0,0) in the global COM frame
    ATOOLS::Vec4D              m_ahat;            //Used if timekeeper = 1 or 2

    std::shared_ptr<Parton>    p_split_recoil_partner;
    std::shared_ptr<Parton>    p_last_scatter_1;  //Pointer to the last particle that this particle scattered with
    std::shared_ptr<Parton>    p_last_scatter_2;  //Pointer to the second to last particle that this particle scattered with
    std::shared_ptr<Parton>    p_last_scatter_3;  //Pointer to the third to last particle that this particle scattered with
    std::shared_ptr<Parton>    p_split_partner;     //If the particle is created from a splitting: pointer to the other particle resulting from that split
    std::pair<double,double>   m_formation_tau;   //[start tau of formation for merge, end tau of formation for merge] used when creating or merged particle
    std::pair<double,double>   m_formation_time;   //[start tau of formation for split, end tau of formation for split] used when finding tau_split
    bool                       m_splitted_merged;        //Bool to tell if the particle has already splitted (and if so has been removed from tau lists, 
                                                  // interesting in the case of the particle being a recoil partner for another splitting process,
                                                  // since there can still exist a shared_ptr to it as a recoil partner)

    int                        m_N_scatter; //Number of times the parton has elastically scattered
    std::vector<std::pair<int, std::pair<momentum_save, ATOOLS::Flavour>>>  m_x_p_history;  //Save previous momenta and flavour after each type of event
    std::vector<double>        m_x_p_tau; //Temporary list to keep track of corresponding taus in XP history
    std::array<double,4>       a_momentum;
    std::array<double,4>       a_initpos;


  public:
    Parton();
    Parton(ATOOLS::Flavour, ATOOLS::Vec4D);
    Parton(ATOOLS::Flavour, ATOOLS::Vec4D, ATOOLS::Vec4D, double);
    Parton(Parton * part);
    ~Parton();

    inline void   SetMomentum(const ATOOLS::Vec4D & mom) { 
      m_momentum = mom;
      a_momentum = {mom[1], mom[2], mom[3], mom[0]};
      if(m_timekeeper == 1 || m_timekeeper == 2){
        m_lambda = 1./(2.*(m_ahat*m_momentum));
      }
    }
    inline void   SetPosition(const ATOOLS::Vec4D & pos) {
      m_initpos = pos;
      a_initpos = {pos[1], pos[2], pos[3], pos[0]};
    }
    inline void   SetTimekeeper(const int timekeeper, const ATOOLS::Vec4D ahat) {m_timekeeper = timekeeper; m_ahat = ahat; m_lambda = 1./(2.*(ahat*m_momentum)); }
    inline void   SetFlav(const ATOOLS::Flavour & flav) {m_flav = flav;}
    inline void   SetFlow(ATOOLS::Flow flow) {m_flow = flow;}
    inline void   SetFlow(const int index,const int code=-1) {m_flow.SetCode(index,code);}
    inline void   SetInittau(const double itau) {m_inittau = itau;}
    inline void   SetLambda(const double lambda) {m_lambda = lambda;}
    inline void   SetSplitMergeStatus(const bool split_merge_status) { m_splitted_merged = split_merge_status; }
    inline void   SetSplitPartner(std::shared_ptr<Parton> split_partner) { p_split_partner = split_partner; }
    inline void   SetSplitRecoilPartner(std::shared_ptr<Parton> split_recoil_partner) { p_split_recoil_partner = split_recoil_partner; }
    void   SetLastScatter(std::shared_ptr<Parton> part);
    void   SetNumber();
    void   SetFormationTime(double t, double formation_t);

    void                                                                                    AddXPHistory(int process_type, ATOOLS::Vec4D x, ATOOLS::Vec4D p, ATOOLS::Flavour flav, double tau);
    inline std::vector<std::pair<int, std::pair<momentum_save, ATOOLS::Flavour>>>           GetAllXPHistory() {return m_x_p_history; }
    double                                                                                  GetLastXPt();
    inline std::pair<int, std::pair<momentum_save, ATOOLS::Flavour>>                        GetFirstXP() {return m_x_p_history.back(); }
    inline std::pair<int, std::pair<momentum_save, ATOOLS::Flavour>>                        GetLastXP() {return m_x_p_history[0]; }
    std::pair<bool, std::pair<std::pair<ATOOLS::Vec4D, ATOOLS::Vec4D>, ATOOLS::Flavour>>    GetXP(double t);
    std::pair<bool, std::pair<std::pair<ATOOLS::Vec4D, ATOOLS::Vec4D>, ATOOLS::Flavour>>    GetXPTau(double tau);


    inline const ATOOLS::Vec4D     Momentum() const {return m_momentum;}
    inline const std::array<double,4> aMomentum() const {return a_momentum; }
    inline const ATOOLS::Vec4D     Position() const {return m_initpos;}
    inline const std::array<double,4> aPosition() const {return a_initpos; }
    inline const double            Inittau() const {return m_inittau;}
    inline double                  E() const {return m_momentum[0];}
    inline ATOOLS::Flavour         Flav() const {return m_flav;}
    inline ATOOLS::Flow            GetFlow() const {return m_flow;}
    inline int                     GetFlow(const int index) const {return m_flow.Code(index);}
    inline int                     Number() const { return m_number; }
    inline double                  GetLambda() const { return m_lambda; }
    inline std::shared_ptr<Parton> GetLastScatter() {return p_last_scatter_1;}
    inline std::shared_ptr<Parton> GetSecondLastScatter() {return p_last_scatter_2;}
    inline std::shared_ptr<Parton> GetThirdLastScatter() {return p_last_scatter_3;}
    inline std::pair<double,double> GetFormationTau() {return m_formation_tau;}
    inline std::pair<double,double> GetFormationTime() {return m_formation_time;}
    inline bool                     GetSplitMergeStatus() {return m_splitted_merged; }
    inline std::shared_ptr<Parton>  GetSplitPartner() {return p_split_partner; }
    inline std::shared_ptr<Parton>  GetSplitRecoilPartner() {return p_split_recoil_partner; }
    inline int                      GetNScatter() {return m_N_scatter; }
    inline void                     AddNScatter() {m_N_scatter = m_N_scatter + 1;}
    inline static void        Reset(const int number=0)    { s_totalnumber = number; }
    inline void               SaveMomentum() {m_oldmom = m_momentum; }

    ATOOLS::Vec4D             Position(const double tau) const;
    ATOOLS::Vec4D             PositionFromt(const double t) const;
    double                    Tau(const double t) const;
    double                    Dist2(const std::shared_ptr<Parton> part, const double tau);
    bool                      operator==(Parton part);
    double                    Abs3Momentum();

    double                    TauBound(size_t i, double bound);

    

    

    double TauBar(std::shared_ptr<Parton> part);
	  double TauBar(std::shared_ptr<Parton> part, std::array<double,3> shift);

    std::array<double,4> getPosition(double t);
    std::array<double,4> getPositionNoBound(double t);
    std::array<double,4> getPositionTau(double t);
    std::array<double,4> getPositionTauNoBound(double t);

    double getInvDistTau(std::shared_ptr<Parton> part, double tau, std::array<double,3> shift);

    std::array<double,3> findStartPoint(std::shared_ptr<Parton> part, double taunow);
    std::pair<std::pair<double,double>,std::array<double,3>> findTaubarDMin(std::shared_ptr<Parton> part, double taumin, double taumax, double dminmax, std::array<double,3> & startp, bool print);
    std::pair<double,double> findDMinPos2D(std::shared_ptr<Parton> part, double m);
    void findIntersection(std::shared_ptr<Parton> part, int & side, double & k, double & l, double & m, std::array<double,3> shat);

    //Temporary math for taubar, should be removed when switching all relevant variables to ATOOLS::Vec4D
    template <typename T>
    inline T sqr(T value){
      return value * value;
    }

    inline std::array<double,4> subtract(std::array<double,4> a, std::array<double,4> b){
      std::array<double,4> res;
      for (size_t i = 0; i < 4; ++i) {
        res[i] = a[i] - b[i];
      }
      return res;
    }

    inline std::array<double,4> add(std::array<double,4> a, std::array<double,4> b){
      std::array<double,4> res;
      for (size_t i = 0; i < 4; ++i) {
        res[i] = a[i] + b[i];
      }
      return res;
    }

    inline std::array<double,4> scaltimes(double a, std::array<double,4> b){
      std::array<double,4> res;
      for (size_t i = 0; i < 4; ++i) {
        res[i] = a*b[i];
      }
      return res;
    }

    inline double abs2(std::array<double,4> x){
      return x[3]*x[3] - x[0]*x[0] - x[1]*x[1] - x[2]*x[2];
    }

    inline double abs_vec(std::array<double,4> x) {
      return sqrt(x[3]*x[3] - x[0]*x[0] - x[1]*x[1] - x[2]*x[2]);
    }

    inline double dotpr(std::array<double,4> x, std::array<double,4> y){
      return x[3]*y[3] - x[0]*y[0] - x[1]*y[1] - x[2]*y[2];
    }


    inline std::array<double,3> subtract3(std::array<double,3> a, std::array<double,3> b){
      std::array<double,3> res;
      for (size_t i = 0; i < 3; ++i) {
        res[i] = a[i] - b[i];
      }
      return res;
    }

    inline std::array<double,3> add3(std::array<double,3> a, std::array<double,3> b){
      std::array<double,3> res;
      for (size_t i = 0; i < 3; ++i) {
        res[i] = a[i] + b[i];
      }
      return res;
    }

    inline std::array<double,3> scaltimes3(double a, std::array<double,3> b){
      std::array<double,3> res;
      for (size_t i = 0; i < 3; ++i) {
        res[i] = a*b[i];
      }
      return res;
    }

    inline double abs32(std::array<double,3> x){
      return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    }

    inline double abs3(std::array<double,3> x){
      return sqrt(abs32(x));
    }

    inline double dotpr3(std::array<double,3> x, std::array<double,3> y){
      return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
    }

    inline double abs342(std::array<double,4> x){
      return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    }
  };

  std::ostream& operator<<(std::ostream &,const Parton &);
  
}

#endif
