#ifndef Hadron_List_H
#define Hadron_List_H

#include "Cluster.H"
#include "Flavour.H"
#include <map>
#include <string>

namespace AHADIC {
  class ConstituentCharacteristic {
  public:
    double m_mass, m_flweight, m_spweight, m_smearing;
    int    m_ispin;
    ConstituentCharacteristic(double _m=0., int _is=-1, 
			      double _fl=0., double _sp=0., double _sm=0) : 
      m_mass(_m), m_flweight(_fl), m_spweight(_sp), m_ispin(_is), 
      m_smearing(_sm/_m) {}
    double Mass()      { return m_mass;     }
    double FlWeight()  { return m_flweight; }
    double SpWeight()  { return m_spweight; }
    double TotWeight() { return m_spweight*m_flweight; }
    double Smearing()  { return m_smearing; }
    int    ISpin()     { return m_ispin;    }
    void   SetMass(double _mass)         { m_mass     = _mass;     }
    void   SetFlWeight(double _flweight) { m_flweight = _flweight; }
    void   SetSpWeight(double _spweight) { m_spweight = _spweight; }
    void   SetSmearing(double _smearing) { m_smearing = _smearing; }
    void   SetISpin(int _ispin)          { m_ispin    = _ispin;    }
  };

  typedef std::map<ATOOLS::Flavour,ConstituentCharacteristic *> FlavCCMap;
  typedef FlavCCMap::iterator FlavCCMap_Iterator;

  class ConstituentMasses {
  public:
    FlavCCMap CCMap;
    ConstituentMasses();
    ~ConstituentMasses();

    double Mass(ATOOLS::Flavour & flav);
    double FlWeight(ATOOLS::Flavour & flav);
    double SpWeight(ATOOLS::Flavour & flav);
    double TotWeight(ATOOLS::Flavour & flav);
    double Smearing(ATOOLS::Flavour & flav);
    int    ISpin(ATOOLS::Flavour & flav);

    void   PrintConstituents();
  };

  //##############################################################################
  //##############################################################################
  //##############################################################################
  //##############################################################################
  //##############################################################################

  typedef std::map<FlavourPair,double>               WFcomponent;
  typedef WFcomponent::iterator                      WFcompiter;
  class Wave_Function {
  private:
    ATOOLS::Flavour m_hadron;
    int             m_spin, m_kfcode;
    WFcomponent     m_waves;
    bool            m_barrable;
  public:
    Wave_Function();
    ~Wave_Function() { m_waves.clear(); }
    Wave_Function(const ATOOLS::Flavour &);    

    ATOOLS::Flavour GetFlavour() const  { return m_hadron;    }
    void SetSpin(const int _spin)       { m_spin   = _spin;   }
    int  Spin() const                   { return m_spin;      }
    void SetKfCode(const int _kfcode)   { m_kfcode = _kfcode; }
    int  KfCode() const                 { return m_kfcode;    }

    WFcomponent * GetWaves()            { return &m_waves;    }
    void AddToWaves(FlavourPair,double);

    bool ExistAnti()                    { return m_barrable;  }
    Wave_Function * Anti();
    friend std::ostream & operator<<(std::ostream & s, Wave_Function & wf);
  };



  struct Hadron_Multiplet {
    std::string        m_name;
    ATOOLS::FlavourSet m_elements;
    double             m_weight;
  };


  typedef std::map<ATOOLS::Flavour,Wave_Function *> WFMap; 
  typedef WFMap::iterator                           WFmapiter;
  typedef std::map<int,Hadron_Multiplet *>          MultipletMap;
  typedef MultipletMap::iterator                    Mpletmapiter;

  class Hadron_List {
  private:
    WFMap             * p_wavefunctions;
    MultipletMap      * p_multiplets;

    void ConstructWaveFunctions();
    void ConstructAntiWaveFunctions();
    void CreateMultiplets();
    void LookUpAngles(const int angular,const int spin,double & costh,double & sinth);
    bool LookUpWaveFunction(const int kfcode,const int fl1,const int fl2,
			    const Wave_Function * wf);
  public:
    Hadron_List();
    ~Hadron_List();
    WFMap * GetWaveFunctions() { return p_wavefunctions; }
    void    PrintWaveFunctions();
  };


  //##############################################################################
  //##############################################################################
  //##############################################################################
  //##############################################################################
  //##############################################################################


  typedef std::map<std::string,double> StringDoubleMap;
  typedef StringDoubleMap::iterator    SDiter;

  struct CChannel_Element {
    FlavourPair     m_pairs[2];
    ATOOLS::Flavour m_hadrons[2];
    double          m_mass, m_weight;
  };
  std::ostream & operator<<(std::ostream & s, const CChannel_Element & ce);

  class Channel_Sorting_Mass {
  public :
    bool operator() (const CChannel_Element * ch1,const CChannel_Element * ch2) {
      if (ch1->m_mass>ch2->m_mass) return 0;
      return 1;
    }
  };

  typedef std::set<CChannel_Element *,Channel_Sorting_Mass> ChannelSet;
  typedef ChannelSet::iterator                              Channeliter;
  std::ostream & operator<<(std::ostream & s, const ChannelSet & cs);


  class Flavour_Sorting_Mass {
  public :
    bool operator() (const ATOOLS::Flavour & fl1,const ATOOLS::Flavour & fl2) {
      if (fl1.Mass()<fl2.Mass()) return 0;
      return 1;
    }
  };
  typedef std::map<ATOOLS::Flavour,double,Flavour_Sorting_Mass>  WeightMap;
  typedef WeightMap::iterator                                    Weightiter;
  typedef std::map<FlavourPair,WeightMap *>                      AllWeights;           
  typedef AllWeights::iterator                                   Allweightiter; 
  typedef std::map<FlavourPair,ChannelSet *>                     AllChannels;
  typedef AllChannels::iterator                                  Allchanneliter; 


  class Hadron_Parameters {
  private:
    ConstituentMasses * p_constituents;
    Hadron_List       * p_hadrons;
    AllWeights        * p_allweights;
    AllChannels       * p_allchannels;

    StringDoubleMap m_parametermap;
    SDiter          m_piter;


    void   ReadParameters(std::string,std::string);
    void   CreateWeightLists();
    void   AddToWeights(FlavourPair,ATOOLS::Flavour,double);
    void   CreateChannelLists();
  public:
    Hadron_Parameters();
    ~Hadron_Parameters();

    void   Init(std::string,std::string);
    double Get(std::string keyword);
    double MinimalMass(FlavourPair *);

    double Mass(ATOOLS::Flavour & flav)   { return p_constituents->Mass(flav); }         
    AllWeights        * GetAllWeights()   { return p_allweights;   }
    AllChannels       * GetAllChannels()  { return p_allchannels;  }
    ConstituentMasses * GetConstituents() { return p_constituents; }

    void   PrintWeightLists();
    void   PrintChannelLists();
  };

  extern Hadron_Parameters hadpars;

}
#endif
