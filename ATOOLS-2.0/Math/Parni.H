#ifndef Parni_H
#define Parni_H

#include <list>

namespace ATOOLS {
  template <int t_dim>
  class N_Tuple {
    double m_x[t_dim];
  public:
    N_Tuple();
    N_Tuple(const N_Tuple<t_dim>&);
    inline double& operator[] (int i);
    inline double operator[] (int i) const;
    N_Tuple<t_dim> & operator=(const N_Tuple<t_dim> &);
  };

  template <int t_dim>
  double& N_Tuple<t_dim>::operator[](int i)
  {
    return m_x[i];
  }

  template <int t_dim>
  double N_Tuple<t_dim>::operator[](int i) const
  {
    return m_x[i];
  }


  template <int t_dim>
  bool operator<(const N_Tuple<t_dim> &, const N_Tuple<t_dim> &);

  template <int t_dim>
  bool operator<=(const N_Tuple<t_dim> &, const N_Tuple<t_dim> &);

  template <int t_dim>
  bool operator==(const N_Tuple<t_dim> &, const N_Tuple<t_dim> &);


  template <int t_dim>
  class Parni_Channel {
    int m_id;

    N_Tuple<t_dim> m_lower;
    N_Tuple<t_dim> m_upper;

    unsigned m_points;
    double m_alpha;
    double m_sum;
    double m_sum2;
    double m_min;
    double m_max;
    N_Tuple<t_dim> m_point;
    bool m_include;
  public:
    Parni_Channel();
    Parni_Channel(const Parni_Channel<t_dim> &);
    Parni_Channel & operator=(const Parni_Channel<t_dim> &);
    bool Includes(const N_Tuple<t_dim> &);
    bool Includes() const { return m_include; };
    N_Tuple<t_dim> & GeneratePoint();
    void AddPoint(double,double);

    int Id() const { return m_id; }
    void SetId(int id) { m_id=id; }
    const N_Tuple<t_dim> & Lower() const { return m_lower; }
    const N_Tuple<t_dim> & Upper() const { return m_upper; }
    void SetLower(const N_Tuple<t_dim> & low) { m_lower=low; }
    void SetUpper(const N_Tuple<t_dim> & up)  { m_upper=up; }
    
    double Alpha() const { return m_alpha;}
    double NPoints() const { return m_points; }
    void SetAlpha(const double a) { m_alpha=a;}
    double OptimizeAlpha();
    void Reset();
    double Integral() const;
  };

  template <int t_dim>
  bool operator==(const Parni_Channel<t_dim> &, const Parni_Channel<t_dim> &);

  template <int t_dim>
  std::ostream & operator<<(std::ostream &, const Parni_Channel<t_dim> &);

  template <int t_dim>
  class Parni {
    typedef std::list<Parni_Channel<t_dim> > Channel_List;
    typedef typename Channel_List::iterator Channel_Iterator;
    typedef typename Channel_List::const_iterator Const_Channel_Iterator;

    std::string m_name;
    int m_mode;
    Channel_List m_channels;
    Channel_Iterator m_selected;
    unsigned long m_nevt;
    double m_weight;

    void SplittChannel();
    void SplittChannel(Channel_Iterator);
    void SplittChannel(const Parni_Channel<t_dim> &,int);
    void Merge(int =12);
    N_Tuple<t_dim> MinDelta() const;
  public:
    Parni(const std::string &,int);
    N_Tuple<t_dim> & GeneratePoint();
    double GenerateWeight(const N_Tuple<t_dim> &);
    void AddPoint(double);
    void AddPoint(const N_Tuple<t_dim> &, double);

    void Optimize();

    void WriteOut(const std::string &);
    void WriteOutHisto(const std::string &);
    void ReadIn(const std::string &);
  };

}

#endif
