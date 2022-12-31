#include <iostream>
#include <cmath>

#include <btwxt.h>
#include <stdlib.h> 
#include <fstream>

using namespace Btwxt;


class Grid{ 
  int m_N;
  double m_delta_x;
  double m_as;
  double m_bs;
  double m_at;
  double m_bt;
  int m_region;
  double m_mass_power;
  std::string m_file_name;
  
  std::vector<double> first_axis;
  std::vector<double> second_axis;
  std::vector<double> third_axis;
  std::vector<double> fourth_axis;
  std::vector<double> target={0.,0.,0.,0.}; 

  RegularGridInterpolator *my_interpolator;   

  public:

  Grid():my_interpolator(NULL){};

  Grid(const int N, const std::string file_name)
      : m_N(N),my_interpolator(NULL)
  {
    m_delta_x = 1./N;
    m_file_name = file_name;
    first_axis.resize(N);
    second_axis.resize(N);
    third_axis.resize(N);
    fourth_axis.resize(N);
  };
  
  ~Grid()
  {
    delete my_interpolator;
  };

  
  void set_param(double as, double bs, double at, double bt, int i)
  {
    m_as = as;
    m_bs = bs;
    m_at = at;
    m_bt = bt;
    m_region = i;
  };

  
  
  double k(double s, double p32, double p42)
  {
    return sqrt( s*s + p32*p32 + p42*p42 - 2.*( s*p32 + p32*p42 + s*p42) );
  };
  
  
  
  double xp(double p2)
  {
    return pow((p2-40.)/(1690000.-40.),0.2);
  };
  
  
  
  double p2(double x)
  {
    return (1690000-40)*pow(x,5)+40;
  };
  
  
  
  double beta3(double s, double p32, double p42)
  {
    return k(s,p32,p42)/(s + p32 - p42);
  };
  


  double costheta3(double s, double t, double p32, double p42)
  {
    return (2*t + s - p32 - p42)/k(s,p32,p42);
  };


  
  double s(double beta, double p32, double p42)
  {
    return ( p32*(1+beta*beta) + p42*(1-beta*beta) + 2*sqrt( p32*p42*(1-beta*beta) + p32*p32*beta*beta ) )/(1-beta*beta);
  };
  
  
  
  double t(double cos, double beta, double s, double p32, double p42)
  {
    return 0.5*(s*(beta*cos-1)+p32*(beta*cos+1)+p42*(1-beta*cos));
  };

  void init_grid()
  {
    for(int i=0;i<m_N;i++)
    {
      first_axis[i] = i*m_delta_x;
      second_axis[i] = i*m_delta_x;
      third_axis[i] = i*m_delta_x;
      fourth_axis[i] = i*m_delta_x;
    }
  };
  

  bool check_grid_site(double ma2, double mb2, double beta, double cos)
  {
    double int_1,int_2,int_3,int_4;
    double rest_1 = std::modf(xp(ma2)*m_N,&int_1);
    double rest_2 = std::modf(xp(mb2)*m_N,&int_2);
    double rest_3 = std::modf(((m_as-beta)/(m_as-m_bs))*m_N,&int_3);
    double rest_4 = std::modf((0.5*(cos - 1 + 2*m_at)/(m_at-m_bt))*m_N,&int_4);

    std::cout<<xp(ma2)<<' '<<xp(mb2)<<' '<<(m_as-beta)/(m_as-m_bs)<<' '<<0.5*(cos - 1 + 2*m_at)/(m_at-m_bt)*m_N<<std::endl;
    std::cout<<int_1<<' '<<int_2<<' '<<int_3<<' '<<int_4<<std::endl;
    std::cout<<rest_1<<' '<<rest_2<<' '<<rest_3<<' '<<rest_4<<std::endl;
    
    
    return ((rest_1<1.e-12)&&(rest_2<1.e-12)&&(rest_3<1.e-12)&&(rest_4<1.e-12));    
  };
 
  void build_interpolator(std::vector<std::vector<double>> vals)
  {
    std::vector<GridAxis> my_grid{first_axis, second_axis, third_axis, fourth_axis};
  
    GriddedData gridded_data(my_grid,vals);

    gridded_data.set_axis_interp_method(0, Method::CUBIC);
    gridded_data.set_axis_interp_method(1, Method::CUBIC);
    gridded_data.set_axis_interp_method(2, Method::CUBIC);
    gridded_data.set_axis_interp_method(3, Method::CUBIC);
  
    gridded_data.set_axis_extrap_method(0, Method::CUBIC);
    gridded_data.set_axis_extrap_method(1, Method::CUBIC);
    gridded_data.set_axis_extrap_method(2, Method::CUBIC);
    gridded_data.set_axis_extrap_method(3, Method::CUBIC);
  
    my_interpolator = new RegularGridInterpolator(gridded_data);   
  };


  bool check_point(double s_inv,double t_inv,double p32,double p42)
  {
    double p3,p4,cos,beta;
    double x1,x2,x3,x4;
  


    p3=sqrt(p32);
    p4=sqrt(p42);

    cos = costheta3(s_inv,t_inv,p32,p42);
    beta = beta3(s_inv,p32,p42);

    x1 = xp(p32);
    x2 = xp(p42);
    x3 = (m_as-beta)/(m_as-m_bs); 
    x4 = 0.5*(cos - 1 + 2*m_at)/(m_at-m_bt);
    
    if(m_region==0) return (beta>=0.8&&cos<=0);
    if(m_region==1) return (beta>=0.8&&cos>0);
    if(m_region==2) return (beta<0.8&&cos<=0);
    if(m_region==3) return (beta<0.8&&cos>0);
   
    return false;
  };


  void set_target(double s_inv, double t_inv, double p32, double p42)
  {
    double p3,p4,cos;
    double x1,x2,x3,x4;

    p3=sqrt(p32);
    p4=sqrt(p42);

    cos = costheta3(s_inv,t_inv,p32,p42);
  
    x1 = xp(p32);
    x2 = xp(p42);
    x3 = (m_as-beta3(s_inv,p32,p42))/(m_as-m_bs); 
    x4 = 0.5*(cos - 1 + 2*m_at)/(m_at-m_bt);
    
    target[0]=x1;
    target[1]=x2;
    target[2]=x3;
    target[3]=x4;
  };

  std::vector<double> evaluate(double s_inv,double t_inv,double p32,double p42)
  { 
    std::vector<double> result;
    result = (*my_interpolator)(target); 
    return result;
  };

  double evaluate(double s_inv,double t_inv,double p32,double p42, double mass_power)
  { 
    std::vector<double> result;
    result = (*my_interpolator)(target);
    result[0]=result[0]/pow((p32+p42+s_inv+t_inv),mass_power); 
    return result[0];
  };

  void SetMassPower(double mass_power) {m_mass_power=mass_power;}
  std::string strGetFileName() {return m_file_name;}

  double Get_as(){return m_as;};
  double Get_bs(){return m_bs;};
  double Get_at(){return m_at;};
  double Get_bt(){return m_bt;};


};
