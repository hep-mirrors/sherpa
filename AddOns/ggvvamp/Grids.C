#include "Grids.H"
#include <string>
#include <fstream>

namespace Grids{
    
  Grids::Grids()
  {
    grid.resize(4);
  }

  void Grids::Init_Grids(int ngrids)
  {
    m_ngrids = ngrids;
    std::cout<<"Starting to initialize grids...\n";
    int index_map;
    std::string dir, file_name;
    dir = "/scratch/villani1/sherpa-db/AddOns/ggvvamp/"; 
      
    grid[0] = Grid(50,"/scratch/villani1/sherpa-db/AddOns/ggvvamp/new_grid_1.dat");
      grid[0].set_param(0.8,1-0.0001,0.5,1-0.01,0);
      grid[0].init_grid();
     
      
      grid[1] = Grid(50,"/scratch/villani1/sherpa-db/AddOns/ggvvamp/new_grid_2.dat");
      grid[1].set_param(0.8,1-0.0001,0.01,0.5,1);
      grid[1].init_grid();
      
      grid[2] = Grid(50,"/scratch/villani1/sherpa-db/AddOns/ggvvamp/new_grid_3.dat");
      grid[2].set_param(0.0001,0.8,0.5,1-0.01,2);
      grid[2].init_grid();
      
      grid[3] = Grid(50,"/scratch/villani1/sherpa-db/AddOns/ggvvamp/new_grid_4.dat");
      grid[3].set_param(0.0001,0.8,0.01,0.5,3);
      grid[3].init_grid();
    std::cout<<"Grid initialized. Time to initialize the interpolator...\n";

    Initialize_Interpolator();
  }

    void Grids::Initialize_Interpolator()
    {
       for(int i=0;i<4;i++)
       {
           Fill_Grid(i); 
       }
    }


    void Grids::Fill_Grid(int i)
    {
      int grid_site=0,index_map;
      std::ifstream grid_file(grid[i].strGetFileName());
      std::vector<std::vector<double>> values;
      values.resize(72);
      for(int j=0;j<72;j++) values[j].resize(50*50*50*50);

      std::cout<<"start reading\n";
      
      int grid_map;
      int coeff_map;
      double ma2,mb2,s,t;
      double as=grid[i].Get_as();
      double bs=grid[i].Get_bs();
      double at=grid[i].Get_at();
      double bt=grid[i].Get_bt();
      double beta, cos;
      
      for(int x1=0;x1<50;x1++)
        for(int x2=0;x2<50;x2++)
          for(int x3=0;x3<50;x3++)
            for(int x4=0;x4<50;x4++)
            { 
              grid_map=x1*50*50*50+x2*50*50+x3*50+x4;
              ma2 = grid[i].p2(x1/50.);
              mb2 = grid[i].p2(x2/50.);
              beta = as + x3*(bs-as)/50.;
              cos = 1 + 2*(at-bt)*x4/50. - 2*at;
              s = grid[i].s(beta,ma2,mb2);
              t = grid[i].t(cos,beta,s,ma2,mb2);
              for(int l=0;l<nLoop;++l)
                for(int hel=0;hel<nEhel;hel++)
                  for(int coeff=0;coeff<nEcoeff;coeff++)
                    for(int reim=0;reim<ReIm;reim++)
                    {
                      coeff_map=l*nEhel*nEcoeff*ReIm + hel*nEcoeff*ReIm + coeff*ReIm + reim;        
                      grid_file>>values[coeff_map][grid_map];
                      if(coeff>0&&coeff<5) values[coeff_map][grid_map]*=(pow((s+t+ma2+mb2),3));
                      else values[coeff_map][grid_map]*=(pow((s+t+ma2+mb2),2));   
                    }
            }
    
      std::cout<<"Reading completed\n";
      grid_file.close();

      std::cout<<"Building Interpolator\n";
      grid[i].build_interpolator(values);

    }
    
   int Grids::Interpolate(double s, double t, double ma2, double mb2)
   {
    int region(0);
    for(int i=0;i<4;i++)
      if(grid[i].check_point(s,t,ma2,mb2))
      {
          region=i;
          break;
      }
 

    std::vector<double> res;
    int index_map;
    double beta, cos;

    beta = grid[region].beta3(s,ma2,mb2);
    cos = grid[region].costheta3(s,t,ma2,mb2);
        
    grid[region].set_target(s,t,ma2,mb2);
           
    res = grid[region].evaluate(s,t,ma2,mb2);
    
    for(int l=0;l<nLoop;++l)
      for(int hel=0;hel<nEhel;hel++)
        for(int coeff=0;coeff<nEcoeff;coeff++)
          for(int reim=0;reim<ReIm;reim++)
          {
            index_map=l*nEhel*nEcoeff*ReIm + hel*nEcoeff*ReIm + coeff*ReIm + reim;        
            if(coeff>0&&coeff<5) E[hel][coeff+1][l][reim] = res[index_map]/(pow((s+t+ma2+mb2),3));
            else  E[hel][coeff+1][l][reim] = res[index_map]/(pow((s+t+ma2+mb2),2));
          }
    
    return region;
   }

}
