#include <iostream>
#include "Apacic.H"

using namespace APACIC;

int main() {  
  bool test=0; // testmode

//   std::cout<<"Run Apacic in Testmode ? (Yes = 1, No = 0)"<<std::endl;
//   std::cin>>test;

  Apacic Generator(test);
  Generator.Init();
  if (!test) {
    Generator.CrossSections();
    Generator.GenerateEvents();
  }
  return 1;
}


