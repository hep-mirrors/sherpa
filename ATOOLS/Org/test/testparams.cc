#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <string>
#include <vector>
#include <cmath>
#include "ATOOLS/Org/SherpaParameter.H"
#include "ATOOLS/Org/SherpaSwitch.H"
#include "ATOOLS/Org/SherpaVector.H"
#include "ATOOLS/Org/SherpaMatrix.H"

using std::vector;
int main() {
    
    SherpaParameter <double>p("PARAM", 1.0, 1.0, 0.0, 4.0);

    vector<int> k;
    k+=1,2,3,4,5,6,7;
    SherpaSwitch <int> s("SWITCH", 1,2, k);
    s.print();

    vector<std::string> v1, v2;
    v1+="a","b","c","d";
    v2+="a","b","e","d";
    SherpaVector <std::string>v("VECTOR", v1, v2);
    v.print();
    p.lock();
    p.print();
    p.setValue(90);

    std::vector< std::vector<std::string> >m1, m2;
    m1+=v1,v2;
    m2+=v1,v1;

    SherpaMatrix <std::string >m("MATRIX", m1, m2);
    m.print();


    return 0;
}
