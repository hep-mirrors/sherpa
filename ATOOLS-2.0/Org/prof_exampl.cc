#include <iostream>
#include <fstream>
#include "prof.hh"

unsigned int n = 0;

void f0()
{
    PROFILE_HERE;
    ++n;
}
template <int b>
void f()
{
    PROFILE_HERE;
    f<b-1>();
    f<b-2>();
    f<b-1>();
    f<b-2>();
    f<b-1>();
    f<b-2>();
    f0();
}


template <>
void f<1>()
{
    PROFILE_HERE;
    f0();
}

template <>
void f<0>()
{
    PROFILE_HERE;
    f0();
}

int main()
{
    set_prof();
    f<13>();
    cout << "n= " << n <<endl;
    ofstream file("profile.out");
    print_profile( file );
    //    print_profile( cout );
}
