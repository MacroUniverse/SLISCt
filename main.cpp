// comprehensive test of nr3plus.h

#include "nr3.h"
#include "nr3plus.h"

using namespace std;

int main()
{
  VecDoub v(3,0.); VecComplex vc(3,0.);
  linspace(v, 0., PI);
  disp(v, 8);
  disp(v, 0, 2, 8);
  linspace(vc, 0., PI + 2.*PI*I);
  disp(vc, 8);
  disp(vc, 0, 2, 8);
  Complex c, c1(3.,5.);
  c = c1 + 1;
  cout << c << endl;
  c = c1 - 1;
  cout << c << endl;
  c = 1 - c1;
  cout << c << endl;
  c = c1 * 2;
  cout << c << endl;
  c = 2 * c1;
  cout << c << endl;
  c = 2/c1;
  cout << c << endl;
  MatDoub a(3,3,0.);
  linspace(a, 0., 8.);
  disp(a, 5);
}
