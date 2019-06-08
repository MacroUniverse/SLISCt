#pragma once
#include <iostream>
#include <cmath>
#include <functional>
#include <utility>
#include <complex>
using std::cout; using std::endl; //using std::complex;

namespace shit
{
	struct data
	{
		double x;
	};

	void fun(data s)
	{
		cout << " fun(s) = " << s.x << endl;
	}

	void fun2(double x)
	{
		cout << "x = " << x << endl;
	}
}

void fun(shit::data s)
{
	cout << " fun(s) = " << s.x + 1. << endl;
}

void fun2(double x)
{
	cout << "x = " << x + 1. << endl;
}

template <class T>
void move(T x)
{
	cout << "inside move(T)" << endl;
}

void myfun()
{
	cout << "hello, world" << endl;
	//complex<double> s(0.,3.1415926535);
	//cout << "s = " << s << endl;
	//cout << "exp(s) = " << exp(s) << endl;
	//cout << "sin(x) = " << sin(1.) << endl;

	//shit::data s; s.x = 5.5;
	//cout << "s.x = " << s.x << endl;
	//cout << "call fun()..." << endl;
	//shit::fun(s);

	//cout << "call fun2()..." << endl;
	//shit::fun2(3.14);
	move(std::complex<double>(1.3,2.6));
}
