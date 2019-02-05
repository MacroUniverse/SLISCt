// type and class definitions and dependencies
// this header file can be used alone

#pragma once

#ifndef NDEBUG
// this will not check the last index
#define _CHECKBOUNDS_
#endif

// all the system #include's we'll ever need
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <limits>
#include "typedef.h"

namespace slisc
{

// NaN definition
static const Doub NaN = std::numeric_limits<Doub>::quiet_NaN();

// === constants ===

const Doub PI = 3.14159265358979323;
const Doub E = 2.71828182845904524;
const Comp I(0., 1.);

// report error and pause execution
#define error(str) do{std::cout << "error: " << __FILE__ << ": line " << __LINE__ << ": " << str << std::endl; getchar();} while(0)

#define warning(str) do{std::cout << "warning: " << __FILE__ << ": line " << __LINE__ << ": " << str << std::endl;} while(0)

}

#include "vector.h"
#include "matrix.h"
#include "mat3d.h"
