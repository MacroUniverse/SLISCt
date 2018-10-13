// for basic matrix arithmatics
#pragma once
#include "../slisc.h"

void mul(MatDoub_O &v, MatDoub_I &v1, MatDoub_I &v2);

void mul(MatComp_O &v, MatComp_I &v1, MatComp_I &v2);

void mul(MatComp_O &v, MatComp_I &v1, MatDoub_I &v2);

void mul(MatComp_O &v, MatDoub_I &v1, MatComp_I &v2);
