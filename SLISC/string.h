// string utilities

#pragma once
#include "global.h"

namespace slisc {

// string utilities

template <typename T>
inline void num2str(Str_O str, T s)
{
	str = to_string(s);
	if (str.find('.') != Str::npos)
		str.erase(str.find_last_not_of('0') + 1);
}

} // namespace slisc
