#ifndef _TEST_H_
#define _TEST_H_

template <class T>
class Test
{
public:
	T a, b;
	explicit Test(const T a, const T b) : a(a), b(b) {}
	T plus() const;
};

template <class T>
T Test<T>::plus() const
{
	return a + b;
}
#endif