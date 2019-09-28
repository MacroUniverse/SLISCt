// === Vector Class ===

// copy constructor
template <class T>
Vector<T>::Vector(const Vector<T> &rhs) : Vector(rhs.nn)
{
    memcpy(v, rhs.v, nn*sizeof(T));
}

// move constructor
template <class T>
Vector<T>::Vector(Vector<T> &&rhs) : nn(rhs.nn), v(rhs.v)
{
    rhs.v = nullptr;
}

// move assignment
template <class T>
Vector<T> & Vector<T>::operator=(Vector &&rhs)
{
    if (v != nullptr) delete[] v;
    nn = rhs.nn;
    v = rhs.v;
    rhs.v = nullptr;
    return *this;
}

// === Matrix Class ===

// copy constructor
template <class T>
Matrix<T>::Matrix(const Matrix<T> &rhs) : Matrix(rhs.nn, rhs.mm)
{
    memcpy(v[0], rhs.v[0], nn*mm*sizeof(T));
}

// move constructor
template <class T>
Matrix<T>::Matrix(Matrix<T> &&rhs) : nn(rhs.nn), mm(rhs.mm), v(rhs.v)
{
    rhs.v = nullptr;
}

// move assignment
template <class T>
Matrix<T> & Matrix<T>::operator=(Matrix<T> &&rhs)
{
    data_free(v);
    nn = rhs.nn; mm = rhs.mm; v = rhs.v;
    rhs.v = nullptr;
    return *this;
}

// === 3D Matrix Class ===

// copy constructor
template <class T>
NRMat3d<T>::NRMat3d(const NRMat3d<T> &rhs) : NRMat3d(rhs.nn, rhs.mm, rhs.kk)
{
    memcpy(v[0][0], rhs.v[0][0], nn*mm*kk*sizeof(T));
}

// move constructor
template <class T>
NRMat3d<T>::NRMat3d(NRMat3d<T> &&rhs) : nn(rhs.nn), mm(rhs.mm), kk(rhs.kk), v(rhs.v)
{
    rhs.v = nullptr;
}

// move assignment operator
template <class T>
NRMat3d<T> &NRMat3d<T>::operator=(NRMat3d<T> &&rhs)
{
    data_free(v);
    nn = rhs.nn; mm = rhs.mm; kk = rhs.kk; v = rhs.v;
    rhs.v = nullptr;
    return *this;
}

// === functions that returns matrixs or vectors ===

template <class T>
inline Vector<T> operator+(const Vector<T> &v1, Doub_I s)
{
    Sizet i, N{ v1.size() };
    Vector<T> v(N);
    for (i = 0; i < N; ++i)
        v[i] = v1[i] + s;
    return v;
}

template <class T>
inline Vector<T> operator+(Doub_I s, const Vector<T> &v1)
{
    Sizet i, N{ v1.size() };
    Vector<T> v(N);
    for (i = 0; i < N; ++i)
        v[i] = v1[i] + s;
    return v;
}

template <class T>
inline Vector<T> operator-(const Vector<T> &v1, Doub_I s)
{
    Sizet i, N{ v1.size() };
    Vector<T> v(N);
    for (i = 0; i < N; ++i)
        v[i] = v1[i] - s;
    return v;
}

template <class T>
inline Vector<T> operator-(Doub_I s, const Vector<T> &v1)
{
    Sizet i, N{ v1.size() };
    Vector<T> v(N);
    for (i = 0; i < N; ++i)
        v[i] = s - v1[i];
    return v;
}

template <class T>
inline Vector<T> operator*(const Vector<T> &v1, Doub_I s)
{
    Sizet i, N{ v1.size() };
    Vector<T> v(N);
    for (i = 0; i < N; ++i)
        v[i] = v1[i] - s;
    return v;
}

template <class T>
inline Vector<T> operator*(Doub_I s, const Vector<T> &v1)
{
    Sizet i, N{ v1.size() };
    Vector<T> v(N);
    for (i = 0; i < N; ++i)
        v[i] = v1[i] - s;
    return v;
}

template <class T>
inline Vector<T> operator/(const Vector<T> &v1, Doub_I s)
{
    Int i, N{ v1.size() };
    Doub sInv = 1./s;
    Vector<T> v(N);
    for (i = 0; i < N; ++i)
        v[i] = v1[i] * sInv;
    return v;
}

template <class T>
inline Vector<T> operator/(Doub_I s, const Vector<T> &v1)
{
    Sizet i, N{ v1.size() };
    Vector<T> v(N);
    for (i = 0; i < N; ++i)
        v[i] = s / v1[i];
    return v;
}