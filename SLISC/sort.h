#pragma once
#include "scalar_arith.h"
#include "vector.h"

namespace slisc {
// ======== declarations =========

// sort v in ascending order, and rearrange v1 at the same time
template<class T, class U>
void sort2(Vector<T> &v, Vector<U> &v1);

template<class T, class U>
void sort2(vector_IO<T> v, vector_IO<U> v1);

// sort v, and output the original index of each element


// ======== implementations =======
template<class T>
void sort(Vector<T> &arr, Long m = -1)
{
    static const Int M = 7, NSTACK = 64;
    Long i, ir, j, k, jstack = -1, l = 0, n = arr.size();
    T a;
    VecLong istack(NSTACK);
    if (m>0) n = MIN(m, n);
    ir = n - 1;
    for (;;) {
        if (ir - l < M) {
            for (j = l + 1; j <= ir; j++) {
                a = arr[j];
                for (i = j - 1; i >= l; i--) {
                    if (arr[i] <= a) break;
                    arr[i + 1] = arr[i];
                }
                arr[i + 1] = a;
            }
            if (jstack < 0) break;
            ir = istack[jstack--];
            l = istack[jstack--];
        }
        else {
            k = (l + ir) >> 1;
            swap(arr[k], arr[l + 1]);
            if (arr[l] > arr[ir]) {
                swap(arr[l], arr[ir]);
            }
            if (arr[l + 1] > arr[ir]) {
                swap(arr[l + 1], arr[ir]);
            }
            if (arr[l] > arr[l + 1]) {
                swap(arr[l], arr[l + 1]);
            }
            i = l + 1;
            j = ir;
            a = arr[l + 1];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                swap(arr[i], arr[j]);
            }
            arr[l + 1] = arr[j];
            arr[j] = a;
            jstack += 2;
            if (jstack >= NSTACK) throw("NSTACK too small in sort.");
            if (ir - i + 1 >= j - l) {
                istack[jstack] = ir;
                istack[jstack - 1] = i;
                ir = j - 1;
            }
            else {
                istack[jstack] = j - 1;
                istack[jstack - 1] = l;
                l = i;
            }
        }
    }
}

template<class T, class U>
void sort2_vv(T *v, U *v1, Long_I N)
{
    const Long M = 7, NSTACK = 64;
    Long i, ir, j, k, jstack = -1, l = 0, n = N;
    T a;
    U b;
    VecLong istack(NSTACK);
    ir = n - 1;
    for (;;) {
        if (ir - l < M) {
            for (j = l + 1; j <= ir; j++) {
                a = v[j];
                b = v1[j];
                for (i = j - 1; i >= l; i--) {
                    if (v[i] <= a) break;
                    v[i + 1] = v[i];
                    v1[i + 1] = v1[i];
                }
                v[i + 1] = a;
                v1[i + 1] = b;
            }
            if (jstack < 0) break;
            ir = istack[jstack--];
            l = istack[jstack--];
        }
        else {
            k = (l + ir) >> 1;
            swap(v[k], v[l + 1]);
            swap(v1[k], v1[l + 1]);
            if (v[l] > v[ir]) {
                swap(v[l], v[ir]);
                swap(v1[l], v1[ir]);
            }
            if (v[l + 1] > v[ir]) {
                swap(v[l + 1], v[ir]);
                swap(v1[l + 1], v1[ir]);
            }
            if (v[l] > v[l + 1]) {
                swap(v[l], v[l + 1]);
                swap(v1[l], v1[l + 1]);
            }
            i = l + 1;
            j = ir;
            a = v[l + 1];
            b = v1[l + 1];
            for (;;) {
                do i++; while (v[i] < a);
                do j--; while (v[j] > a);
                if (j < i) break;
                swap(v[i], v[j]);
                swap(v1[i], v1[j]);
            }
            v[l + 1] = v[j];
            v[j] = a;
            v1[l + 1] = v1[j];
            v1[j] = b;
            jstack += 2;
            if (jstack >= NSTACK) throw("NSTACK too small in sort2.");
            if (ir - i + 1 >= j - l) {
                istack[jstack] = ir;
                istack[jstack - 1] = i;
                ir = j - 1;
            }
            else {
                istack[jstack] = j - 1;
                istack[jstack - 1] = l;
                l = i;
            }
        }
    }
}

template<class T, class U>
void sort2(Vector<T> &v, Vector<U> &v1)
{ sort2_vv(v.ptr(), v1.ptr(), v1.size()); }

template<class T, class U>
void sort2(vector_IO<T> v, vector_IO<U> v1)
{ sort2_vv(v.data(), v1.data(), v1.size()); }

template<class T>
void shell(Vector<T> &a, Int m = -1)
{
    Int i, j, inc, n = a.size();
    T v;
    if (m>0) n = MIN(m, n);
    inc = 1;
    do {
        inc *= 3;
        inc++;
    } while (inc <= n);
    do {
        inc /= 3;
        for (i = inc; i<n; i++) {
            v = a[i];
            j = i;
            while (a[j - inc] > v) {
                a[j] = a[j - inc];
                j -= inc;
                if (j < inc) break;
            }
            a[j] = v;
        }
    } while (inc > 1);
}

namespace hpsort_util {
    template<class T>
    void sift_down(Vector<T> &ra, const Int l, const Int r)
    {
        Int j, jold;
        T a;
        a = ra[l];
        jold = l;
        j = 2 * l + 1;
        while (j <= r) {
            if (j < r && ra[j] < ra[j + 1]) j++;
            if (a >= ra[j]) break;
            ra[jold] = ra[j];
            jold = j;
            j = 2 * j + 1;
        }
        ra[jold] = a;
    }
}

template<class T>
void hpsort(Vector<T> &ra)
{
    Int i, n = ra.size();
    for (i = n / 2 - 1; i >= 0; i--)
        hpsort_util::sift_down(ra, i, n - 1);
    for (i = n - 1; i>0; i--) {
        swap(ra[0], ra[i]);
        hpsort_util::sift_down(ra, 0, i - 1);
    }
}

template<class T>
void piksrt(Vector<T> &arr)
{
    Int i, j, n = arr.size();
    T a;
    for (j = 1; j<n; j++) {
        a = arr[j];
        i = j;
        while (i > 0 && arr[i - 1] > a) {
            arr[i] = arr[i - 1];
            i--;
        }
        arr[i] = a;
    }
}

template<class T, class U>
void piksr2(Vector<T> &arr, Vector<U> &brr)
{
    Int i, j, n = arr.size();
    T a;
    U b;
    for (j = 1; j<n; j++) {
        a = arr[j];
        b = brr[j];
        i = j;
        while (i > 0 && arr[i - 1] > a) {
            arr[i] = arr[i - 1];
            brr[i] = brr[i - 1];
            i--;
        }
        arr[i] = a;
        brr[i] = b;
    }
}

struct Indexx {
    Int n;
    VecInt indx;

    template<class T> Indexx(const Vector<T> &arr) {
        index(&arr[0], arr.size());
    }
    Indexx() : indx(0) {}

    template<class T> void sort(Vector<T> &brr) {
        if (brr.size() != n) throw("bad size in Index sort");
        Vector<T> tmp(brr);
        for (Int j = 0; j<n; j++) brr[j] = tmp[indx[j]];
    }

    template<class T> inline const T & el(Vector<T> &brr, Int j) const {
        return brr[indx[j]];
    }
    template<class T> inline T & el(Vector<T> &brr, Int j) {
        return brr[indx[j]];
    }

    template<class T> void index(const T *arr, Int nn);

    void rank(VecInt_O irank) {
        irank.resize(n);
        for (Int j = 0; j<n; j++) irank[indx[j]] = j;
    }

};

template<class T>
void Indexx::index(const T *arr, Int nn)
{
    const Int M = 7, NSTACK = 64;
    Int i, indxt, ir, j, k, jstack = -1, l = 0;
    T a;
    VecInt istack(NSTACK);
    n = nn;
    indx.resize(n);
    ir = n - 1;
    for (j = 0; j<n; j++) indx[j] = j;
    for (;;) {
        if (ir - l < M) {
            for (j = l + 1; j <= ir; j++) {
                indxt = indx[j];
                a = arr[indxt];
                for (i = j - 1; i >= l; i--) {
                    if (arr[indx[i]] <= a) break;
                    indx[i + 1] = indx[i];
                }
                indx[i + 1] = indxt;
            }
            if (jstack < 0) break;
            ir = istack[jstack--];
            l = istack[jstack--];
        }
        else {
            k = (l + ir) >> 1;
            swap(indx[k], indx[l + 1]);
            if (arr[indx[l]] > arr[indx[ir]]) {
                swap(indx[l], indx[ir]);
            }
            if (arr[indx[l + 1]] > arr[indx[ir]]) {
                swap(indx[l + 1], indx[ir]);
            }
            if (arr[indx[l]] > arr[indx[l + 1]]) {
                swap(indx[l], indx[l + 1]);
            }
            i = l + 1;
            j = ir;
            indxt = indx[l + 1];
            a = arr[indxt];
            for (;;) {
                do i++; while (arr[indx[i]] < a);
                do j--; while (arr[indx[j]] > a);
                if (j < i) break;
                swap(indx[i], indx[j]);
            }
            indx[l + 1] = indx[j];
            indx[j] = indxt;
            jstack += 2;
            if (jstack >= NSTACK) throw("NSTACK too small in index.");
            if (ir - i + 1 >= j - l) {
                istack[jstack] = ir;
                istack[jstack - 1] = i;
                ir = j - 1;
            }
            else {
                istack[jstack] = j - 1;
                istack[jstack - 1] = l;
                l = i;
            }
        }
    }
}

template<class T>
T select(const Int k, Vector<T> &arr)
{
    Int i, ir, j, l, mid, n = arr.size();
    T a;
    l = 0;
    ir = n - 1;
    for (;;) {
        if (ir <= l + 1) {
            if (ir == l + 1 && arr[ir] < arr[l])
                swap(arr[l], arr[ir]);
            return arr[k];
        }
        else {
            mid = (l + ir) >> 1;
            swap(arr[mid], arr[l + 1]);
            if (arr[l] > arr[ir])
                swap(arr[l], arr[ir]);
            if (arr[l + 1] > arr[ir])
                swap(arr[l + 1], arr[ir]);
            if (arr[l] > arr[l + 1])
                swap(arr[l], arr[l + 1]);
            i = l + 1;
            j = ir;
            a = arr[l + 1];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                swap(arr[i], arr[j]);
            }
            arr[l + 1] = arr[j];
            arr[j] = a;
            if (j >= k) ir = j - 1;
            if (j <= k) l = i;
        }
    }
}

struct Heapselect {
    Int m, n, srtd;
    VecDoub heap;

    Heapselect(Int mm) : m(mm), n(0), srtd(0), heap(mm, 1.e99) {}

    void add(Doub val) {
        Int j, k;
        if (n<m) {
            heap[n++] = val;
            if (n == m) sort(heap);
        }
        else {
            if (val > heap[0]) {
                heap[0] = val;
                for (j = 0;;) {
                    k = (j << 1) + 1;
                    if (k > m - 1) break;
                    if (k != (m - 1) && heap[k] > heap[k + 1]) k++;
                    if (heap[j] <= heap[k]) break;
                    swap(heap[k], heap[j]);
                    j = k;
                }
            }
            n++;
        }
        srtd = 0;
    }

    Doub report(Int k) {
        Int mm = MIN(n, m);
        if (k > mm - 1) throw("Heapselect k too big");
        if (k == m - 1) return heap[0];
        if (!srtd) { sort(heap); srtd = 1; }
        return heap[mm - 1 - k];
    }
};

} // namespace slisc
