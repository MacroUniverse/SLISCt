#pragma once
#include "../global.h"
#include "../scalar_arith.h"

namespace slisc {

// ======== translation details ============
// variables that subtracted 1 comparing to F77 version
// arguments: iexph
// internal variables: i, j, k, icoef, ih2, ip, iq, iput, iget, ifree, iused

void ZGPADM(Int_I ideg, Int_I m, Doub_I t, const Comp *H, Int_I ldh, Comp *wsp, Int_I lwsp,
    Int *ipiv, Int_O iexph, Int_O ns, Int_O iflag)
{
    Int i, j, k, icoef, mm, ih2, iodd, iused, ifree, iq, ip, iput, iget;
    Doub hnorm;
    Comp cp, cq, scale, scale2, temp;
    const Comp zero = 0., one = 1.;

    mm = m*m;
    iflag = 0;
    if (ldh < m) iflag = -1;
    if (lwsp < 4 * mm + ideg + 1) iflag = -2;
    if (iflag != 0)
        SLS_ERR("bad sizes (in input of ZGPADM)");

    icoef = 0;
    ih2 = icoef + ideg + 1;
    ip = ih2 + mm;
    iq = ip + mm;
    ifree = iq + mm;

    for (i = 0; i < m; ++i) {
        wsp[i] = zero;
    }
    for (j = 0; j < m; ++j) {
        for (i = 0; i < m; ++i) {
            wsp[i] = wsp[i] + abs(H[i + ldh*j]);
        }
    }

    hnorm = 0.;
    for (i = 0; i < m; ++i) {
        hnorm = MAX(hnorm, real(wsp[i]));
    }

    hnorm = abs(t*hnorm);
    if (hnorm == 0.)
        SLS_ERR("Error - null H in input of ZGPADM.");
    ns = MAX(0, Int(log(hnorm) / log(2.)) + 2);
    scale = Comp(t / pow(2, ns), 0.);
    scale2 = scale*scale;

    i = ideg;
    j = 2 * ideg;
    wsp[icoef] = one;
    for (k = 0; k < ideg; ++k) {
        wsp[icoef + k + 1] = (wsp[icoef + k] * Doub(i - k)) / Doub((k + 1)*(j - k));
    }

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, &scale2, H, ldh,
        H, ldh, &zero, wsp + ih2, m);

    cp = wsp[icoef + ideg - 1];
    cq = wsp[icoef + ideg];
    for (j = 0; j < m; ++j) {
        for (i = 0; i < m; ++i) {
            wsp[ip + j*m + i] = zero;
            wsp[iq + j*m + i] = zero;
        }
        wsp[ip + j*(m + 1)] = cp;
        wsp[iq + j*(m + 1)] = cq;
    }

    iodd = 1;
    k = ideg - 2;

    /*100*/
    do {
        iused = iodd*(iq - ip) + ip;
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, &one, wsp + iused, m,
            wsp + ih2, m, &zero, wsp + ifree, m);
        for (j = 0; j < m; ++j)
            wsp[ifree + j*(m + 1)] = wsp[ifree + j*(m + 1)] + wsp[icoef + k];

        ip = ifree + iodd*(ip - ifree);
        iq = iodd*(ifree - iq) + iq;
        ifree = iused;
        iodd = 1 - iodd;
        --k;
    } while (k + 1 > 0);


    if (iodd != 0) {
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, &scale, wsp + iq, m,
            H, ldh, &zero, wsp + ifree, m);
        iq = ifree;
    }
    else {
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, &scale, wsp + ip, m,
            H, ldh, &zero, wsp + ifree, m);
        ip = ifree;
    }
    temp = -one;
    cblas_zaxpy(mm, &temp, wsp + ip, 1, wsp + iq, 1);
    iflag = LAPACKE_zgesv(LAPACK_COL_MAJOR, m, m, (double _Complex*)(wsp + iq), m, ipiv,
        (double _Complex*)(wsp + ip), m);
    if (iflag != 0)
        SLS_ERR("Problem in ZGESV (within ZGPADM)");
    cblas_zdscal(mm, 2., wsp + ip, 1);
    for (j = 0; j < m; ++j)
        wsp[ip + j*(m + 1)] = wsp[ip + j*(m + 1)] + one;

    iput = ip;
    if (ns == 0 && iodd != 0) {
        cblas_zdscal(mm, -1., wsp + ip, 1);
    }
    else {
        iodd = 1;
        for (k = 0; k < ns; ++k) {
            iget = iodd*(ip - iq) + iq;
            iput = ip + iodd*(iq - ip);
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, &one, wsp + iget, m,
                wsp + iget, m, &zero, wsp + iput, m);
            iodd = 1 - iodd;
        }
    }

    iexph = iput;
}

} // namespace slisc
