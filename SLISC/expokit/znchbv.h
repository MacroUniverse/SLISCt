// the translated ZNCHBV is not complete yet!
#pragma once
#include "../global.h"

namespace slisc {

// translate as directly as possible for now:
// * any variable values should not change
// * all pointer indexing (including +, -) is subtracted by 1 in place

void ZNCHBV(Int_I m, Doub_I t, const Comp *H, Int_I ldh, Comp *y, Comp *wsp)
{
    SLS_ERR("the translated ZNCHBV is not complete yet, do it now!");
    const Comp zero = 0.;
    const Int ndeg = 7;
    Int i, j, k, ip, ih, iy, iz;
    Doub alpha0;
    Comp alpha[ndeg], theta[ndeg];
    Comp tmpc;
      
    ih = 1;
    iy = ih + m*m;
    iz = iy + m;


    alpha0 = 0.183216998528140087e-11;
    alpha[0] = Comp(0.557503973136501826e+02, -0.204295038779771857e+03);
    alpha[1] = Comp(-0.938666838877006739e+02, 0.912874896775456363e+02);
    alpha[2] = Comp(0.469965415550370835e+02, -0.116167609985818103e+02);
    alpha[3] = Comp(-0.961424200626061065e+01, -0.264195613880262669e+01);
    alpha[4] = Comp(0.752722063978321642e+00, 0.670367365566377770e+00);
    alpha[5] = Comp(-0.188781253158648576e-01, -0.343696176445802414e-01);
    alpha[6] = Comp(0.143086431411801849e-03, 0.287221133228814096e-03);

    theta[0] = Comp(-0.562314417475317895e+01, 0.119406921611247440e+01);
    theta[1] = Comp(-0.508934679728216110e+01, 0.358882439228376881e+01);
    theta[2] = Comp(-0.399337136365302569e+01, 0.600483209099604664e+01);
    theta[3] = Comp(-0.226978543095856366e+01, 0.846173881758693369e+01);
    theta[4] = Comp(0.208756929753827868e+00, 0.109912615662209418e+02);
    theta[5] = Comp(0.370327340957595652e+01, 0.136563731924991884e+02);
    theta[6] = Comp(0.889777151877331107e+01, 0.166309842834712071e+02);

    for (ip = 1; ip <= ndeg; ++ip) {
        theta[ndeg + ip - 1] = conj(theta[ip - 1]);
        alpha[ndeg + ip - 1] = conj(alpha[ip - 1]);
    }

    for (j = 1; j <= m; ++j) {
        wsp[iz + j - 2] = y[j - 1];
        y[j - 1] = y[j - 1] * alpha0;
    }

    for (ip = 1; ip <= 2 * ndeg; ++ip) {
        alpha[ip - 1] = 0.5*alpha[ip - 1];
        for (j = 1; j <= m; ++j) {
            wsp[iy + j - 2] = wsp[iz + j - 2];
            for (i = 1; i <= MIN(j + 1, m); ++i) {
                wsp[ih + (j - 1)*m + i - 2] = -t*H[i + ldh*(j-1) - 1];
            }
            wsp[ih + (j - 1)*m + j - 2] = wsp[ih + (j - 1)*m + j - 2] - theta[ip - 1];
            for (k = i; k <= m; ++k) {
                wsp[ih + (j - 1)*m + k - 2] = zero;
            }
        }
        for (j = 1; j <= m; ++j) {

        }
        for (i = 1; i <= m - 1; ++i) {
            if (abs(wsp[ih + (i - 1)*m + i - 2]) < abs(wsp[ih + (i - 1)*m + i - 1])) {
                cblas_zswap(m - i + 1, wsp + ih + (i - 1)*m + i - 2, m,
                    wsp + ih + (i - 1)*m + i - 1, m);
                cblas_zswap(1, wsp + iy + i - 2, 1, wsp + iy + i - 1, 1);
            }

            tmpc = wsp[ih + (i - 1)*m + i - 1] / wsp[ih + (i - 1)*m + i - 2];
            Comp temp = -tmpc;
            cblas_zaxpy(m - i, &temp, wsp + ih + i*m + i - 2, m, wsp + ih + i*m + i - 1, m);
            wsp[iy + i - 1] = wsp[iy + i - 1] - tmpc*wsp[iy + i - 2];
        }

        for (i = m; i >= 1; --i) {
            tmpc = wsp[iy + i - 2];
            for (j = i + 1; j <= m; ++j) {
                tmpc = tmpc - wsp[ih + (j - 1)*m + i - 2]*wsp[iy + j - 2];
            }
            wsp[iy + i - 2] = tmpc / wsp[ih + (i - 1)*m + i - 2];
        }
        for (j = 1; j <= m; ++j) {
            y[j - 1] = y[j - 1] + alpha[ip - 1]*wsp[iy + j - 2];
        }
    }
}

} // namespace slisc
