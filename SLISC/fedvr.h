#pragma once
#include "cmat.h"
#include "matcooh.h"

namespace slisc {

// generate Gauss-Lobatto abscissas x and weights w
// x is in [-1,1]
// data from https://keisan.casio.com/exec/system/1280801905
// use 38 digits data in case long double is needed
inline void GaussLobatto(VecDoub_O x, VecDoub_O w)
{
#ifdef SLS_CHECK_SHAPE
    if (x.size() != w.size())
        SLS_ERR("wrong shape!");
#endif
    Long_I N = x.size();
    Long N2 = N/2;
    x[0] = -1; x[N - 1] = 1;
    // only set 0 <= x < 1 and w for x >= 0! Others will be calculated. index is 0 based!
    if (N == 4) {
        x[2] = 0.44721359549995793928183473374625524709;

        w[2] = 0.83333333333333333333333333333333333333;
        w[3] = 0.16666666666666666666666666666666666667;
    }
    else if (N == 6) {
        x[3] = 0.28523151648064509631415099404087907192;
        x[4] = 0.76505532392946469285100297395933815037;

        w[3] = 0.55485837703548635301672052512130868086;
        w[4] = 0.37847495629784698031661280821202465248;
        w[5] = 0.066666666666666666666666666666666666667;
    }
    else if (N == 8) {
        x[4] = 0.2092992179024788687686572603453512553;
        x[5] = 0.59170018143314230214451073139795318995;
        x[6] = 0.8717401485096066153374457612206634381;

        w[4] = 0.41245879465870388156705297140220978948;
        w[5] = 0.34112269248350436476424067710774817178;
        w[6] = 0.2107042271435060393829920657757563245;
        w[7] = 0.035714285714285714285714285714285714286;
    }
    else if (N == 10) {
        x[5] = 0.16527895766638702462621976595817353323;
        x[6] = 0.47792494981044449566117509273125799789;
        x[7] = 0.73877386510550507500310617485983072502;
        x[8] = 0.91953390816645881382893266082233813415;

        w[5] = 0.32753976118389745665651052791689314471;
        w[6] = 0.29204268367968375787558225737444389221;
        w[7] = 0.22488934206312645211945782173104784275;
        w[8] = 0.13330599085107011112622717075539289811;
        w[9] = 0.022222222222222222222222222222222222222;
    }
    else if (N == 12) {
        x[6] = 0.1365529328549275548640618557396938969;
        x[7] = 0.3995309409653489322643497915669669005;
        x[8] = 0.63287615303186067766240485444365585824;
        x[9] = 0.81927932164400667834864158171690266069;
        x[10] = 0.94489927222288222340758013830321871361;

        w[6] = 0.2714052409106961770002883384996028404;
        w[7] = 0.25127560319920128029324441214759618333;
        w[8] = 0.2125084177610211453583020773668662702;
        w[9] = 0.15797470556437011516467106270034026825;
        w[10] = 0.09168451741319613066834259413407928631;
        w[11] = 0.015151515151515151515151515151515151515;
    }
    else if (N == 14) {
        x[7] = 0.11633186888370386765877670973616016794;
        x[8] = 0.34272401334271284504390340364167464483;
        x[9] = 0.55063940292864705531662270585908063446;
        x[10] = 0.72886859909132614058467240052088159566;
        x[11] = 0.86780105383034725100022020290826421325;
        x[12] = 0.95993504526726090135510016201542438907;

        w[7] = 0.2316127944684570588896283572926391903;
        w[8] = 0.21912625300977075487116252395416761296;
        w[9] = 0.1948261493734161186403317783758845121;
        w[10] = 0.1600218517629521424128209979875946405;
        w[11] = 0.11658665589871165154099667065465023;
        w[12] = 0.0668372844976812846340706607460528251;
        w[13] = 0.010989010989010989010989010989010989011;
    }
    else if (N == 16) {
        x[8] = 0.10132627352194944784303300504591776253;
        x[9] = 0.29983046890076320809835345472230064782;
        x[10] = 0.48605942188713761178189078584687469689;
        x[11] = 0.65238870288249308946788321964058148032;
        x[12] = 0.79200829186181506393108827096314570581;
        x[13] = 0.89920053309347209299462826151984947675;
        x[14] = 0.96956804627021793295224273836745924139;

        w[8] = 0.20195830817822987148919912541094012055;
        w[9] = 0.19369002382520358431691359885352203031;
        w[10] = 0.17749191339170412530107566952835777016;
        w[11] = 0.15402698080716428081564494048499451545;
        w[12] = 0.12425538213251409834953633265731320077;
        w[13] = 0.0893936973259308009910520801660837149;
        w[14] = 0.0508503610059199054032449195654553145;
        w[15] = 0.0083333333333333333333333333333333333333;
    }
    else {
        SLS_ERR("no data!");
    }
    for (Long i = 1; i < N2; ++i)
        x[i] = -x[N - i - 1];
    for (Long i = 0; i < N2; ++i)
        w[i] = w[N - i - 1];
}


// get derivatives of legendre interpolation polynomials t abscissas
// df is an NxN matrix
// f'_i(x_j) = df(j,i)
inline void legendre_interp_der(CmatDoub_O df, VecDoub_I x)
{
    Long i, j, k, N{ x.size() };
    Doub t;
#ifdef SLS_CHECK_SHAPE
    if (df.n1() != N || df.n2() != N)
        SLS_ERR("wrong shape!");
#endif
    for (i = 0; i < N; ++i)
        for (j = 0; j < N; ++j) {
            if (j != i) {
                t = 1.;
                for (k = 0; k < i; ++k) {
                    if (k == j) t /= x[i] - x[k];
                    else t *= (x[j] - x[k]) / (x[i] - x[k]);
                }
                for (k = i + 1; k < N; ++k) {
                    if (k == j) t /= x[i] - x[k];
                    else t *= (x[j] - x[k]) / (x[i] - x[k]);
                }
            }
            else {
                t = 0.;
                for (k = 0; k < i; ++k)
                    t += 1. / (x[i] - x[k]);
                for (k = i + 1; k < N; ++k)
                    t += 1. / (x[i] - x[k]);
            }
            df(j,i) = t;
        }
}

// calculate FEDVR global index from FE index i and DVR index j
inline Long indFEDVR(Long_I i, Long_I j, Long_I Ngs)
{ return (Ngs-1) * i + j - 1; }

// generate FEDVR grid and weight
// 'Nfe' is the number of finite elements
// 'Ngs' is the number of abscissas in each element (including 2 boundaries)
// 'bounds' size = Nfe+1, are the boundary points of finite elements
// 'x0','w0' are grid points and weights in [-1,1]
// 'x','w' are the global grid points and weights

inline void FEDVR_grid(VecDoub_O x, VecDoub_O w, VecDoub_I wFE, VecDoub_I xFE, VecDoub_I x0, VecDoub_I w0)
{
    Long i, j, k = 0;
    Long Ngs = x0.size();
    Long Nfe = wFE.size();
    Long Nx = Nfe * (Ngs - 1) - 1;
#ifdef SLS_CHECK_SHAPE
    if (x.size() != Nx || w.size() != Nx)
        SLS_ERR("wrong shape!");
#endif
    Doub a, b;
    for (i = 0; i < Nfe; ++i) {
        a = wFE(i);
        b = xFE(i);
        if (i > 0)
            w[k-1] += a*w0(0);
        for (j = 1; j < Ngs; ++j) {
            x[k] = a*x0[j] + b;
            w[k] = a*w0[j];
            if (++k == Nx) break;
        }
    }
}

// number of non-zero elements in fedvr second derivative matrix
inline Long fedvr_d2_nnz(Long_I Ngs, Long_I Nfe)
{
    return (Ngs*Ngs - 1)*Nfe - 4 * Ngs + 3;
}

// sparse second derivative matrix for normalized FEDVR basis
inline void D2_matrix(McooDoub_O D2, VecDoub_I w0, VecDoub_I wFE, CmatDoub_I df)
{
    Long i, j, k, m, n, mm, nn;
    Long Nfe = wFE.size();
    Long Ngs = w0.size();
#ifdef SLS_CHECK_SHAPE
    Long Nx = Nfe * (Ngs - 1) - 1;
    if (D2.n1() != Nx || D2.n1() != Nx)
        SLS_ERR("wrong shape!");
#endif

    Doub s, coeff;

    // prepare basic block
    CmatDoub block(Ngs, Ngs);
    for (j = 0; j < Ngs; ++j) {
        for (i = 0; i <= j; ++i) {
            s = 0.;
            for (k = 0; k < Ngs; ++k)
                s += w0(k) * df(k, i) * df(k, j);
            block(i, j) = s;
        }
    }

    // calculate Kinetic matrix T
    D2.reserve(fedvr_d2_nnz(Ngs, Nfe)); // # non-zero elements
    for (i = 0; i < Nfe; ++i) {
        // blocks without boundary
        for (n = 1; n < Ngs - 1; ++n) {
            coeff = -1. / SQR(wFE(i));
            for (m = 1; m <= n; ++m) {
                mm = indFEDVR(i, m, Ngs); nn = indFEDVR(i, n, Ngs);
                s = coeff * block(m,n);
                D2.push(s, mm, nn);
                if (mm != nn) {
                    D2.push(s, nn, mm);
                }
            }
        }
    }

    for (i = 0; i < Nfe - 1; ++i) {
        // block right boundary
        n = Ngs - 1;
        coeff = -1. / (pow(wFE(i), 1.5) * ::sqrt(wFE(i) + wFE(i + 1)));
        for (m = 1; m < n; ++m) {
            mm = indFEDVR(i, m, Ngs); nn = indFEDVR(i, n, Ngs);
            s = coeff * block(m, n);
            D2.push(s, mm, nn); D2.push(s, nn, mm);
        }

        // block lower right corner
        m = Ngs - 1;
        coeff = -1. / (wFE(i) + wFE(i + 1));
        mm = indFEDVR(i, m, Ngs);
        s = coeff*(block(m,n)/wFE(i) + block(0,0)/wFE(i+1));
        D2.push(s, mm, mm);

        // block upper boundary
        coeff = -1. / (pow(wFE(i + 1), 1.5) * ::sqrt(wFE(i) + wFE(i + 1)));
        for (n = 1; n < Ngs - 1; ++n) {
            mm = indFEDVR(i, m, Ngs); nn = indFEDVR(i + 1, n, Ngs);
            s = coeff * block(0, n);
            D2.push(s, mm, nn); D2.push(s, nn, mm);
        }
    }

    for (i = 0; i < Nfe - 2; ++i) {
        // block upper right corner
        coeff = -1. / (wFE(i + 1) * ::sqrt((wFE(i) + wFE(i + 1))*(wFE(i + 1) + wFE(i + 2))));
        mm = indFEDVR(i, Ngs-1, Ngs); nn = indFEDVR(i + 1, Ngs-1, Ngs);
        s = coeff * block(0, Ngs - 1);
        D2.push(s, mm, nn); D2.push(s, nn, mm);
    }

    D2.sort_r();
}

// bounds: FE boundaries, size = Nfe + 1
// Ngs: grid points per finite element (including boundaries)
// `x`, `w` are the global grid points and weights
// `u` is the maximum value of each basis
void D2_matrix(McooDoub_O D2, VecDoub_O x, VecDoub_O w, VecDoub_O u, VecDoub_I bounds, Long_I Ngs)
{
    Long Nfe = bounds.size() - 1; // number of finite elements
#ifdef SLS_CHECK_SHAPE
    Long Nx = Nfe * (Ngs - 1) - 1;
    if (x.size() != Nx || w.size() != Nx || u.size() != Nx ||
        D2.n1() != Nx || D2.n1() != Nx)
        SLS_ERR("wrong shape!");
#endif

    // grid points, weights, base function values in [-1, 1]
    VecDoub x0(Ngs), w0(Ngs), f0(Ngs);
    GaussLobatto(x0, w0);
    pow(f0, w0, -0.5);
    CmatDoub df(Ngs, Ngs); // df(i, j) = f_j(x_i)
    legendre_interp_der(df, x0);
    for (Long i = 0; i < Ngs; ++i)
        for (Long j = 0; j < Ngs; ++j)
            df(j, i) *= f0(i);

    // midpoints, widths and bounds of finite elements
    VecDoub xFE(Nfe), wFE(Nfe);

    for (Long i = 0; i < Nfe; ++i) {
        wFE(i) = 0.5*(bounds(i + 1) - bounds(i));
        xFE(i) = 0.5*(bounds(i) + bounds(i + 1));
    }

    FEDVR_grid(x, w, wFE, xFE, x0, w0);
    pow(u, w, -0.5);

    // Sparse Hamiltonian
    D2_matrix(D2, w0, wFE, df);
}
} // namespace slisc
