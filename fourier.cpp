#include "fourier.h"

// if isign = 1, replaces data[0..2*n-1] by its ifft(), exponent is exp(ikx) .
// if isign = -1, replaces data[0..2*n-1] by n times its fft, exponent is exp(-ikx).
// data is a complex array of length n stored as a real array of length 2*n.
// n must be an integer power of 2.
void four1(Doub *data, Int_I n, Int_I isign) {
	Int nn,mmax,m,j,istep,i;
	Doub wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
	if (n<2 || n&(n-1)) error("n must be power of 2 in four1")
	bit_inv((Comp*)data, n);
	nn = n << 1;
	mmax=2;
	while (nn > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=nn;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void fft(VecComp_IO &data)
{ four1((Doub*)data.ptr(), data.size(), -1); }

void ifft(VecComp_IO &data)
{ four1((Doub*)data.ptr(), data.size(), 1); }

// fft for each column of matrix
// not optimized, very slow
void fft(MatComp_IO &data)
{
	Long i, j, m{ data.nrows() }, n{ data.ncols() };
	VecComp column(m);
	for (j = 0; j < n; ++j) {
		for (i = 0; i < m; ++i)
			column[i] = data[i][j];
		fft(column);
		for (i = 0; i < m; ++i)
			data[i][j] = column[i];
	}
}

void ifft(MatComp_IO &data)
{
	Long i, j, m{ data.nrows() }, n{ data.ncols() };
	VecComp column(m);
	for (j = 0; j < n; ++j) {
		for (i = 0; i < m; ++i)
			column[i] = data[i][j];
		ifft(column);
		for (i = 0; i < m; ++i)
			data[i][j] = column[i];
	}
}

// internal function, don't use
void four1(VecDoub_IO &data, const Int isign)
{ four1(&data[0], data.size() / 2, isign); }

void realft(VecDoub_IO &data, Int_I isign) {
	Int i,i1,i2,i3,i4,n=data.size();
	Doub c1=0.5,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp;
	Doub theta=3.141592653589793238/Doub(n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	for (i=1;i<(n>>2);i++) {
		i2=1+(i1=i+i);
		i4=1+(i3=n-i1);
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r= -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4]= -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[0] = (h1r=data[0])+data[1];
		data[1] = h1r-data[1];
	} else {
		data[0]=c1*((h1r=data[0])+data[1]);
		data[1]=c1*(h1r-data[1]);
		four1(data,-1);
	}
}

void sinft(VecDoub_IO &y) {
	Int j,n=y.size();
	Doub sum,y1,y2,theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;
	theta=3.141592653589793238/Doub(n);
	wtemp=sin(0.5*theta);
	wpr= -2.0*wtemp*wtemp;
	wpi=sin(theta);
	y[0]=0.0;
	for (j=1;j<(n>>1)+1;j++) {
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
		y1=wi*(y[j]+y[n-j]);
		y2=0.5*(y[j]-y[n-j]);
		y[j]=y1+y2;
		y[n-j]=y1-y2;
	}
	realft(y,1);
	y[0]*=0.5;
	sum=y[1]=0.0;
	for (j=0;j<n-1;j+=2) {
		sum += y[j];
		y[j]=y[j+1];
		y[j+1]=sum;
	}
}

void cosft1(VecDoub_IO &y) {
	const Doub PI=3.141592653589793238;
	Int j,n=y.size()-1;
	Doub sum,y1,y2,theta,wi=0.0,wpi,wpr,wr=1.0,wtemp;
	VecDoub yy(n);
	theta=PI/n;
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	sum=0.5*(y[0]-y[n]);
	yy[0]=0.5*(y[0]+y[n]);
	for (j=1;j<n/2;j++) {
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
		y1=0.5*(y[j]+y[n-j]);
		y2=(y[j]-y[n-j]);
		yy[j]=y1-wi*y2;
		yy[n-j]=y1+wi*y2;
		sum += wr*y2;
	}
	yy[n/2]=y[n/2];
	realft(yy,1);
	for (j=0;j<n;j++) y[j]=yy[j];
	y[n]=y[1];
	y[1]=sum;
	for (j=3;j<n;j+=2) {
		sum += y[j];
		y[j]=sum;
	}
}

void cosft2(VecDoub_IO &y, Int_I isign) {
	const Doub PI=3.141592653589793238;
	Int i,n=y.size();
	Doub sum,sum1,y1,y2,ytemp,theta,wi=0.0,wi1,wpi,wpr,wr=1.0,wr1,wtemp;
	theta=0.5*PI/n;
	wr1=cos(theta);
	wi1=sin(theta);
	wpr = -2.0*wi1*wi1;
	wpi=sin(2.0*theta);
	if (isign == 1) {
		for (i=0;i<n/2;i++) {
			y1=0.5*(y[i]+y[n-1-i]);
			y2=wi1*(y[i]-y[n-1-i]);
			y[i]=y1+y2;
			y[n-1-i]=y1-y2;
			wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
			wi1=wi1*wpr+wtemp*wpi+wi1;
		}
		realft(y,1);
		for (i=2;i<n;i+=2) {
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
			y1=y[i]*wr-y[i+1]*wi;
			y2=y[i+1]*wr+y[i]*wi;
			y[i]=y1;
			y[i+1]=y2;
		}
		sum=0.5*y[1];
		for (i=n-1;i>0;i-=2) {
			sum1=sum;
			sum += y[i];
			y[i]=sum1;
		}
	} else if (isign == -1) {
		ytemp=y[n-1];
		for (i=n-1;i>2;i-=2)
			y[i]=y[i-2]-y[i];
		y[1]=2.0*ytemp;
		for (i=2;i<n;i+=2) {
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
			y1=y[i]*wr+y[i+1]*wi;
			y2=y[i+1]*wr-y[i]*wi;
			y[i]=y1;
			y[i+1]=y2;
		}
		realft(y,-1);
		for (i=0;i<n/2;i++) {
			y1=y[i]+y[n-1-i];
			y2=(0.5/wi1)*(y[i]-y[n-1-i]);
			y[i]=0.5*(y1+y2);
			y[n-1-i]=0.5*(y1-y2);
			wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
			wi1=wi1*wpr+wtemp*wpi+wi1;
		}
	}
}

// interpolation using sampling theorem (sinc interp)
// not recommended, use ifft2x(fft()) instead
// assuming y.size() is double
// x must be linspaced.
Comp fft_interp(Doub_I x1, VecDoub_I &x, VecComp_I &y)
{
	Comp sum{};
	Doub dx{ (x.end() - x[0]) / (x.size() - 1.) }, a{ PI/dx };
	Long i, N{ y.size() };
	for (i = 0; i < N; ++i)
		sum += y[i] * sinc(a*(x1 - x[i]));
	return sum;
}

// vector version of fft_interp()
void fft_interp(VecComp_O &y1, VecDoub_I &x1, VecDoub_I &x, VecComp_I &y)
{
	Doub dx{ (x.end() - x[0]) / (x.size() - 1.) }, a{ PI / dx };
	Long i, j, N{ y.size() }, N1{ x1.size() };
	y1.resize(x1);
	for (j = 0; j < N1; ++j)
		for (i = 0; i < N; ++i)
			y1[j] += y[i] * sinc(a*(x1[j] - x[i]));
}

// optimized double zero padding: four1([data, 0, 0, data])
void four2x(Doub *data2, const Doub *data, Int_I n, Int_I isign) {
	Int nn,mmax,m,j,istep,i;
	Doub wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
	if (n<2 || n&(n-1)) error("n must be power of 2 in four1")
	nn = n << 1;
	// get bit inverse order to the end of data2
	Doub *pwork = data2+nn;
	bit_inv((Comp*)pwork, (Comp*)data, n);
	// do two element dft
	Doub *p = data2;
	for (i = 0; i < n>>1; ++i) {
		p[2] = p[0] = pwork[0];
		p[3] = p[1] = pwork[1];
		p[6] = -(p[4] = pwork[2]);
		p[7] = -(p[5] = pwork[3]);
		pwork += 4; p += 8;
	}
	// do the rest
	nn <<= 1;
	mmax=4;
	while (nn > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=nn;i+=istep) {
				j=i+mmax;
				tempr=wr*data2[j-1]-wi*data2[j];
				tempi=wr*data2[j]+wi*data2[j-1];
				data2[j-1]=data2[i-1]-tempr;
				data2[j]=data2[i]-tempi;
				data2[i-1] += tempr;
				data2[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void fft2x(VecComp_O &data2, VecComp_I &data)
{
	data2.resize(data.size()*2);
	four2x((Doub*)data2.ptr(), (Doub*)data.ptr(), data.size(), -1);
}

void ifft2x(VecComp_O &data2, VecComp_I &data)
{
	data2.resize(data.size()*2);
	four2x((Doub*)data2.ptr(), (Doub*)data.ptr(), data.size(), 1);
}

// optimized 4x zero padding: four1([data, 0, 0 ,0, 0, 0, 0, data]);
void four4x(Doub *data2, const Doub *data, Int_I n, Int_I isign) {
	Int nn,mmax,m,j,istep,i;
	Doub wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
	if (n<2 || n&(n-1)) error("n must be power of 2 in four1")
	nn = n << 1;
	// get bit inverse order to the end of data2
	Doub *pwork = data2+3*nn;
	bit_inv((Comp*)pwork, (Comp*)data, n);
	// do two element dft
	Doub *p = data2, temp;
	if (isign < 0)
		for (i = 0; i < n>>1; ++i) {
			p[6] = p[4] = p[2] = p[0] = pwork[0];
			p[7] = p[5] = p[3] = p[1] = pwork[1];
			p += 8; temp = pwork[3];
			p[7] = p[4] = -(p[3] = p[0] = pwork[2]);
			p[5] = p[2] = -(p[6] = p[1] = temp);
			pwork += 4; p += 8;
		}
	else
		for (i = 0; i < n>>1; ++i) {
			p[6] = p[4] = p[2] = p[0] = pwork[0];
			p[7] = p[5] = p[3] = p[1] = pwork[1];
			p += 8; temp = pwork[3];
			p[4] = p[3] = -(p[7] = p[0] = pwork[2]);
			p[6] = p[5] = -(p[2] = p[1] = temp);
			pwork += 4; p += 8;
		}
	// do the rest
	nn <<= 2;
	mmax=8;
	while (nn > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=nn;i+=istep) {
				j=i+mmax;
				tempr=wr*data2[j-1]-wi*data2[j];
				tempi=wr*data2[j]+wi*data2[j-1];
				data2[j-1]=data2[i-1]-tempr;
				data2[j]=data2[i]-tempi;
				data2[i-1] += tempr;
				data2[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void fft4x(VecComp_O &data4, VecComp_I &data)
{
	data4.resize(data.size()*4);
	four4x((Doub*)data4.ptr(), (Doub*)data.ptr(), data.size(), -1);
}

void ifft4x(VecComp_O &data4, VecComp_I &data)
{
	data4.resize(data.size()*4);
	four4x((Doub*)data4.ptr(), (Doub*)data.ptr(), data.size(), 1);
}

void dft(MatComp_O &Y, Doub kmin, Doub kmax, Long_I Nk, MatComp_I &X, Doub xmin, Doub xmax)
{
	Long i, j, k, Nx = X.nrows(), Nc = X.ncols();
	Doub dk = (kmax - kmin) / (Nk - 1), dx = (xmax - xmin) / (Nx - 1);
	const Comp *pxi;
	Comp *pyj, factor, expo, dexpo;
	Y.resize(Nk, Nc); Y = 0.;
	for (j = 0; j < Nk; ++j) {
		pyj = Y[j];
		expo = exp(Comp(0, -(kmin + dk*j)*(xmin - dx)));
		dexpo = exp(Comp(0, -(kmin + dk*j)*dx));
		for (i = 0; i < Nx; ++i) {
			pxi = X[i];
			expo *= dexpo;
			for (k = 0; k < Nc; ++k)
				pyj[k] += expo*pxi[k];
		}
	}
}

// parallel version
void dft_par(MatComp_O &Y, Doub kmin, Doub kmax, Long_I Nk, MatComp_I &X, Doub xmin, Doub xmax)
{
	Long j, Nx = X.nrows(), Nc = X.ncols();
	Doub dk = (kmax - kmin) / (Nk - 1), dx = (xmax - xmin) / (Nx - 1);
	Y.resize(Nk, Nc); Y = 0.;
#pragma omp parallel for
	for (j = 0; j < Nk; ++j) {
		Long i, k;
		const Comp *pxi;
		Comp *pyj, factor, expo, dexpo;
		pyj = Y[j];
		expo = exp(Comp(0, -(kmin + dk*j)*(xmin - dx)));
		dexpo = exp(Comp(0, -(kmin + dk*j)*dx));
		for (i = 0; i < Nx; ++i) {
			pxi = X[i];
			expo *= dexpo;
			for (k = 0; k < Nc; ++k)
				pyj[k] += expo*pxi[k];
		}
	}
}

void idft(MatComp_O &X, Doub xmin, Doub xmax, Long_I Nx, MatComp_I &Y, Doub kmin, Doub kmax)
{
	Long i, j, k, Nk = Y.nrows(), Nc = Y.ncols();
	Doub dk = (kmax - kmin) / (Nk - 1), dx = (xmax - xmin) / (Nx - 1);
	const Comp *pyi;
	Comp *pxj, factor, expo, dexpo;
	X.resize(Nx, Nc); X = 0.;
	for (j = 0; j < Nx; ++j) {
		pxj = X[j];
		expo = exp(Comp(0, (xmin + dx*j)*(kmin - dk)));
		dexpo = exp(Comp(0, (xmin + dx*j)*dk));
		for (i = 0; i < Nk; ++i) {
			pyi = Y[i];
			expo *= dexpo;
			for (k = 0; k < Nc; ++k)
				pxj[k] += expo*pyi[k];
		}
	}
}

void idft_par(MatComp_O &X, Doub xmin, Doub xmax, Long_I Nx, MatComp_I &Y, Doub kmin, Doub kmax)
{
	Long j, Nk = Y.nrows(), Nc = Y.ncols();
	Doub dk = (kmax - kmin) / (Nk - 1), dx = (xmax - xmin) / (Nx - 1);
	X.resize(Nx, Nc); X = 0.;
#pragma omp parallel for
	for (j = 0; j < Nx; ++j) {
		Long i, k;
		const Comp *pyi;
		Comp *pxj, factor, expo, dexpo;
		pxj = X[j];
		expo = exp(Comp(0, (xmin + dx*j)*(kmin - dk)));
		dexpo = exp(Comp(0, (xmin + dx*j)*dk));
		for (i = 0; i < Nk; ++i) {
			pyi = Y[i];
			expo *= dexpo;
			for (k = 0; k < Nc; ++k)
				pxj[k] += expo*pyi[k];
		}
	}
}
