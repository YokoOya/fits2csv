#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <ctime>
#include "fitsio.h"

const int lstr = 256;
const char DEG[] = "deg";
const double INF = 1e+308;


using namespace std;

void printerror(int status) {
	if (status) {
		fits_report_error(stderr, status);
		exit(status);
	}
	return;
}


void reorderImage2(double *buffin, double *buffout, long *naxis) {
	for (int i = 0; i < naxis[0]; ++i) for (int j = 0; j < naxis[1]; ++j) buffout[i * naxis[1] + j] = buffin[i + j * naxis[0]];
	return ;
}

void reorderImage3(double *buffin, double *buffout, long *naxis) {
	for (int i = 0; i < naxis[0]; ++i) for (int j = 0; j < naxis[1]; ++j) for (int k = 0; k < naxis[2]; ++k) buffout[(i * naxis[1] + j) * naxis[2] + k] = buffin[(k * naxis[1] + j) * naxis[0] + i];
	return ;
}

void reorderImage4(double *buffin, double *buffout, long *naxis) {
	for (int i = 0; i < naxis[0]; ++i) for (int j = 0; j < naxis[1]; ++j) for (int k = 0; k < naxis[2]; ++k) for (int l = 0; l < naxis[3]; ++l) buffout[((i * naxis[1] + j) * naxis[2] + k) * naxis[3] + l] = buffin[((naxis[2] * l + k) * naxis[1] + j) * naxis[0] + i];
	return ;
}


void readheader(fitsfile *fptrin, FILE *fptrout, int hdu, int *statusptr) {
	
	int nkeys, keypos, hdutype;
	char card[FLEN_CARD];
	
	if (fits_movabs_hdu(fptrin, hdu, &hdutype, statusptr)) printerror(*statusptr);
	if (fits_get_hdrpos(fptrin, &nkeys, &keypos, statusptr)) printerror(*statusptr);
	
	for (int ihead = 0; ihead < nkeys; ++ihead) {
		if (fits_read_record(fptrin, ihead, card, statusptr)) printerror(*statusptr);
		fprintf(fptrout, "#%s\n", card);
	}
	
	return ;
}




void readimage(fitsfile *fptrin, FILE *fptrout, int hdu, int *statusptr) {
	
	int hdutype, nfound, anynull;
	char dummy_char[lstr], comment[lstr], btype[lstr], bunit[lstr];
	long naxis, npixel = 1, fpixel = 1;
	double nullval = 0., bscale, bzero, peakIntensity = -INF;
	
	if (fits_movabs_hdu(fptrin, hdu, &hdutype, statusptr)) printerror(*statusptr);
	
	
	if (fits_read_key(fptrin, TLONG, "NAXIS", &naxis, comment, statusptr)) printerror(*statusptr);
	
	long naxes[naxis], iaxes[naxis];
	double crpix[naxis], crval[naxis], cdelt[naxis];
	char cunit[naxis][lstr];
	if (fits_read_keys_lng(fptrin, "NAXIS", 1, naxis, naxes, &nfound, statusptr)) printerror(*statusptr);
	
	for (long i = 0; i < naxis; ++i) {
		sprintf(dummy_char, "CUNIT%ld", i + 1);
		if (fits_read_key(fptrin, TSTRING, dummy_char, &cunit[i], comment, statusptr)) printerror(*statusptr);
		sprintf(dummy_char, "CRPIX%ld", i + 1);
		if (fits_read_key(fptrin, TDOUBLE, dummy_char, &crpix[i], comment, statusptr)) printerror(*statusptr);
		sprintf(dummy_char, "CRVAL%ld", i + 1);
		if (fits_read_key(fptrin, TDOUBLE, dummy_char, &crval[i], comment, statusptr)) printerror(*statusptr);
		sprintf(dummy_char, "CDELT%ld", i + 1);
		if (fits_read_key(fptrin, TDOUBLE, dummy_char, &cdelt[i], comment, statusptr)) printerror(*statusptr);
		
		--crpix[i];
		
		for (int j = 0; j < int(strlen(DEG)); ++j) {
			if (cunit[i][j] != DEG[j]) break;
			if (j == strlen(DEG) - 1) {
				sprintf(cunit[i], "ARCSEC  ");
				crval[i] *= 3600.;
				cdelt[i] *= 3600.;
			}
		}
	}
	
	
	if(fits_read_key(fptrin, TDOUBLE, "BSCALE", &bscale, comment, statusptr)) printerror(*statusptr);
	if(fits_read_key(fptrin, TDOUBLE, "BZERO", &bzero, comment, statusptr)) printerror(*statusptr);
	if(fits_read_key(fptrin, TSTRING, "BTYPE", &btype, comment, statusptr)) printerror(*statusptr);
	if(fits_read_key(fptrin, TSTRING, "BUNIT", &bunit, comment, statusptr)) printerror(*statusptr);

	
	
	for (int i = 0; i < naxis; ++i) npixel *= naxes[i];
	
	
	double *buffin, *buffout;
	buffin = (double*) malloc(sizeof(double) * npixel);
	buffout = (double*) malloc(sizeof(double) * npixel);

	
	memset(buffin, 0., sizeof(double) * npixel);
	
	
	
	if (fits_read_img(fptrin, TDOUBLE, fpixel, npixel, &nullval, buffin, &anynull, statusptr)) printerror(*statusptr);
	
	if (naxis == 2) {reorderImage2(buffin, buffout, naxes);}
	else if (naxis == 3) {reorderImage3(buffin, buffout, naxes);}
	else if (naxis == 4) {reorderImage4(buffin, buffout, naxes);}
    else if (naxis < 2) {printf("\n\n!! ERROR :: naxis > 2.  The values for intensity will not be copied correctly !!\n\n");}
    else {printf("\n\n!! ERROR :: naxis > 4.  The values for intensity will not be copied correctly !!\n\n");}

	
	
	for (long i = 0; i < naxis; ++i) fprintf(fptrout, ",AXIS%ld id", i + 1);
	fprintf(fptrout, ",%s\n", btype);
	
	memset(iaxes, 0, sizeof(long) * naxis);
    
    for (long i = 0; i < npixel; ++i) {
        if (peakIntensity < buffout[i]) peakIntensity = buffout[i];
    }
	
    printf("\n\tPeak intensity = %le\t\t-> Intensities are normalized by this value.", peakIntensity);
    
	for (long i = 0; i < npixel; ++i) {
        fprintf(fptrout, "%ld", i);
		for (long j = 0; j < naxis; ++j) {
			fprintf(fptrout, ",%ld", iaxes[j]);
		}
		
		fprintf(fptrout, ",%le\n", buffout[i] / peakIntensity);
	
		++iaxes[naxis - 1];
		for (long j = naxis - 1; j > 0; --j) if (iaxes[j] >= naxes[j]) {
				iaxes[j] -= naxes[j];
				++iaxes[j - 1];
		}
	}
	
	
	
	free(buffin);
	free(buffout);

	
	return ;
}



int main() {
	
	puts("\n----------------------");
	clock_t start, end;
	start = clock();
	//
	
	
	int status = 0, hdu = 1, nfitsfile;
	char infitsname[lstr], outcsvname[lstr];
	FILE *outcsvfile;
	fitsfile *infitsfile;
	
	scanf("%d ", &nfitsfile);
	
	
	for (int i = 0; i < nfitsfile; ++i) {
		
		printf("---ReadFits (%d / %d)---\n", i + 1, nfitsfile);

		scanf("%s ", infitsname);
		printf("infitsname = %s\n", infitsname);
		sprintf(outcsvname, "%s.normalized.csv", infitsname);
		
		
		
		if (fits_open_file(&infitsfile, infitsname, READONLY, &status)) printerror(status);
		outcsvfile = fopen(outcsvname, "w");
		
		
		
		readimage(infitsfile, outcsvfile, hdu, &status);
		
		
		
		
		if (fits_close_file(infitsfile, &status)) printerror(status);
		fclose(outcsvfile);
		printf("\n\nSaved: '%s'\n\n", outcsvname);
		
		end = clock();
		printf("Duration time: %lf [sec]\n\n", double(end - start) / CLOCKS_PER_SEC);
	}
	
	end = clock();
	puts("\n----------------------\n");
	printf("%d fits files are converted to csv files.\n", nfitsfile);
	printf("Duration time: %lf [sec]\n", double(end - start) / CLOCKS_PER_SEC);
	puts("Done.\n\n");
	
	puts("\007");
}


