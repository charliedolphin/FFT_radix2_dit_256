// ditradix2_256
// 2-256 point FFT radix 2 (Decimation in Time) Algorithm
// 25.12.2015 by Dmitry Chistyakov 
//
// Updates:
//      25.12.2015 - release.
//      01.04.2017 - optimized binary-inversion function, reduced calculation time.
//  Time of calculation (between Frame forming stage and calculation amplitudes)
//  is 2.2 ms (64 point, clock = 160MHz (RISC)).
//      12.04.2017 - reduced number of math operations.
//*****************************************************************************
//*****************************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

// To set frame size (2-256), replace size of arrays and N
float fft_r[64]; 
float fft_j[64]; 
float exp_r[64]; 
float exp_j[64]; 

unsigned short int N = 64;

unsigned short int k, n, step;  
unsigned char u;  
unsigned short int i = 0;  

// Intermediate results
float r_snk, j_snk, r_sns2k, j_sns2k, number, rimd, jimd;  

// Indexes
unsigned short int index0, index1;

// Binary-inversion function for sorting data into a frame
unsigned char reversebits (unsigned short int value) { 
    unsigned char n, ss, ss2;                                        
    unsigned char bith, bitl;                                                                                      
    for (n=0; n<4; n++) {

        bith = value >> 7-n;
        bith = bith << 7;
        bith = bith >> 7;
        bitl = value << 7-n;
        bitl = bitl >> 7;

        switch (n){
            case 0:
                ss=1;
                ss2=128;
                break;
            case 1:
                ss=2;
                ss2=64;
                break;
            case 2:
                ss=4;
                ss2=32;
                break;
            case 3:
                ss=8;
                ss2=16;
                break;
            }

            if (bith > bitl) {
                value=value-ss2;
                value=value+ss;
            }
            if (bith < bitl) {
                value=value+ss2;
                value=value-ss;
            }
    }
    return value; 
}

// Return REAL part of multiplying
// x1, y1 should be float, else missing accuracy 
float m_real (float x1, float y1, float x2, float y2) { 
    return x1*x2-y1*y2; 
}
// Return IM part of multiplying
float m_im (float x1, float y1, float x2, float y2) { 
    return x1*y2+y1*x2; 
} 
// Function for getting table of complex exponents
void FFT_ExpCalculation(unsigned short int N){
    float arg;
    float pi = 3.1415926535897932384626553633832795;
    unsigned short int k, step;
    unsigned short int i = 0;
    for (step=2; step<=N; step*=2) { 
        for (k=0; k<step/2; k++) {
            arg = 2*pi*k/step;
            exp_r[i]=cos(arg); 
            exp_j[i]=sin(-arg);
            i++;  
       }
    }
}

int main() { 
    // Get table of exponents
    FFT_ExpCalculation(N);

    // Frame forming stage (getting data for analysis)
    for (i=0; i<N; i++){
        fft_r[i]=i; // Put sample into array
    }
    // Copying frame
    for (i=0; i<N; i++) { 
        number=fft_r[i]; 
        fft_j[i]=number; 
    }
    // Applying binary-inversion function
    for (i=0; i<N; i++) { 
        u=reversebits(i); 
        fft_r[i]=fft_j[u*N/256]; 
    }
    // Clearing the copy
    for (i=0; i<N; i++) { 
        fft_j[i]=0; 
    } 

    i=0;
    for (step=2; step<=N; step*=2) {
        for (k=0; k<step/2; k++)  {
            for (n=0; n<N/step; n++) {
                // Perform butterly's
                index0 = step*n+k;
                index1 = (step*n)+(step/2)+k;

                rimd = m_real(exp_r[i], exp_j[i], 
                    fft_r[index1], fft_j[index1]);
                jimd = m_im(exp_r[i], exp_j[i], 
                    fft_r[index1], fft_j[index1]);

                r_snk=fft_r[index0]+rimd;
                j_snk=fft_j[index0]+jimd;
                r_sns2k=fft_r[index0]-rimd;
                j_sns2k=fft_j[index0]-jimd;
                
                fft_r[index0]=r_snk; 
                fft_j[index0]=j_snk;
                fft_r[index1]=r_sns2k;                      
                fft_j[index1]=j_sns2k; 
            } 
            i++; 
        } 
    }
    // FFT end.
    // Calculate amplitude spectrum
    for (i=0; i<N; i++) {
        fft_r[i]=sqrt(fft_r[i]*fft_r[i]+fft_j[i]*fft_j[i]);
    }  

    // See results
    for (i=0; i<N; i++) {
        printf("    %f\n", fft_r[i]);
    }  
    return 0;
}
