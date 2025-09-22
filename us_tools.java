/**
V. 1.0
Ultrasound tools (us_tools) organized as java methods for calling using MATLAB.

Syntax in MATLAB:

	javaaddpath(pwd);
	Obj = us_tools;
	Result = javaMethod('method_name',Obj,arg1,arg2,...,argN);

Author: João Henrique Uliana
Institution: Department of Physics, FFCLRP-USP, Ribeirão Preto - SP, Brazil
Created: 08-2023
Updated: 10-2023

*/

import java.util.*;
import java.lang.Math;

public class us_tools{
//------------------------------------------Delay Map-------------------------------------------------------
	
	//Delay map for conventional US
	public static double[][] delays (int Z, int X, double[]prm){
		
		/**--------------------------------------------------
		Calling in matlab:
			
			Z = size(rfdata,1);
			X = size(rfdata,2);
			prm(1) = dx; %[mm]
			prm(2) = c; %[m/s]
			prm(3) = sf; %[Hz]
		
			delays = javaMethod('delays',Obj,Z,X,prm);
			
		----------------------------------------------------*/
		
		//parameters
		double dx = prm[0]; //transducer's element pitch[mm]
		double c = prm[1]; //speed of sound [m/s]
		double sf = prm[2]; //sampling frequency [Hz]		
		
		double[][] delays = new double[Z][2*X]; //delays with dimensions [Z,2X]
		
		double dd = c*(1/sf)*1000; //[mm]
		double dz = dd/2;
		
		double z_st;
		double x_pos;
		double d;
		for(int z = 0; z < Z; z++) {
			z_st = z*dz;
			for(int x = 0; x < 2*X; x++) {
				x_pos = (X-x)*dx;
				d = Math.sqrt(z_st*z_st + x_pos*x_pos) - z_st;
				delays[z][x] = (int)(d/dd); //delays in pts
			}
		}
		return delays;
	}
	
	//Delay map for photoacoustic image
	public static double[][] delaysPA (int Z, int X, double[]prm){
		
		/**--------------------------------------------------
		Calling in matlab:
			
			Z = size(rfdata,1);
			X = size(rfdata,2);
			prm(1) = dx; %[mm]
			prm(2) = c; %[m/s]
			prm(3) = sf; %[Hz]
		
			delays = javaMethod('delaysPA',Obj,Z,X,prm);
			
		----------------------------------------------------*/
		
		//parameters
		double dx = prm[0]; //transducer's element pitch[mm]
		double c = prm[1]; //speed of sound [mm]
		double sf = prm[2]; //sampling frequency [mm]
		
		double[][] delaysPA = new double[Z][2*X]; //delays with dimensions [Z,2X]
		
		double dz = c*(1/sf)*1000; //[mm]
		
		double z_st;
		double x_pos;
		double d;
		for(int z = 0; z < Z; z++) {
			z_st = z*dz;
			for(int x = 0; x < 2*X; x++) {
				x_pos = (X-x)*dx;
				delaysPA[z][x] = (int)((Math.sqrt(z_st*z_st + x_pos*x_pos) - z_st)/dz);//delays in pts
			}
		}
		return delaysPA;
	}


//------------------------------------Decimated Beamformers-------------------------------------------------

	// Delay and Sum (das) with decimation of factor n
	public static double[][] das (double[][] rfdata, double[][]delays, double[]prm, int n){
		
		/**--------------------------------------------------
		Calling in matlab:
		
			prm(1) = dz; %[mm]
			prm(2) = dx; %[mm]
			prm(3) = aperture; %[channels]
			prm(4) = cut_angle; %[°]
		
			dasdata = javaMethod('das',Obj,rfdata,delays,prm, n);
			
		----------------------------------------------------*/
		
		//parameters
		double dz = prm[0]; //axial resolution [mm]
		double dx = prm[1]; //lateral resolution [mm]
		int aperture = (int)prm[2]; //aperture [channels]
		double ctag = prm[3]; //cut angle [°]
		
		//Matrix (image) dimensions [px]
		int X = rfdata[1].length; //Lateral
		int Z = rfdata.length; //Axial
		int Z2 = (int)(Z/n); //decimation in Z (n)
		
		double[][] das = new double[Z2][X];
		double[][] aux_das = new double[Z2+10][X];
		
		//ZeroPad for depth
		int Zeroadd = (int)(Math.sqrt((Z*dz)*(Z*dz) + (X*dx)*(X*dx))/dz); //maxdepth
		double[][] rfpad = new double[Zeroadd][X];
		for (double[] row: rfpad){
			Arrays.fill(row, 0.0);
		}
		for(int z = 0; z < Z; z++){
			for(int x = 0; x < X; x++) {
				rfpad[z][x] = rfdata[z][x];
			}
		}
		
		//Directivity of element
		double acut = (double)(dz/dx)*Math.tan((ctag*Math.PI)/180.0);
		int zlim = (int)(aperture/acut);
		
		//Recon variables
		int hap = (int)(aperture/2);
		int realap = 0;
		int i_d = 0;
		int z = 0;
		int m = 0;
		
		//Delay and sum
		while(z < Z) {
			if(z < zlim)
				realap = (int)(((acut/2)*(z))+1);
			else
				realap = hap;
			for(int x = 0; x < X; x++) {
				double count = 0;
				for (int i = (x - realap); i <= (x + realap); i++){
					if(i >= 0 && i < X){
						i_d = z + (int)delays[z][((X-1)-x)+i];
						aux_das[m][x] = aux_das[m][x] + rfpad[i_d][i];
						count = count +1;
					}				
				}
			aux_das[m][x] = aux_das[m][x]/count;
            }
			m = m+1;
			z = z+n;
		}
		for(z = 0; z < Z2; z++) {
			for(int x = 0; x < X; x++) {
				das[z][x] = aux_das[z][x];
			}
		}			
		return das;
	}
	
	
	// Delay Multiply and Sum (dmas) with decimation of n
	public static double[][] dmas (double[][] rfdata, double[][]delays, double[]prm, int n){
		
		/**--------------------------------------------------
		Calling in matlab:
		
			prm(1) = dz; %[mm]
			prm(2) = dx; %[mm]
			prm(3) = aperture; %[channels]
			prm(4) = cut_angle; %[°]
		
			dasdata = javaMethod('dmas',Obj,rfdata,delays,prm, n);
			
		----------------------------------------------------*/
		
		//parameters
		double dz = prm[0]; //axial resolution [mm]
		double dx = prm[1]; //lateral resolution [mm]
		int aperture = (int)prm[2]; //aperture [channels]
		double ctag = prm[3]; //cut angle [°]
		
		//dimensions (px)
		int X = rfdata[1].length; //Lateral
		int Z = rfdata.length; //Axial
		int Z2 = (int)(Z/n); //decimation in Z (n)
		
		double[][] dmas = new double[Z2][X];
		double[][] aux_das = new double[Z2+10][X];
		
		//ZeroPad
		int Zeroadd = (int)(Math.sqrt((Z*dz)*(Z*dz) + (X*dx)*(X*dx))/dz); //maxdepth
		double[][] rfpad = new double[Zeroadd][X];
		for (double[] row: rfpad){
			Arrays.fill(row, 0.0);
		}
		for(int z = 0; z < Z; z++){
			for(int x = 0; x < X; x++) {
				rfpad[z][x] = rfdata[z][x];
			}
		}
		
		//Directivity of element
		double acut = (double)(dz/dx)*Math.tan((ctag*Math.PI)/180.0);
		int zlim = (int)(aperture/acut);
		
		//Recon variables
		int hap = (int)(aperture/2);
		int realap = 0;
		int i_d = 0;
		int i_dj = 0;
		int z = 0;
		int m = 0;
		
		//Delay Multiply and Sum
		while(z < Z) {
			if(z < zlim)
				realap = (int)(((acut/2)*(z))+1);
			else
				realap = hap;
			for(int x = 0; x < X; x++) {
				double count = 1;
				for (int i = (x - realap); i <= (x + realap); i++){
					if(i >= 0 && i < (X-1)){
						i_d = z + (int)delays[z][((X-1)-x)+i];
						for (int j = (i+1); j <= (x + realap); j++){
							if(j < X){
								i_dj = z + (int)delays[z][((X-1)-x)+j];
								aux_das[m][x] = aux_das[m][x] + rfpad[i_d][i]*rfpad[i_dj][j];
								count = count +1;
							}
						}
					}				
				}
			aux_das[m][x] = aux_das[m][x]/count;
            }
			m = m+1;
			z = z+n;
		}
		for(int z2 = 0; z2 < Z2; z2++){
			for(int x = 0; x < X; x++){
				dmas[z2][x] = aux_das[z2][x];
			}
		}		
		return dmas;
	}
	
	// Filtered Delay Multiply and Sum (fdmas) with decimation of n
	public static double[][] fdmas (double[][] rfdata, double[][]delays, double[]prm, int n){
		
		/**--------------------------------------------------
		Calling in matlab:
		
			prm(1) = dz; %[mm]
			prm(2) = dx; %[mm]
			prm(3) = aperture; %[channels]
			prm(4) = cut_angle; %[°]
		
			dasdata = javaMethod('fdmas',Obj,rfdata,delays,prm, n);
			
		----------------------------------------------------*/
		
		//parameters
		double dz = prm[0]; //axial resolution [mm]
		double dx = prm[1]; //lateral resolution [mm]
		int aperture = (int)prm[2]; //aperture [channels]
		double ctag = prm[3]; //cut angle [°]
		
		//dimensions (px)
		int X = rfdata[1].length; //Lateral
		int Z = rfdata.length; //Axial
		int Z2 = (int)(Z/n); //decimation in Z (n)
		
		double[][] fdmas = new double[Z2][X];
		double[][] aux_das = new double[Z2+10][X];
		
		//ZeroPad
		int Zeroadd = (int)(Math.sqrt((Z*dz)*(Z*dz) + (X*dx)*(X*dx))/dz); //maxDepth
		double[][] rfpad = new double[Zeroadd][X];
		for (double[] row: rfpad){
			Arrays.fill(row, 0.0);
		}
		for(int z = 0; z < Z; z++){
			for(int x = 0; x < X; x++) {
				rfpad[z][x] = rfdata[z][x];
			}
		}
		
		//Directivity of element
		double acut = (double)(dz/dx)*Math.tan((ctag*Math.PI)/180.0);
		int zlim = (int)(aperture/acut);
		
		//recon variables
		int hap = (int)(aperture/2);
		int realap = 0;
		int i_d = 0;
		int i_dj = 0;
		int z = 0;
		int m = 0;
		
		//Filtered Delay Multiply and Sum
		while(z < Z) {
			if(z < zlim)
				realap = (int)(((acut/2)*(z))+1);
			else
				realap = hap;
			for(int x = 0; x < X; x++) {
				double count = 0;
				for (int i = (x - realap); i <= (x + realap); i++){
					if(i >= 0 && i < (X-1)){
						i_d = z + (int)delays[z][((X-1)-x)+i];
						for (int j = (i+1); j <= (x + realap); j++){
							if(j < X){
								i_dj = z + (int)delays[z][((X-1)-x)+j];
								aux_das[m][x] = aux_das[m][x] + Math.signum(rfpad[i_d][i]*rfpad[i_dj][j])
								*Math.sqrt(Math.abs(rfpad[i_d][i]*rfpad[i_dj][j]));
								count = count +1;
							}
						}
					}				
				}
			aux_das[m][x] = aux_das[m][x]/count;
            }
			m = m+1;
			z = z+n;
		}
		for(int z2 = 0; z2 < Z2; z2++){
			for(int x = 0; x < X; x++){
				fdmas[z2][x] = aux_das[z2][x];
			}
		}
		
		return fdmas;
	}
	
	// Short lag spatial coherence with decimation of N
	public static double[][] slsc (double[][] rfdata, double[][]delays, double[]prm, int N){
		
		/**--------------------------------------------------
		Calling in matlab:
		
			prm(1) = dz; %[mm]
			prm(2) = dx; %[mm]
			prm(3) = aperture; %[channels]
			prm(4) = quality factor; %[%]
		
			dasdata = javaMethod('slsc',Obj,rfdata,delays,prm,N);
			
		----------------------------------------------------*/
		
		//parameters
		double dz = prm[0];//axial resolution [mm]
		double dx = prm[1];//lateral resolution [mm]
		int aperture = (int)prm[2];//aperture [channels]
		double Q = prm[3]; //Quality factor (%)
		
		//dimensions (px)
		int X = rfdata[1].length; //Lateral
		int Z = rfdata.length; //Axial
		int Z2 = (int)(Z/N); //decimation in Z (N)
		int M = (int)(aperture*(Q/100)); //M factor = aperture*(Q/100)
		
		double[][] slsc = new double[Z2][X];
		double[][] aux_das = new double[Z2+10][X];
		
		//ZeroPad
		int Zeroadd = (int)(Math.sqrt((Z*dz)*(Z*dz) + (X*dx)*(X*dx))/dz); //MaxDepth
		double[][] rfpad = new double[Zeroadd][X];
		for (double[] row: rfpad){
			Arrays.fill(row, 0.0);
		}
		for(int z = 0; z < Z; z++){
			for(int x = 0; x < X; x++) {
				rfpad[z][x] = rfdata[z][x];
			}
		}
		
		//Recon variables
		int hap = (int)(aperture/2);
		int realap = 0;
		int i_d = 0;
		int i_m = 0;
		int z = 0;
		int n = 0;
		
		//SLSC
		double s1,s2,s3,s4,slags,count;
		while(z < Z) {
			for(int x = 0; x < X; x++) {
				slags = 0;
				for(int m = 1; m <= M; m++) {
					s1 = 0;
					s2 = 0;
					s3 = 0;
					s4 = 0;
					for (int i = (x - hap); i <= (x + hap); i++){
						if(i >= 0 && i < ((X-1)-m)){
							for (int j = z; j <= z+4; j++){
								i_d = j + (int)delays[z][((X-1)-x)+i];
								i_m = j + (int)delays[z][((X-1)-x)+i+m];
								s1 = s1 + rfpad[i_d][i]*rfpad[i_m][i+m];
								s2 = s2 + rfpad[i_d][i]*rfpad[i_d][i];
								s3 = s3 + rfpad[i_m][i+m]*rfpad[i_m][i+m];
							}
						}
					}				
					if(s2*s3 != 0)
						s4 = s4 + (s1/Math.sqrt(s2*s3));
					else
						s4 = s4 + 0;
					slags = slags + s4/((2*hap)-m);
				}
				aux_das[n][x] = slags;
            }
			n = n+1;
			z = z+N;
		}
		for(int z2 = 0; z2 < Z2; z2++){
			for(int x = 0; x < X; x++){
				slsc[z2][x] = aux_das[z2][x];
			}
		}
		
		return slsc;
	}
	
	// Anguled Plane Wave Delay and Sum (apwdas) for Coherent Plane Wave Compound (CPWC) 
	public static double[][] apwdas (double[][] rfdata, double[]prm, int n){
		//this method does'nt need the delays map
		/**--------------------------------------------------
		Calling in matlab:
		
			prm(1) = dz; %[mm]
			prm(2) = dx; %[mm]
			prm(3) = aperture; %[channels]
			prm(4) = cut_angle; %[°]
			prm(5) = plane wave angle; %[rads]
		
			dasdata = javaMethod('apwdas',Obj,rfdata,prm,n);
			
		----------------------------------------------------*/
		
		
		//parameters
		double dz = prm[0];//axial resolution [mm]
		double dx = prm[1];//lateral resolution [mm]
		int aperture = (int)prm[2];//aperture [channels]
		double ctag = prm[3];//cut angle [°]
		double wag = prm[4]; //plane wave angle [rads]
		
		//dimensions (px)
		int X = rfdata[1].length; //Lateral
		int Z = rfdata.length; //Axial
		int Z2 = (int)(Z/n); //decimation in Z
		
		double[][] apwdas = new double[Z2][X];
		double[][] aux_das = new double[Z2+10][2*X];
		
		//ZeroPad
		int Zeroadd = (int)((Math.sqrt((Z*dz)*(Z*dz) + (X*dx)*(X*dx))/dz)); //maxDepth
		double[][] rfpad = new double[Zeroadd][X];
		for (double[] row: rfpad){
			Arrays.fill(row, 0.0);
		}
		for(int z = 0; z < Z; z++){
			for(int x = 0; x < X; x++) {
				rfpad[z][x] = rfdata[z][x];
			}
		}
		
		//Directivity of element
		double acut = (double)(dz/dx)*Math.tan((ctag*Math.PI)/180.0);
		int zlim = (int)(aperture/acut);
		
		//Recon variables
		int hap = (int)(aperture/2);
		int realap;
		int i_d = 0;
		int z = 0;
		int m = 0;
		int i_x;
		double count;
		double cos_val = Math.cos(wag);
		double sin_val = Math.sin(Math.abs(wag));
		double z_pos;
		double x_pos;
		double d;
		double z_st;
		
		//Anguled Plane Wave Delay and Sum
		while(z < Z) {
			z_st = z*dz;
			if(z < zlim)
				realap = (int)(((acut/2)*z)+1);
			else
				realap = hap;
			for(int x = 0; x < X; x++){
				count = 0;
				if(wag >= 0)
					z_pos = (z_st*cos_val + x*sin_val*dx);
				else
					z_pos = (z_st*cos_val + (X-x)*sin_val*dx);
				for (int i = (x - realap); i <= (x + realap); i++){
					if(i >= 0 && i < X){
						x_pos = (i-x)*dx;
						i_d = z + (int)((Math.sqrt(z_st*z_st + x_pos*x_pos) + z_pos - 2*z_st)/(2*dz));
						aux_das[m][x] = aux_das[m][x] + rfpad[i_d][i];
						count = count + 1;
					}				
				}
			aux_das[m][x] = aux_das[m][x]/(count);
		    }
			m = m+1;
			z = z+n;
		}
		for(int z2 = 0; z2 < Z2; z2++){
			for(int x = 0; x < X; x++){
				apwdas[z2][x] = aux_das[z2][x];
			}
		}		
		return apwdas;
	}

	
//-------------------------------------Image processing-----------------------------------------------------	
	
	//Envelope calculation using quadrature demodulation(IQ)
	public static double[][] iqenv (double[][] rfdata, double sf, double cf){
		
		//parameters
		double wc = cf*2.0*Math.PI; //central frequency [rads]
		
		int X = rfdata[1].length;
		int Z = rfdata.length;
				
		double [] sinw = new double[Z];
		double [] cosw = new double[Z];
		//double [][] Idata = new double [Z][X];
		//double [][] Qdata = new double [Z][X];
		double [][] iqenv = new double [Z][X];
		double [][] inter_iqenv = new double [Z][X];
		double [][] I = new double [Z][X];
		double [][] Q = new double [Z][X];
		
		double t =0;
		double dt = 1/sf;
		for(int i = 0; i < Z; i++){
			sinw[i] = Math.sin(wc*t);
			cosw[i] = Math.cos(wc*t);
			t = t+dt;
		}
		
		//Demodulation	
		for(int z = 0; z < Z; z++){
			for(int x = 0; x < X; x++){
				I[z][x] = sinw[z]*rfdata[z][x];
				Q[z][x] = cosw[z]*rfdata[z][x];
				iqenv[z][x] = Math.sqrt(I[z][x]*I[z][x] + Q[z][x]*Q[z][x]);
			}
			z = z+1;
		}
		//apply low pass filter
		/*double[] coefs = lpfircoefs (sf,cf);
		//iqenv = filtrf(iqenv,coefs);
		I = filtrf(I,coefs);
		Q = filtrf(Q,coefs);
		for(int z = 0; z < Z; z++){
			for(int x = 0; x < X; x++){
			iqenv[z][x] = Math.sqrt(I[z][x]*I[z][x] + Q[z][x]*Q[z][x]);
			}
		}
		*/
		return iqenv;
	}
	
	// 2nd order low pass filter coefficients
	public static double[] lpfircoefs (double samplerate, double cutoff){
		
		double wc  = 2.0 * Math.PI * (cutoff / samplerate);
		double K = Math.tan(wc/2.0);
		double W = K*K;
		double Q = 100.0;
		double DE = 1.0+(1.0/Q)*K+W;
		double[] lpfircoefs = new double [6];
		
		lpfircoefs[3] = 1.0;
		lpfircoefs[4] = 2.0*(W-1.0)/DE;
		lpfircoefs[5] = (1.0-(1.0/Q)*K+W)/DE;
		lpfircoefs[0] = W/DE;
		lpfircoefs[1] = 2.0*(W/DE);
		lpfircoefs[2] = W/DE;
		
		/*
		double QcRaw  = (2 * Math.PI * cutoff) / samplerate; // Find cutoff frequency in [0..PI]
		double QcWarp = Math.tan(QcRaw); // Warp cutoff frequency
		double[] lpfircoefs = new double [6];
		
		double gain = 1 / (1+Math.sqrt(2)/QcWarp + 2/(QcWarp*QcWarp));
		lpfircoefs[3] = (1 - Math.sqrt(2)/QcWarp + 2/(QcWarp*QcWarp)) * gain;
		lpfircoefs[4] = (2 - 2 * 2/(QcWarp*QcWarp)) * gain;
		lpfircoefs[5] = 1;
		lpfircoefs[0] = 1 * gain;
		lpfircoefs[1] = 2 * gain;
		lpfircoefs[2] = 1 * gain;
		*/
		return lpfircoefs;
	}
	
	// apply filter to axial dimension over all channels
	public static double[][] filtrf(double [][] rfdata, double[] coefs){
		int X = rfdata[1].length;
		int Z = rfdata.length;
		double [] xv = new double [3];
		double [] yv = new double [3];
		double [][] filtrf = new double [Z][X];
		for (int x = 0; x < X; x++){
			for (int z = 0 ; z < 3 ;z++){
				xv[z] = 0;
				yv[z] = 0;
			}
			
			for (int z = 0 ; z < Z ;z++){
			xv[2] = xv[1]; 
			xv[1] = xv[0];
			xv[0] = rfdata[z][x];
			
			yv[2] = yv[1];
			yv[1] = yv[0];

			yv[0] = (coefs[0]*xv[0] + coefs[1]*xv[1] + coefs[2]*xv[2]
					-coefs[4]*yv[0] - coefs[5]*yv[1]);

			filtrf[z][x] = yv[0];
			}
		}
		return filtrf;
	}

	
	public static void main( String args[] ){
	
	}
}