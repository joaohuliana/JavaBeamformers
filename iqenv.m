function [iq_data,I,Q] = iqenv(data,prm)
%%%%This script performs envelope detection by quadrature demodulation
%Autor: J. H. Uliana
%10-2023
%iqenv(data, prm)
%data: beamformed RF data (RF image)
%prm: 
    %sampling frequency (prm.sf) [Hz]
    %transducer's central frequency (prm.cf) [Hz]
    
sf = prm.sf;
cf = prm.cf;
wc = cf*2*pi;

%dimensions
Z = size(data,1);
X = size(data,2);

%IQBuffer
I = zeros(Z,X);
Q = zeros(Z,X);

t = 0:1/sf:(Z-1)*(1/sf);
cosw = cos(wc*t);
sinw = sin(wc*t);

[b,a] = butter(5,0.5,'low');
for x = 1:X
     I(:,x) = filtfilt(b,a,sinw'.*data(:,x));
     Q(:,x) = filtfilt(b,a,cosw'.*data(:,x));
end
iq_data = sqrt(I.^2+Q.^2);
    