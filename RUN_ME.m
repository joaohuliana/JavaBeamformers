clc
clear
close all
%%%%This script demonstrate the beamformers
%Autor: J. H. Uliana
%09-2025
%prm: 
    %sampling frequency (prm.sf) [Hz]
    %transducer's central frequency (prm.cf) [Hz]
    %speed of sound (prm.c) [m/s]
    %element pitch (prm.dx) [mm]
    %axial diferential (prm.dz) [mm]
    %Aperture (prm.aperture) [number of channels]
    %directional function angle limit (prm.ctag) [º]
    %quality factor (for SLSC) (prm.Q)
    %Photoacoustic or BMode reconstruction (prm.PA) [boolean]


%setting Java
javaaddpath(pwd);
O = us_tools;

%load .mat file
[file.name, file.path] = uigetfile('*.mat','Load raw data...');
load([file.path file.name]);

%dimensions
Z = size(rawData,1);
X = size(rawData,2);
dimx = -X/2*prm.dx:prm.dx:X/2*prm.dx;
dimz = 0:prm.dz:Z*prm.dz;

%delay map parameters
dprm(1) = prm.dx; %[mm]
dprm(2) = prm.c; %[m/s]
dprm(3) = prm.sf; %[Hz]

%delay map
if prm.PA
    delays = javaMethod('delaysPA',O,Z,X,dprm);
else
    delays = javaMethod('delays',O,Z,X,dprm);
end

%beamformer parameters
jprm(1) = prm.dz; %[mm]
jprm(2) = prm.dx; %[mm]
jprm(3) = prm.aperture; %[channels]
jprm(4) = prm.ctag; %[°]

slscprm = jprm;
slsc(4) = prm.Q;

%beamformers
bf1 = javaMethod('das',O,rawData,delays,jprm,1);
disp('DaS complete')
bf2 = javaMethod('fdmas',O,rawData,delays,jprm,1);
disp('fDMaS complete')
bf3 = javaMethod('slsc',O,rawData,delays,slscprm,1);
disp('SLSC complete')
bf4 = javaMethod('dmas',O,rawData,delays,jprm,1);
disp('DMaS complete')

%envelope detection
bf1 = abs(iqenv(bf1,prm));
bf2 = abs(iqenv(bf2,prm));
bf3 = abs(bf3);
bf4 = abs(iqenv(bf4,prm));

%normalize
bf1 = bf1/max(bf1(:));
bf2 = bf2/max(bf2(:));
bf3 = bf3/max(bf3(:));
bf4 = bf4/max(bf4(:));


%% Show Images
close all
ax1 = subplot(2,2,1);
if prm.PA
    imagesc(dimx,dimz,bf1)
    colormap(ax1,'hot');
else
    imagesc(dimx,dimz,sqrt(bf1))
    colormap(ax1,'gray');
end
axis(ax1, 'image');
title(ax1,'Delay and Sum');
xlabel(ax1,'Lateral [mm]')
ylabel(ax1,'Axial [mm]')
c1 = colorbar(ax1);
ylabel(c1, 'Intensity [a.u.]')

ax2 = subplot(2,2,2);
if prm.PA
    imagesc(dimx,dimz,bf2)
    colormap(ax2,'hot');
else
    imagesc(dimx,dimz,sqrt(bf2))
    colormap(ax2,'gray');
end
axis(ax2, 'image');
title(ax2,'filtered Delay Multiply and Sum');
xlabel(ax2,'Lateral [mm]')
ylabel(ax2,'Axial [mm]')
c2 = colorbar(ax2);
ylabel(c2, 'Intensity [a.u.]')

ax3 = subplot(2,2,3);
imagesc(dimx,dimz,bf3)
if prm.PA
    colormap(ax3,'hot');
else
    colormap(ax3,'gray');
end
axis(ax3, 'image');
title(ax3,'Short-Lag Spatial Coherence');
xlabel(ax3,'Lateral [mm]')
ylabel(ax3,'Axial [mm]')
c3 = colorbar(ax3);
ylabel(c3, 'Intensity [a.u.]')

ax4 = subplot(2,2,4);
if prm.PA
    imagesc(dimx,dimz,bf4)
    colormap(ax4,'hot');
else
    imagesc(dimx,dimz,sqrt(bf4))
    colormap(ax4,'gray');
end
axis(ax4, 'image');
title(ax4,'Delay Multiply and Sum');
xlabel(ax4,'Lateral [mm]')
ylabel(ax4,'Axial [mm]')
c4 = colorbar(ax4);
ylabel(c4, 'Intensity [a.u.]')