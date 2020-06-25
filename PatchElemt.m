% This promgration is part of a study in the 5G network article, supported by
% Research Tecchnology privated Company.


%% Parameter desing of the frecuency band and light.
fc = 77e9;
fmin = 73e9;
fmax = 80e9;
vp = physconst('lightspeed');
lambda = vp/fc;

% Array (2x4) FMCW Antenna desing.
cosineElement = phased.CosineAntennaElement;
cosineElement.FrequencyRange = [fmin fmax];
cosinePattern = figure;
pattern(cosineElement,fc)

% Array Pattenr and the bumper rectangular Array (2x4).
Nrow = 2;
Ncol = 4;
fmcwCosineArray = phased.URA;
fmcwCosineArray.Element = cosineElement;
fmcwCosineArray.Size = [Nrow Ncol];
fmcwCosineArray.ElementSpacing = [0.5*lambda 0.5*lambda];
cosineArrayPattern = figure;
pattern(fmcwCosineArray,fc);

% Design realistic Antenna Patch this half-wavelet is 60GHz and 66 Ghz.
patchElement = design(patchMicrostrip, fc);
patchElement.Tilt = 90;
patchElement.TiltAxis = [0 1 0];
figure
show(patchElement)
axis tight
view(140,20)

% Directivity Pattenr in 3D for the 60 GGHz and 66GHz Antenna peak
% directivity  6-9 dBi
pattern(patchElement,fc)

% Resonances patch (azimut=0°) with respect the impendace behavior.
Numfreqs = 21;
freqsweep = unique([linspace(fmin,fmax,Numfreqs) fc]);
impedance(patchElement,freqsweep);

% Stabilish Bandwith recpecto to the reflection coefficient S11=-10dBi.
s = sparameters(patchElement,freqsweep);
figure
rfplot(s,'m-.')
hold on
line(freqsweep,ones(1,numel(freqsweep))*-10,'LineWidth',1.5)
hold off

% Pattern confirmation at center to corner frecuencies, patterns plots
% 76.5-77.6 GHz.
fc2 = 77.6e9;
lambda_fc2 = vp/77.6e9;
fmcwPatchArray = phased.URA;
fmcwPatchArray.Element = patchElement;
fmcwPatchArray.Size = [Nrow Ncol];
fmcwPatchArray.ElementSpacing = [0.5*lambda_fc2 0.5*lambda_fc2];
az = -180:5:180;
el = -90:5:90;
patchArrayPattern = figure;
pattern(fmcwPatchArray,fc,az,el);

% Array from isolated radiators of the plots pattern by a frecuency(lambda/2).
[Dcosine_az_zero,~,eln] = pattern(fmcwCosineArray,fc,0,el);
[Dcosine_el_zero,azn] =  pattern(fmcwCosineArray,fc,az,0);
[Dpatch_az_zero,~,elp] = pattern(fmcwPatchArray,fc,0,el);
[Dpatch_el_zero,azp] =  pattern(fmcwPatchArray,fc,az,0);
elPattern = figure;
plot(eln,Dcosine_az_zero,eln,Dpatch_az_zero,'LineWidth',1.5)
axis([min(eln) max(eln) -40 17])
grid on
xlabel('Elevation (deg.)')
ylabel('Directivity (dBi)')
title('Array Directivity Variation-Azimuth = 0 deg.')
legend('Cosine element','Patch Antenna','Location','best')
azPattern = figure;
plot(azn,Dcosine_el_zero,azn,Dpatch_el_zero,'LineWidth',1.5)
axis([min(azn) max(azn) -40 17])
grid on
xlabel('Azimuth (deg.)')
ylabel('Directivity (dBi)')
title('Array Directivity Variation-Elevation  = 0 deg.')
legend('Cosine element','Patch Antenna','Location','best')

