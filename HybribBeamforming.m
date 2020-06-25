% This promgration is part of a study in the 5G network article, supported by
% Research Tecchnology privated Company.
% //Large Finite Array For Hybrid Beamforming//*

%% Array parameters.
N = 11;
fc = 28e9;
az = 30;
el = 20;

% Find phase shifts for azimuth control.
l = design(linearArray,fc);
elem = l.Element;
elem.Tilt = 90;
l.NumElements = N;
figure
show(l)
ps_az = phaseShift(l,fc,[az;0]);

% Find phase shifts for elevation control.
elem.Tilt = 90;
elem.TiltAxis = [0 1 0];
l.Tilt = 90;
l.TiltAxis = [0 1 0];
l.ElementSpacing = 1.05*(elem.Length) ;
figure
show(l)
ps_el = phaseShift(l,fc,[0;el]);

% Create subarrays (nxn).
l.Tilt = 0;
elem.Tilt = 0;
l.PhaseShift = ps_az;
c = conformalArray;
zposn = fliplr((-N+1)/2:1:(N-1)/2);
for i = 1:N
   c.Element{i} = l;
   c.ElementPosition(i,:) = [0,0,zposn(i)*l.ElementSpacing];
end
figure
show(c)

% Assign phaseshifts and plot pattern.
c.PhaseShift = ps_el;
figure
pattern(c,fc);
figure
patternElevation(c,fc,az);

% Array with Large Reflector Backing.
lambda = physconst('lightspeed')/fc;
ref_offset = lambda/4;
p = platform;
p.FileName = 'Groundplane.stl';
p.Units = 'm';
p.Tilt = 90;
f = installedAntenna;
f.Platform = p;
f.Element = c.Element;
f.ElementPosition = c.ElementPosition;
f.ElementPosition(:,2) = ref_offset;
f.FeedPhase = ps_el;
figure
show(f)

% Approximate Array Pattern.
figure
pattern(f,fc)


