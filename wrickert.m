function w = wrickert(fpeak,t)
% w = wrickert(fpeak,t)
% Ricker wavelet as a function of time t for peak frequency fpeak in Hz
% output: amplitudes ww(t)
%
% Ricker wavelet:
%  A*(1-(t/s)^2)*exp(-1/2*(t/s)^2), with s = 1/(pi*sqrt(2)*fpeak);
%  A=1 instead of 1/(s^3*sqrt(2*Pi))
%
s = 1/(pi*sqrt(2)*fpeak);
w = (t/s).^2; w = (1.-w).*exp(-0.5*w);
%EOF

