% addremove_series.m
% illustrates series approach to dispersion correction
% - first time run may take long because a list of finite-difference
%   operators is computed
%
% TODO: make dispcor_series.m 
%         consistent with max derivative and stencil width and nhextra 
%         and kmax=halfordermax ...
%       reuse data from file?
%       cleanup everywhere
% _____________________________________________________________________________
close all; clear all; % overdone
% _____________________________________________________________________________
% after Koene et al., fig 3 with add&rm, 8 Hz Ricker, time (0:0.015:0.4)
tmin = 0; tmax = 0.4; doplot = 1;
type_wavelet = 0; % 0 (Ricker) or 1 (bump), details below
% for series
halfordermax = 4; % up to dt^(2*halfordermax)
nhextra = 3; % extra points for stencil to suppress noise (2*nhextra added)
reduce_order = 1; % reduce derivative order when dt^(2k) increases (1 is OK)
% choose nsub for subsampling, integer >=1; 1 is no subsampling
nsub = 1;
if 0 % 1 or 0 
  nsub = 2;
end

switch type_wavelet
  case 1 % bump
    Tw = 0.18; ww = @(t0) wavelet16(t0-tmax/2,Tw)';
    dt = 0.015/4; % timestep has to be small enough for convergence
  otherwise % Ricker
    fpeak = 8; ww = @(t0) wrickert(fpeak,t0-tmax/2)';
    dt = 0.015/4; % timestep has to be small enough for convergence
end

if nsub>2, doplot = 0; end
% _____________________________________________________________________________
t0 = (tmin:dt/nsub:tmax); wt0 = ww(t0);
if nsub>1
  wt0 = wt0(1:nsub:end); t0 = t0(1:nsub:end);
end

for halforder=1:halfordermax
  % add dispersion
  ac2 = dispcor_series(wt0,t0-t0(1),nsub,'f',halforder,reduce_order,nhextra);
  ws1 = wt0+ac2;
  % remove dispersion  
  bc2 = dispcor_series(ws1,t0-t0(1),nsub,'i',halforder,reduce_order,nhextra);
  ws2 = ws1+bc2;
  
  fprintf(1,'    nsub=%d; order/2=%d; nextra/2=%d: error %g\n',...
    nsub,halforder,nhextra,max(abs(ws2-wt0)));

  subplot(2,1,1); plot(t0,wt0,'k-',t0,ws1,'r--'); axis tight;
  legend('original','added','Location','Best');

  subplot(2,1,2); plot(t0,ws2-wt0,'k-');axis tight;
  title(sprintf('error (nsub=%d, order=2*%d, 2*%d extra points)',...
    nsub,halforder,nhextra));
  pause(1);  
end
% _____________________________________________________________________________
%EOF
