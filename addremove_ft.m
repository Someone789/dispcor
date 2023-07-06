% addremove_ft.m
% orginal approach from Koene et al. and alternatives
% _____________________________________________________________________________
close all; clear all; % overdone
froot = mfilename();
% _____________________________________________________________________________
% Ricker, form Koene et al. fig 3 add/rm 8 Hz  (0:0.015:0.4)
tmin = 0; tmax = 0.4; dt = 0.015; doplot = 1;
type_wavelet = 0+1; % or 1

switch type_wavelet
  case 1, Tw = 0.18; ww = @(t0) wavelet16(t0-Tw,Tw,0)';
  otherwise, fpeak = 8; ww = @(t0) wrickert(fpeak,t0-tmax/2)';
end

% choose nsub for subsampling, integer >=1; 1 is no subsampling
nsub = 1;
if 1 % or if 0 
  nsub = 2; tmax = 2*tmax;
end
if nsub>2, doplot = 0; end
% _____________________________________________________________________________
t0 = (tmin:dt/nsub:tmax); wt0 = ww(t0);

% Koene et al., recoded
wt1 = FTDT(wt0); % add
wt2 = ITDT(wt1); % remove
if nsub>1 % subsampled
  t0sub = (tmin:dt:tmax); wt0sub = wt0(1:nsub:end);
  wt1sub = FTDT(wt0sub,nsub); % add
  wt2sub = ITDT(wt1sub,nsub); % remove
end

if doplot
  close all; clf; hfig1 = figure(1); r = groot; pos = get(hfig1,'Position');  
  pos(4) = min(2*pos(3),round(0.8*r.ScreenSize(4))); 
  pos(2) = r.ScreenSize(4)-round(1.1*pos(4));
  set(hfig1,'Position',pos);
 
  subplot(3,2,1); plot(t0,wt0,'ro-',t0,wt2,'k+-');
  legn = {'input','add&remove'};
  if nsub>1
    hold on; plot(t0sub,wt2sub,'b.-'); hold off; legn = [legn 'subsampled'];
  end
  xlabel('Time (s)'); legend(legn); axis tight;
  subplot(3,2,2); plot(t0,wt2-wt0,'k.'); legn = {'difference'};
  if nsub>1
    hold on; plot(t0sub,wt2sub-wt0sub,'bo'); hold off; legn = [legn 'subsampled'];
  end
  xlabel('Time (s)'); legend(legn); axis tight; pause(0.1);
end
% nt0=length(wt0); nt=2*nt0; bt=wvlt2beta(t0,wt0,nt0,true);
% if mod(nt0,2)==0, assert(length(bt.freq)== nt0+1); end % if nt0 even
% assert(2*bt.freq(end)*dt == 1); %  (Nyquist)
wr1 = FITDTr('add',wt0); % add dispersion with alternative approach
if nsub>1
  wr1sub = FITDTr('add',wt0sub,nsub);
end

if doplot
  subplot(3,2,3); plot(t0,wt1,'ro',t0,wr1,'k+-'); legn={'input','add'};
  if nsub>1
    hold on; plot(t0sub,wt1sub,'r*',t0sub,wr1sub,'b+'); hold off;
    legn = [legn 'subsampled1' 'subsampled2'];
  end
  xlabel('Time (s)'); legend(legn); axis tight;
  
  subplot(3,2,4); plot(t0,wr1-wt1,'k.');
  if nsub>1, hold on; plot(t0sub,wr1sub-wt1sub,'bo'); hold off; end
  xlabel('Time (s)'); legend('difference'); axis tight; pause(0.1);
end

wr2 = FITDTr('rem',wr1); % remove dispersion with alternative approach
if nsub>1
  ws2sub = FITDTr('rem',wr1sub,nsub);
end

if doplot
  subplot(3,2,5); plot(t0,wt0,'r-',t0,wr2,'k--'); legn = {'input','add&remove'};
  if nsub>1
    hold on; plot(t0sub,wt2sub,'b+',t0sub,ws2sub,'go'); hold off;
    legn = [legn {'subsampled 1','subsampled 2'}];
  end
  xlabel('Time (s)'); legend(legn); axis tight;
  if nsub>1
    subplot(3,2,6); plot(t0sub,wt2sub-wt0sub,'ro',t0sub,ws2sub-wt0sub,'g+');
    xlabel('Time (s)'); legend(legn(end-1:end)); axis tight;
  end
  pause(0.1);
end

%EOF
