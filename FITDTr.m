function g = FITDTr(fiflag,f,nsubsampling)
% g = FITDTr(fiflag,f,[nsubsampling]) 
%   forward or inverse time dispersion transform
%   variant with slow Fourier transform to get back to time
% fiflag = 0 or 'f' or 'a' adds    time dispersion
% fifalg = 1 or 'i' or 'r' removes time dispersion (default if empty)
% f and g are column vectors
% nsubsampling (integer>=1, default 1) is the amount of subsampling

if nargin<3, ss = 1; 
else, ss = max(1,nsubsampling); 
end

doremove = true;
if ~isempty(fiflag)
  if ischar(fiflag(1)), doremove = ~( fiflag(1)=='f' || fiflag(1)=='a' ); % add
  else, doremove = (fiflag~=0);
  end
end

if size(f,1)==1, f = transpose(f); end % make column vector
% zero-padd to twice the size
nf0 = length(f); f = [f(:); zeros(nf0,1)];

% df = 1/(nf*dt) after padding
% omega*dt = 2*pi*(0:nf0)/(nf*dt)*dt=(0:nf0)*2*pi/nf=omdt
% t*sin(omdt/2)/(dt/2) = t*fn(omdt)/dt with t=(0:nf-1)*dt

% g has fmin=0 at index 1, fmax=nf0*dt at index nf0+1, 
% 0==max(abs(g(2:nf0)-conj(flipud(g(end-nf0+2:end))))) for even nf
g = fft(f);

nf = 2*nf0;
if doremove
  jmax = nf0;
  omdt2 = (0:jmax)*(pi/nf/ss); % omega dt/2
  fdt2 = sin(omdt2)*(ss/pi);   % phase shift function times 1/(2*pi)
else % add
  jmax = floor(nf/pi);         % drop too high frequencies with omdt2>1
  omdt2 = (0:jmax)*(pi/nf/ss); % omega dt/2
  fdt2 = asin(omdt2)*(ss/pi);  % phase shift function times 1/(2*pi)
end
g = local_iftslow(g(1:jmax+1),fdt2,(0:nf0-1));

end
% _____________________________________________________________________________
function b = local_iftslow(a,freqa,timeb)
% b = iftslow(a,freqa,timeb) 
% Fourier transform from complex-valued a as a function of
% frequencies freqa>=0 to times timeb and real-valued b

if max(abs(imag(freqa)))>0
  warning('iftslow cannot handle complex frequencies; taking real part');
  freqa = real(freqa);
end

dfr = [0 diff(freqa) 0];
% 0 at end produces [1/2,1,...,1,1/2] for the trapezoid rule
% intervals and integration weights:
%   dfr = 0.5*(dfr(1:end-1)+dfr(2:end)); cc = 2*dfr;
%   factor 2 in cc=2*dfr for conjugate symmetry at dropped negative frequencies

% ca = cc'.*a; for computational purposes
ca = (dfr(1:end-1)+dfr(2:end))' .* a;

if exist('nufft')
  % matlab fft has exp( -2*pi*1i*(j-1)*(k-1)/n ) and nufft also has - sign
  % -freqa to get opposite of sign of fft and nufft
  b = real( nufft(ca,-freqa,timeb) ); 
else
  tpi2 = 2*pi*timeb;
  b = zeros(length(timeb),1); % column vector
  for j=1:length(b)
    % freqa is row vector, ca column vector, result should be scalar
    b(j) = real( exp(1i*tpi2(j)*freqa) * ca );
  end
end

end
% _____________________________________________________________________________
%EOF
