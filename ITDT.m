function [ itdt ] = ITDT( f, varargin )
%ITDT Inverse Time Dispersion Transform
%   This function artificially REMOVES time dispersion from a vector of samples
%   in time. It models the phase shift error. It uses the fact that every
%   sample equals 1 time step (assuming equidistant sampling), removing the dt dependency.
%   f        [VEC]        measurements
%   varargin [INT] (opt.) subsampling rate
%  USE: 
%      ITDT( data    ); % Removes dispersion from entire trace, assuming 1 sample equals 1 simulation timestep.
%      ITDT( data, 4 ); % Removes dispersion from entire trace, assuming 1 sample equals 4 simulation timesteps.
%
% 10-2016 EK: Basic version (Wang & Xu, 2015)
% 1 -2017 EK: Various optimisations.
% 11-2022 WAM: work-around for memory problem
% 02-2023 WAM: bug (?) fix in old code (jmax+2) instead of (jmax+1)
%              add nufft()

% Amount of subsampling
if length(varargin)==1
  ss = varargin{1};
else
  ss = 1;
end

% Zero-padd to twice the size
f = [        f(:)       ;
  zeros(length(f), 1 )];
nf = length(f);
% The phase shift function
fn = @(omega) 2*asin( omega/2 );
jmax = floor(nf/pi);

% Take altered Fourier Transform
%-% IFTDTuf= exp( -1i * fn([0:floor(nf/pi)]*2*pi/nf/ss)'*ss * [0:nf-1]) * f(:); 
%-% IFTDTuf(ceil(nf/pi):nf)= 0; % Block out omega*dt>2 ==> f>1/pi
% WAM should be IFTDTuf(1+ceil(nf/pi):nf)= 0;  ?

IFTDTuf = zeros(nf,1); dfn = 2*pi/nf/ss; 
if exist('nufft') % less accurate?  
  IFTDTuf(1:jmax+1) = nufft(f,ss*(0:nf-1),fn( (0:jmax)*dfn )/(2*pi) );
else
  s1 = (-1i*ss) * (0:nf-1) ;
  if nf > 5e4 % to avoid memory problems (WAM)
    for j=0:jmax, IFTDTuf(j+1) = exp( fn(j*dfn) * s1 ) * f(:); end
  else % rephrased from old code by Koene, with bug fix
    IFTDTuf(1:jmax+1) = exp( fn( (0:jmax)*dfn )' * s1) * f(:);
  end
end

itdt   = ifft( IFTDTuf, 'symmetric'); % Bring back to the time domain 
itdt   = itdt(1:end/2);               % Remove zero-padded entries 

end
%EOF
