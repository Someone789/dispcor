function [ ftdt ] = FTDT( f, varargin )
%FTDT Forward Time Dispersion Transform
%   This function artificially ADDS time dispersion to a vector of samples
%   in time. It models the phase shift error. Because of the correspondence
%   of 1 sample to 1 dt, and maximum of 1/(2dt) Nyquist frequency also
%   inherent in the FFT, we actually lose all dt dependency.
%   f        [VEC]        measurements
%   varargin [INT] (opt.) subsampling rate
%  USE:
%      FTDT( data    ); % Applies dispersion to entire trace, assuming 1 sample equals 1 simulation timestep.
%      FTDT( data, 4 ); % Applies dispersion to entire trace, assuming 1 sample equals 4 simulation timesteps.
%
% 10-2016 EK: Basic version (Wang & Xu, 2015)
% 1 -2017 EK: Different function fn, true inverse filter operation, various optimisations.
% 02-2023 WAM: bug(?) fix in old code (jmax+2) instead of (jmax+1)
%              add nufft()

% Amount of subsampling
if length(varargin)==1
  ss = varargin{1};
else
  ss = 1;
end

% Zero-padd to twice the size
nf0 = length(f); nf = 2*nf0; f = [f(:); zeros(nf0,1)];
% The phase shift function
fn = @(omega) 2*sin( omega/2 );
jmax = nf0;

% Take altered Fourier Transform
%-% IFTDTuf= exp( -1i * fn([0:nf/2]*2*pi/nf/ss)'*ss * [0:nf-1] ) * f(:);
%-% IFTDTuf((nf/2+1):nf) = 0; % WAM: should be nf/2+2?

IFTDTuf = zeros(nf,1); dfn = 2*pi/nf/ss; 
if exist('nufft')
  IFTDTuf(1:jmax+1) = nufft(f,ss*(0:nf-1),fn((0:jmax)*dfn)/(2*pi) );
else
  s1 = (-1i*ss) * (0:nf-1) ;
  if nf > 5e4 % to avoid memory problems (WAM)
    for j=0:jmax, IFTDTuf(j+1) = exp( fn(j*dfn) * s1 ) * f(:); end
  else % rephrased from old code by Koene, with bug fix
    IFTDTuf(1:jmax+1) = exp( fn( (0:jmax)*dfn )' * s1) * f(:);
  end
end

% Take back to the time domain
ftdt   = ifft( IFTDTuf, 'symmetric'); % Bring back to the time domain 
ftdt   = ftdt(1:end/2);               % Remove zero-padded entries 

end
