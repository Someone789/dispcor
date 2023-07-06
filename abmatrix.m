function [A,B] = abmatrixsym(kmax,doprint)
% matrices A and B with cofficients for adding or removing dispersion
% up to kmax rows, default 10
% if doprint, print result, default 1
% _____________________________________________________________________________
if nargin<2, doprint = 1; end
if nargin<1, kmax = 10; end
kmax = max(2,kmax);

% function x_{2j}=f^{(2j)}(0) with j2=2*j for f(x)=sin(x)/x
funx2ja = @(j2) (-1)^(j2/2) / (j2+1);
funcfactora = @(k) (-1)^k * (2*k+1); % scales first column to 1

% function x_{2j}=f^{(2j)}(0) with j2=2*j for f(x)=arcsin(x)/x 
funx2jb = @(j2) ( factorial(j2) / factorial(j2/2) )^2 / ( 2^j2*(j2+1) );
funcfactorb = @(k) (2*k+1) * ( 2^k * factorial(k)/factorial(2*k) )^2;
  % scales first column to 1

if doprint, fprintf(1,'\nA=\n'); end
A = createAB(kmax,funx2ja,funcfactora,doprint);
if doprint, fprintf(1,'\nB=\n'); end
B = createAB(kmax,funx2jb,funcfactorb,doprint);

end
% _____________________________________________________________________________
function AB = createAB(kmax,funx2j,funcfactor,doprint)
% AB = createAB(kmax,funx2j,funcfactor,doprint)
% creates matrix of size kmax*kmax using function handle funx2j(j2), j2=2*j,
% that evaluates x_{2 j} = f^{(2 j)}(0) for even function f(x)
% funcfactor(k) is  scaling factor

if ~isa(funx2j,'function_handle')
  error('  second argument should be a function handle for x(j2)=x_{2*j}');
end
if ~isa(funcfactor,'function_handle')
  error('  third argument should be a function handle for scalefactor(k)');
end
if nargin<4, doprint=1; end

AB = zeros(kmax,kmax);
for k=1:kmax
  for ell=1:k
    AB(k,ell) = funcfactor(k)*Btilde(k,ell,funx2j);
    if doprint, fprintf(1,' %g',AB(k,ell)); end
  end
  if doprint, fprintf(1,'\n'); end
end

end % createAB
% _____________________________________________________________________________
function a = Btilde(k,ell,funx2j)
% a = Btilde(k,ell,x2j)
% k and ell are intgers >= 0
% funx2j(j2) should be a function handle that returns x_{2 j}, j2=2*j
a = 0; 
if k==0 
  if ell==0, a=1; end
  return;
end
if ell<=0 || ell>k || k<0, return; end
% this in not the most efficient approach
for j=1:1+k-ell
  a = a + nchoosek(2*k-1,2*j-1) * funx2j(2*j) * Btilde(k-j,ell-1,funx2j);
end

end % Btilde
% _____________________________________________________________________________
%EOF
