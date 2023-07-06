function [A,B] = abmatrixsym(kmax,doprint)
% matrices A and B with cofficients for adding or removing dispersion
% in symbolic form as fractions, up to kmax rows, default 10
% if doprint, print result, default 1
% _____________________________________________________________________________
if nargin<2, doprint = 1; end
if nargin<1, kmax = 10; end
kmax = max(2,kmax);

% function x_{2j}=f^{(2j)}(0) with j2=2*j for f(x)=sin(x)/x
syms funx2ja(j2); funx2ja(j2) = (-1)^(j2/2) / (j2+1);
syms funcfactora(k); % scales first column of matrix to 1
funcfactora(k) = (-1)^k * (2*k+1);

% function x_{2j}=f^{(2j)}(0) with j2=2*j for f(x)=arcsin(x)/x 
syms funx2jb(j2); 
funx2jb(j2) = ( factorial(j2) / factorial(j2/2) )^2 / ( 2^j2*(j2+1) );
syms funcfactorb(k); % scales first column of matrix to 1
funcfactorb(k) = (2*k+1) * ( 2^k * factorial(k) / factorial(2*k) )^2;  

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

if nargin<4, doprint=1; end

syms AB;
for k=1:kmax
  fk = funcfactor(k);
  for ell=1:k
    AB(k,ell) = fk * Btilde(k,ell,funx2j);
    if doprint, fprintf(1,' %s',AB(k,ell)); end
  end
  AB(k,k+1:kmax) = 0;
  if doprint, fprintf(1,'\n'); end
end

end % createAB
% _____________________________________________________________________________
function a = Btilde(k,ell,funx2j)
% a = Btilde(k,ell,x2j)
% k and ell are intgers >= 0
% funx2j(j2) should be a function handle that returns x_{2 j}, j2=2*j

syms a; a = 0; 
if k==0 
  if ell==0, a=1; end % 1 for k=\ell=0 
  return;
end
% 0 outside range 1<=ell<=k, if not k=\ell=0 
if ell<=0 || ell>k || k<0, return; end

% this in not the most efficient approach:
for j=1:1+k-ell
  a = a + nchoosek(2*k-1,2*j-1) * funx2j(2*j) * Btilde(k-j,ell-1,funx2j);
end

end % Btilde
% _____________________________________________________________________________
%EOF
