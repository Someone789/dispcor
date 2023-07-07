function [w,nph1,wsym,wnum,wden0,iorder2] = get_fd_op2(ider,iorder,nhextra)
% [w,nw,wsym,wnum,wden,iorder2] = get_fd_op2(ider,iorder,nhextra)
%   where wsym=wnum/wden
%  w is an row vector with nw weights, symmetric or anti-symmetric
% 
%  w=wnum/wden0 has width nw=2*nph1+1, wsym is it symbolic form
%  FD operator of order iorder for ider-th derivate
%  factor 1/spacing^ider not included
%  2*nhextra is the number of extra points for smoothing
%  iorder2>=iorder is the order actually chosen
%
%  https://doi.org/10.1080/00207160.2012.666348
%  Timothy J. McDevitt (2012) Discrete Lanczos derivatives of noisy data, 
%    International Journal of Computer Mathematics, 89:7, 916-931,
%    DOI: 10.1080/00207160.2012.666348

%   r=ider, p=iorder, N=1+2n stencil width, with r+p-1<=N
%   
%   examples: get_fd_op2(1,2,1) in eq.(32) [-2 -1 0 1 2]/10
%             get_fd_op2(2,2,3) in eq.(36) [28 7 -8 -17 -20 -17 -8 7 28]/462
%             f=get_fd_op2(2,2,25); dot(f,f)
%             f=get_fd_op2(4,4,4); fprintf(1,'%d %g\n',length(f),f'*f)
%                  WRONG not 35/46189=7.5776e-04
% f=get_fd_op2(4,4,0) in eq.(39)

if nargin<3, nhextra = 0; end
if ider<0, error('ider<0 not implemented (integration)'); end
if ider==0, w = 1; nw = 0; return; end
iorder = max(1,iorder); % iorder should be >0, force if not

% stencil width
% exact interpolation for polynomial of degree p-1: p points
%   rth derivative: p+r points, 2n+1>=p+r.
%   if 2n+1>p+r, than order p can be increased by 1

nwidth = 2*floor( (ider+1)/2 )+(iorder-1);
nph0 = ceil((nwidth-1)/2); np0 = 2*nph0+1;
iorder2 = np0-ider+1;
nph1 = nph0+nhextra;
nw = 1+2*nph1; % full stencil width

% moments
syms jj wsym; A = sym(zeros(np0,nw)); b = sym(zeros(np0,1));
jj = sym((-nph1:nph1)); A(1,:) = 1; for k=1:np0-1, A(k+1,:) = jj.^sym(k); end
b(ider+1) = factorial(sym(ider));
% see paper: V'=A, wsym = A' inv(A A') b
wsym = A'*linsolve(A*A',b);
% express wsym as wnum/wden0
wden0 = sym(1); k = 0;
while k<2*length(wsym)
  k = k+1; [wnum,den] = numden(wden0*wsym); den = unique(den);
  if length(den)==1 && den==1, break; end
  wden0 = wden0*max(den);
end
% check result
werr = unique(wsym-wnum/wden0);
if length(werr)>1 || werr~=0
  warning('wsym != wnum/wden0'); wnum=[]; wden0=[];
end

w = double(wsym)'; wnum  = wnum'; % transpose
%EOF
