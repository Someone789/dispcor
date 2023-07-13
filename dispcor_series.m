function [b,npmax] = dispcor_series(a,t,nsub,cdir,halforder,reduce_order,nhextra)
% [b,npmax] = dispcor_series(a,t,nsub,cdir,halforder,reduce_order,nhextra)
% cdir = 'f' or 'i' (add or remove dispersion)
% nsub for subsampling (integer>=1), modelling timestep dt=(t(2)-t(1))/nsub
% nhextra: if nhextra > 0, use smoothing derivatives by enlarging stencil
%   with nhextra points on each side
%
% b      is the corrected version of a
% npmax  is extra number of points, used for central differencing with
%        a stencil of 2*npmax+1 points
persistent Af Ai;

if nargin< 3 || isempty(nsub), nsub = 1; end % subsampling (integer >= 1)
if nargin< 4, cdir = 'i'; end
if nargin< 5, halforder = 2; else, halforder = max(1,halforder); end
if nargin< 6, reduce_order = 1; end
if nargin< 7 || isempty(nhextra), nhextra = 0; end

if isempty(cdir) || ~isa(cdir(1),'char') || ~contains('fair',cdir)
  error('cdir should be ''f'' (or ''a'') for adding or ''i'' (or ''r'') for removing dispersion');
end

if size(a,2)==1, a = transpose(a); atrans = 1; else, atrans = 0; end

if isempty(Af) || halforder > size(Af,1)
  % [Af,Ai] = local_loadaa_10(); % backup if halforder<=10
  kmax = max(halforder,10);
  [Af,Ai] = abmatrix(kmax,0); % up to dt^(2*kmax), kmax rows
  for k=1:kmax
    rk = 1/4^k / factorial(2*k+1); 
    Af(k,:) = rk*Af(k,:);
    sk = 1/(-16)^k * factorial(2*k)/( (2*k+1)*factorial(k)^2 );
    Ai(k,:) = sk*Ai(k,:);
  end
end % r*Af and s*Ai

if cdir(1) == 'f' || cdir(1) == 'a'
  A = Af;
else
  A = Ai;
end

% shift to t(1)=0
t = t-t(1); t(1)= 0; dt2 = t(2); dt = dt2/max(1,nsub); ntstart = 0;

% npextra is half of stencil width, to allow for central differencing
if reduce_order, npextra =   halforder   + floor((1+halforder)/2);
else,            npextra = 2*halforder-1 + floor((1+halforder)/2);
end
npextra = npextra+nhextra;

b = zeros(size(a)); npmax = 0; jorder = 2*halforder;
for k=halforder:-1:1
  % reverse loop assuming values are smaller for larger k
  if reduce_order, jorder = 2*(halforder-(k-1)); end
  for ell=1:k
    ider = 2*k+ell; % -sign in (-dt)^ell
    sfac = A(k,ell)/(-dt)^ell; % (1/dt^(2*k)) already in rf(k)
    % zero padding with npextra points on each side
    if nsub>1
      sfac = sfac/nsub^ider;
    end
    % zero padding with npextra on both ends to enable for central differencing
    a2 = [zeros(1,npextra) (a .* t.^ell) zeros(1,npextra)];
    [da,nda] = apply_fd_op2(a2,ider,jorder,nhextra,sfac);
    npmax = max(npmax,nda); % internal check for stencil width (<=npextra?)
    da = da(npextra+1:end-npextra); % remove padding
    b = b+da;
  end
end

if ntstart>0, b = b(ntstart+1:end); end
if npmax>npextra
  error('nd=%d > npextra=%d; increase initial npextra',npmax,npextra);
end
if atrans, b = transpose(b); end
% _____________________________________________________________________________
function [Af,Ai] = local_loadaa_10()
% [Af,Ai] = local_loadaa_10() backup for abmatrix()
% up to dth^20
Af=[1, 0,0,0,0,0,0,0,0,0;
    1, 5/3, 0,0,0,0,0,0,0,0;
    1, 7, 35/9, 0,0,0,0,0,0,0;
    1, 123/5, 42, 35/3, 0,0,0,0,0,0;
    1, 253/3, 341, 770/3, 385/9, 0,0,0,0,0;
    1, 2041/7, 38324/15, 11869/3, 5005/3, 5005/27, 0,0,0,0;
    1, 1023, 18759, 484484/9, 130130/3, 35035/3, 25025/27, 0,0,0;
    1, 32759/9, 413882/3, 31384873/45, 25967942/27, 4254250/9, 2382380/27, 425425/81, 0,0;
    1, 65531/5, 107585809/105, 133275614/15, 297217817/15, 142908766/9, 47205158/9, 6466460/9, 8083075/243, 0;
    1, 524277/11, 7704576, 112990891, 9860196129/25, 485601753, 6835048220/27, 539949410/9, 56581525/9, 56581525/243];

Ai=[1, 0,0,0,0,0,0,0,0,0;
    1, 5/27, 0,0,0,0,0,0,0,0;
    1, 7/25, 7/405, 0,0,0,0,0,0,0;1, 2067/6125, 6/175, 1/945, 0,0,0,0,0,0;
    1, 4477/11907, 4829/99225, 22/8505, 11/229635, 0,0,0,0,0;
    1, 150761/373527, 990964/16372125, 1261/297675, 13/93555, 13/7577955, 0,0,0,0;
    1, 78103/184041, 48787/693693, 106684/18243225, 314/1216215, 1/173745, 1/19702683, 0,0,0;
    1, 164187887/372683025, 477888394/6087156075, 430651327/58530346875, 4836262/12314176875, 9826/820945125, 68/351833625, 17/13299311025, 0,0;
    1, 335021699/738720125, 309336853817/3621857864625, 1508301434/172469422125, 531060317/995015896875, 1390078/69780335625, 31046/69780335625, 76/13956067125, 19/678264862275, 0;
    1, 10892077437/23467660931, 169784/1859715, 3936101287/393230282445, 118566891343/175549233234375, 15738043/540151486875, 7036/8741712375, 3646/265165275375, 1/7576150725, 1/1841004626175];
% _____________________________________________________________________________
%EOF



