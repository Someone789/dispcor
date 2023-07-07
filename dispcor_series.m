function [b,ndmax] = dispcor_series(a,t,nsub,cdir,halforder,reduce_order,nhextra,npextra)
% [b,npmax] = dispcor_series(a,t,nsub,cdir,halforder,reduce_order,nhextra,npextra)
% cdir = 'f' or 'i' (add or remove dispersion)
% nsub for subsampling (integer>=1)
% nhextra: if nhextra > 0, use smoothing derivatives by enlarging stencil
%   with nhextra points on each side
% npextra: zero padding with npextra points on each side
% 
persistent Af Ai;

if nargin< 3 || isempty(nsub), nsub = 1; end % subsampling (integer >= 1)
if nargin< 4, cdir = 'f'; end
if nargin< 5, halforder = 1; end
if nargin< 6, reduce_order = 1; end
if nargin< 7 || isempty(nhextra), nhextra = 0; end
if nargin< 8 || isempty(npextra)
  % npextra<=half of stencil width
  npextra = 0;
  for k=halforder:-1:1
    if reduce_order, jorder = 2*max(1,halforder-(k-1)); end
    for ell=k:-1:1
      ider = 2*k+ell; % -sign in (-dt)^ell
      nwidth = 2*floor( (ider+1)/2 )+(jorder-1); % from get_fd_opt2.m
      nph0 = ceil((nwidth-1)/2); % np0 = 2*nph0+1; iorder2 = np0-ider+1;
      nph1 = nph0+nhextra; npextra = max(npextra,nph1);
      fprintf(1,'ider=%d, jorder=%d, nph0=%d\n',ider,jorder,nph0);
    end
  end
end
if isempty(cdir) || ~isa(cdir(1),'char') || ~contains('fair',cdir)
  error('cdir should be f or a for adding or i or r for removing dispersion');
end

if size(a,2)==1, a = transpose(a); atrans = 1; else, atrans = 0; end

dt2 = mean(diff(t)); dt = dt2/max(1,nsub); ntstart = 0;
if t(1)>0
  ntstart = round(t(1)/dt2); % rt = t(1)-jt*dt2;
  if ntstart>0 % add samples at start, to let t(1)=0
    t = [ (0:ntstart-1)*dt2 t(:)'];
    a = [zeros(1,ntstart)   a(:)'];
  end
end
if abs(t(1))>1.e-15*max(abs(t))  
  warning('t(1) = %g non-zero\n',t(1));  
else, t(1)=0;
end

if isempty(Af)
  % [Af,Ai] = local_loadaa_10();
  kmax = 10;
  [Af,Ai] = abmatrix(kmax,0); % up to dt^(2*kmax), kmax rows
end % Af and Ai

if cdir(1) == 'f' || cdir(1) == 'a'
  A = Af; cdir(1) = 'f';
else
  A = Ai;
end

if length(A)<halforder
  warning('order/2=%d not implemented; at most %d',halforder,length(A));
  halforder = length(A);
end

rf = zeros(1,halforder); % scale factor / dt^(2*k)
if cdir == 'f'
  for k=1:halforder
    rf(k) = 1/4^k / factorial(2*k+1);
  end
else
  for k=1:halforder
    rf(k) = 1/(-16)^k * factorial(2*k)/( (2*k+1)*factorial(k)^2 );
  end
end

b = zeros(size(a)); ndmax = 0;
jorder = 2*max(1,halforder);
for k=halforder:-1:1
  if reduce_order, jorder = 2*max(1,halforder-(k-1)); end
  for ell=k:-1:1
    ider = 2*k+ell; % -sign in (-dt)^ell
    sfac = rf(k)*A(k,ell)/(-dt)^ell; % (1/dt^(2*k)) already in rf(k)
    % zero padding with npextra points on each side
    if nsub>1
      sfac = sfac/nsub^ider;
    end
    % zero padding with npextra on both ends to enable for central differencing
    a2 = [zeros(1,npextra) (a .* t.^ell) zeros(1,npextra)];
    [da,nda] = apply_fd_op2(a2,ider,jorder,nhextra,sfac);
    ndmax = max(ndmax,nda); % internal check for stencil width (<=npextra?)
    da = da(npextra+1:end-npextra); % remove padding
    b = b+da;
  end
end
if ntstart>0, b = b(ntstart+1:end); end
if ndmax>npextra
  error('nd=%d > npextra=%d; increase initial npextra',ndmax,npextra);
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



