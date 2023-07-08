function [b,n,w] = apply_fd_op2(a,ider,iorder,nhextra,sfactor)
% [b,nskipped,w] = apply_fd_op2(a,ider,iorder,nhextra,sfactor)
% missing factor 1/spacing^ider has to be applied afterwards

if nargin<5, sfactor = 1; end
if nargin<4, nhextra = 0; end

if ider == 0
  n = 0; b = a; if sfactor~=1, b = sfactor*b; end
  return; 
end
[w,n] = fd1_op2(ider,iorder,nhextra);
if sfactor~=1, w = sfactor*w; end

b = zeros(size(a));
if mod(ider,2) == 1
  for j=n+1:length(a)-n
    b(j) = dot( w, a(j+(1:n)) - a(j-(1:n)) );
  end
else
  for j=n+1:length(a)-n
    b(j) = dot( w, (a(j+(1:n))-a(j)) + (a(j-(1:n))-a(j)) );
  end
end

end
% _____________________________________________________________________________
function [b,n,bnum,bden] = fd1_op2(ider,iorder,nhextra)
% [b,n,bnum,bden] = fd1_op2(ider,iorder,nhextra)
% ider-th derivative operator of order iorder with nhextra points on each side
% b = bnum/bden
persistent ww2 hordermax nextramax fname;

if nargin<2 || iorder <3, iorder = 2;
else, iorder = 2*ceil(iorder/2); % should be even
end

% file with all results, pre-allocated with
% nextramax such that half the number of extra stencil points 
% nhextra <= nextramax for all examples considered and 
% hordermax is half of the maximum order of dt
% will be resized if necessary
if isempty(fname), fname = 'loadallww.m'; end
if isempty(ww2) && exist(fname,'file')
  eval( strrep(fname,'.m','') ); % fill ww2 with contents for fname
  if ~isempty(ww2)
    nww2 = size(ww2); hordermax = nww2(2); nextramax = nww2(3);
  end
end
if ~exist(fname,'file')
  nextramax = 2; hordermax = 4;
  local_mkfile(fname,hordermax,nextramax);
elseif isempty(nextramax) || isempty(hordermax)
  [hordermax,nextramax] = local_file_getsize(fname);
end

if isempty(ww2) % pre-allocate
  ww2 = cell(3*hordermax,hordermax,nextramax+1);
end
if ~isempty(ww2) % check size
  if ider>size(ww2,1) || iorder/2>size(ww2,2) || nhextra+1>size(ww2,3)
    hordermax = max([ceil(ider/3),ceil(iorder/2),hordermax]);
    nextramax = max([nhextra+1,nextramax]);
    if exist(fname,'file')
      local_file_changesize(fname,hordermax,nextramax); 
    end
    eval( strrep(fname,'.m','') );  
  end
end

if isempty( ww2{ider,iorder/2,nhextra+1} ) 
  % derive operator and update file fname
  bn = set_fd1_op2(ider,iorder,nhextra,fname);
  ww2{ider,iorder/2,nhextra+1} = bn;
else
  bn = ww2{ider,iorder/2,nhextra+1};
end

b = bn{1}; n = bn{2}; bnum = bn{3}; bden = bn{4}; 
end
% _____________________________________________________________________________
function wn = set_fd1_op2(ider,iorder,nhextra,fnamedif)
%  w = set_fd1_op2(ider,iorder)
%  difference operator is
%     [-fliplr(w)      0  w] for odd ider and
%     [ fliplr(w) -sum(w) w] for even ider

[w,nw,wsym,wnum,wden0,iorder2] = get_fd_op2(ider,iorder,nhextra);
w = w(nw+2:end); ws = wnum(nw+2:end); % use symmetry and ignore central value

append_ww2(fnamedif,w,ws,wden0,ider,iorder,nhextra);
wn = {w,nw,ws,wden0};

end
% _____________________________________________________________________________
function append_ww2(fname,w,ws,wden,ider,iorder,nhextra)
% append_ww2(fname,w, ws,wden, ider,iorder,nhextra)
if isempty(fname)
  fp = 1; % stdout
elseif ~exist(fname,'file')
  error('file %s does not exist',fname);
else 
  fp = fopen(fname,'a');
end
np = length(ws);
if iorder==2
  fprintf(fp,'%% _____________________________________________________________________________\n');
end
fprintf(fp,'ww2{%d,%d,%d}={[',ider,iorder/2,1+nhextra);
for k=1:np-1, fprintf(fp,'%.17g ',w(k)); end
fprintf(fp,'%.17g],%d,[',w(np),np);
for k=1:np-1, fprintf(fp,'%s ',char(ws(k))); end
fprintf(fp,'%s],%s};\n',char(ws(np)),char(wden));
if ~isempty(fname), fclose(fp);end

end
% _____________________________________________________________________________
function local_mkfile(fname,hordermax,nextramax)
% local_mkfile(fname,hordermax,nextramax)
% create file with header
if exist(fname,'file'), warning('%s already exists',fname); return; end

fp = fopen(fname,'w');
fprintf(fp,'%% %s\n%% ww{derivative,order/2,nhextra+1} contains half the operator and its length\n%%   and split into numerator, denominator\n',...
  fname);
fprintf(fp,'ww2 = cell(%d,%d,%d);\n%%\n',3*hordermax,hordermax,nextramax);
fclose(fp);
end
% _____________________________________________________________________________
function [hordermax,nextramax] = local_file_getsize(fname)
% [hordermax,nextramax] = local_file_getsize(fname)
%
hordermax = 0; nextramax = 0; if isempty(fname), return; end
a = readlines(fname); ip = contains(a,'cell'); if ~any(ip), return; end
aip = a(ip); ip2=strfind(aip,'cell'); b = sscanf(aip,"ww2 = cell(%d,%d,%d)");
if length(b)==3
  hordermax = b(2); nextramax = b(3);
end
end
% _____________________________________________________________________________
function local_file_changesize(fname,hordermax,nextramax)
% local_file_changesize(fname,hordermax,nextramax)
if isempty(fname), return; end
a = readlines(fname); 
if isempty(a), warning('problem with file: %s',fname); end
ip = contains(a,'cell'); ip = find(ip); 
if isempty(ip), warning('file % has no line with cell()?',fname); return; end
a(ip(1)) = sprintf('ww2 = cell(%d,%d,%d);',3*hordermax,hordermax,nextramax);
fp = fopen(fname,'w');
for j=1:length(a), fprintf(fp,'%s\n',a(j)); end
fclose(fp);
end
% _____________________________________________________________________________
%EOF
