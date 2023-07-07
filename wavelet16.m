function w = wavelet16(t,Tw,ider)
% w = wavelet16(t,Tw,ider)
%   w = w(t) = [4 (t/Tw) (1-t/tw)]^16 if ider=0 and 0<t<Tw
%   w = w'(t) if ider=1
%   0 outside interval (0,Tw)
if nargin<3, ider=0; end
w = 0*t;
th = t/Tw;
ih = find(th>0 & th<1); if isempty(ih), return; end
th = th(ih);
if ider==0
  w(ih) = ( 4*th.*(1.-th) ).^16;
elseif ider==1
  w(ih) = (64/Tw) * (1-2*th) .* ( 4*th.*(1.-th) ).^15;
else
  error('ider=%d',ider);
end

end
%EOF
