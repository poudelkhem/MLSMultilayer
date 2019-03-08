
function [refp,refs,tefp,tefs] = ref(wl,ep,d,th)
%
% FUNCTION THAT CALCULATES THE REFLECTIVITY OF A MULTILAYER
% wl - wavelength in Angstroms
% ep - dielectric constant of each layer
% d  - thickness in Angstroms of each layer
% th - incident angle in degrees
% BOTH p- AND s- POLARIZED REFLECTIVITIES ARE CALCULATED
%




nlayer = max(size(ep));
th = th*0.01745329;
sinth1 = sin(th);
wki = 2*pi*1e8/wl;

for k=1:nlayer
	costh(k) = sqrt(1.0 - ((ep(1)/ep(k))*sinth1^2));
	wvi(k) = wki*sqrt(ep(k));
end

for m=1:(nlayer-1)
	rp(m) = (sqrt(ep(m+1))*costh(m)-sqrt(ep(m))*costh(m+1))/ ...
		(sqrt(ep(m+1))*costh(m)+sqrt(ep(m))*costh(m+1));
	rs(m) = (sqrt(ep(m))*costh(m)-sqrt(ep(m+1))*costh(m+1))/ ...
		(sqrt(ep(m))*costh(m)+sqrt(ep(m+1))*costh(m+1));
	tp(m) = (2*sqrt(ep(m))*costh(m))/(sqrt(ep(m+1))*costh(m)+ ...
      sqrt(ep(m))*costh(m+1));
    ts(m) =(2*sqrt(ep(m))*costh(m))/(sqrt(ep(m))*costh(m)+...
      sqrt(ep(m+1))*costh(m+1));
end

rrp(nlayer-1) = rp(nlayer-1);
rrs(nlayer-1) = rs(nlayer-1);
ttp(nlayer-1) = tp(nlayer-1);
tts(nlayer-1) = ts(nlayer-1);

for n=nlayer-2:-1:1
	decay = exp(2*i*d(n+1)*costh(n+1)*wvi(n+1)*1e-8);
	decayt = exp(i*d(n+1)*costh(n+1)*wvi(n+1)*1e-8);
	rrp(n)=(rp(n)+rrp(n+1)*decay)/(1.0 + rp(n)*rrp(n+1)*decay);
	rrs(n)=(rs(n)+rrs(n+1)*decay)/(1.0 + rs(n)*rrs(n+1)*decay);
    ttp(n)=(tp(n)*ttp(n+1)*decayt)/(1.0 + rp(n)*rrp(n+1)*decay);
    tts(n)=(ts(n)*tts(n+1)*decayt)/(1.0 + rs(n)*rrs(n+1)*decay);
end
refp = rrp(1);
refs = rrs(1);
tefp = ttp(1);
tefs = tts(1);
