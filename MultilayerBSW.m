
%% MLS MULTILAYER DATE 09/18/2017 KHEM POUDEL,MTSU

 clear all;close all;
 theta=linspace(0,70,20);
 disp(theta);
wl1 =6320;%start wavelength (angs):
wl2 =6320;%end wavelength (angs):
nth = 100000;%input('enter the number of angular steps: ');
delta = (wl2-wl1)/nth;
eta_L0=eta_Sio_2(640);%TiO20=eta_Sio_2(550);
    eta_H0=eta_Tio_2(640);%TiO20=eta_Tio_2(550);
    eta_L0=eta_L0;
    eta_H0=eta_H0;
    d_L=1824.0;
    d_H=1173.0;
    d_L1=1565.0;
   
for k=1:20
    th=theta(k);
    disp(th);
    sum=0;
for i=1:1
	wl = wl1;
    %disp(wl)
    eta_H=eta_Tio_2(wl/10);%TiO2
    %disp(abs(eta_H))
    %eta_H=2.59;
    eps_H=eta_H^2;% Permittivity of High Index material
    eps_H=eps_H+0.0007*1i
    eta_L=eta_Sio_2(wl/10);%%SiO2
    %disp(abs(eta_L))
    %eta_H=1.457;
    eps_L=eta_L^2;% Permittivity of Low Index material
     eps_L=eps_L+0.0001*1i
    ep=[2.25,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,...
        eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L];
d=[d_L*4,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,...
        d_H,d_L,d_H,d_L,d_H,d_L, d_H,d_L,d_H,d_L1];

%     ep=[1,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L, 2.25];
%     d=[0,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L, 0];
    
% d_H=610;% 
% d_L=1082;


	[rp(k),rs(k),tp(k),ts(k)] = ref(wl,ep,d,th);
	REFP(k) = abs(rp(k))^2;
	REFS(k) = abs(rs(k))^2;
    sum=REFS(k)+sum;
    TEFP(k) = abs(tp(k))^2;
    TEFS(k) = abs(ts(k))^2;
	WL(k) = wl;
	PHAP(k) = 57.2957795*atan2((imag(rp(k))),(real(rp(k))));
	PHAS(k) = 57.2957795*atan2((imag(rs(k))),(real(rs(k))));
	if(k>2)
		dw=(2*pi*299792458*2*delta)/((wl-2*delta)^2-(wl-2*delta)*2*delta);
		GPDLYP(k-1)=(PHAP(k)-PHAP(k-2))*0.01745/dw;
        GPDLYP(k-1)=(PHAS(k)-PHAS(k-2))*0.01745/dw;
	end
end
	GPDLY(1) = 0.0;
	GPDLY(nth) = 0.0;
    
  save REFSBSW632.txt REFS -ascii  
    

end
% legend('0 Deg','15 Deg','30 Deg','45 Deg','60 Deg','75 Deg');
figure(1);
 
 plot(theta,REFS,'LineWidth',3,'MarkerSize',8)
 disp(sum/nth)
 
 xlabel('Wavelength [nm]')
 ylabel('Reflection Coeff.');
 set(findall(gcf,'type','text'),'FontSize',16);
 set(gca,'YDir','normal');
 axis tight;
%  figure(2);
%REFSBSW625=load ("REFSBSW625.txt");
% REFSBSW632=load ("REFSBSW632.txt");
%  plot(theta,REFSBSW625,'b--',theta,REFSBSW632,'r',theta,REFS,'k','LineWidth',3,'MarkerSize',8)
%  legend('625 nm','632.8 nm','640 nm')
%  xlabel('Incident angle (Deg)')
%  ylabel('Average Refln. Coeff.');
%  set(findall(gcf,'type','text'),'FontSize',18);
%  set(gca,'YDir','normal');
%  axis tight;