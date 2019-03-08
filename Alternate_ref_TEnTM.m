
%% MLS MULTILAYER DATE 09/18/2017 KHEM POUDEL,MTSU
% origianl code by Dr. Robertson
%plot the reflection coeffcient vs. Wavelength





% th = input('enter the incident angle (degrees): ');
% wl1 = input('enter the start wavelength (angs): ');
% wl2 = input('enter the end wavelength (angstroms): ');
% nth = input('enter the number of wavelength steps: ');
% %ep = input('enter the dielectric constant array: ');
% %d = input('enter the layer thickness array as [0, ..., 0]: ')

 clear all;close all;
 theta=linspace(0,75,11);
 disp(theta);
%th =80;% Incident angle
wl1 =4000;%start wavelength (angs):
wl2 =7000;%end wavelength (angs):
nth = 10000;%input('enter the number of angular steps: ');
delta = (wl2-wl1)/nth;

eta_L0=eta_Sio_2(550);%TiO20=eta_Sio_2(550);
    eta_H0=eta_Tio_2(550);%TiO20=eta_Tio_2(550);
    eta_L0=eta_L0;
    eta_H0=eta_H0;
    d_L=5500/(4*eta_L0);
    d_H=5500/(4*eta_H0);


for i=1:11% 60 degree only
    th=theta(i);
    disp(th);
   sumS=0;
    sumP=0;
for k=1:nth
	wl = wl1 + (k-1)*delta;
    %disp(wl)
    eta_H=eta_Tio_2(wl/10);%TiO2
    %disp(abs(eta_H))
    %eta_H=2.59;
    eps_H=eta_H^2;% Permittivity of High Index material
    eps_H=eps_H+0.0007*1i;
    eta_L=eta_Sio_2(wl/10);%%SiO2
    %disp(abs(eta_L))
    %eta_H=1.457;
    eps_L=eta_L^2;% Permittivity of Low Index material
    eps_L=eps_L+0.0001*1i;
    
        ep=[1,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,...
        eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L, 2.25];
    
        d=[0,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,...
        d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L, 0];
%     ep=[1,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,...
%         eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,...
%         eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,...
%         eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,2.25];
%     
%     d=[0,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,...
%     d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_L,0];
    %disp(ep)
%     d_H=1900/eta_H;% Optical Thickeness
%     d_L=600/eta_L;
% d_H=1579.9/eta_H;% Optical Thickeness d_H=61.08nm at 632 nm
%     d_L=1581/eta_L;% Optical Thickeness d_L=108.57nm at 632 nm
% d_H=610;% Optical Thickeness d_H=61.08nm at 632 nm
%     d_L=1082;% Optical Thickeness d_L=108.57nm at 632 nm
%input('enter the dielectric constant array: ');
%ep=[eps_H,eps_H,eps_H,eps_H,eps_H,eps_L,eps_L,eps_H,eps_H,eps_L,eps_H,eps_L,eps_L,eps_H,eps_L,eps_L,eps_L,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_H,eps_H,eps_L,eps_H,eps_H,eps_L,eps_L,eps_L];
%d = input('enter the layer thickness array as [0, ..., 0]: ')


	[rp(k),rs(k),tp(k),ts(k)] = ref(wl,ep,d,th);
	REFP(k) = abs(rp(k))^2;
	REFS(k) = abs(rs(k))^2;
    sumS=REFS(k)+sumS;
    sumP=REFP(k)+sumP;
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
    RefAvgS(i)=sumS/nth;
   RefAvgP(i)=sumP/nth;
   disp(sumS/nth);
   disp(sumP/nth);
    
  save REFSAltMultilayer_Ang0Deg.txt REFS -ascii  
  save REFSAltMultilayer.txt RefAvgS -ascii  
  save REFPAltMultilayer.txt RefAvgP -ascii  
    figure(1);
 plot(WL/10,REFS,'LineWidth',3,'MarkerSize',8)
% disp(sum/nth);
 
 xlabel('Wavelength [nm]')
 ylabel('Reflection Coeff.');
 set(findall(gcf,'type','text'),'FontSize',16);
 set(gca,'YDir','normal');
 axis tight;
 hold on ;
end
% legend('0 Deg','15 Deg','30 Deg','45 Deg','60 Deg','75 Deg');
legend('0 Deg');