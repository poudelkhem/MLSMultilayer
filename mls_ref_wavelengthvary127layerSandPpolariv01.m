
%% MLS MULTILAYER DATE 09/18/2017 KHEM POUDEL,MTSU
% origianl code by Dr. Robertson
%plot the reflection coeffcient vs. Wavelength for 63 sequence MLS 
% and also fidn average reflection coefficient





% th = input('enter the incident angle (degrees): ');
% wl1 = input('enter the start wavelength (angs): ');
% wl2 = input('enter the end wavelength (angstroms): ');
% nth = input('enter the number of wavelength steps: ');
% %ep = input('enter the dielectric constant array: ');
% %d = input('enter the layer thickness array as [0, ..., 0]: ')

 clear all;close all;
 theta=linspace(0,75,40);
%th =80;% Incident angle
wl1 =3500;%start wavelength (angs):
wl2 =8000;%end wavelength (angs):
nth = 10000;%input('enter the number of angular steps: ');
delta = (wl2-wl1)/nth;
eta_L0=eta_Sio_2(550);%TiO20=eta_Sio_2(550);
    eta_H0=eta_Tio_2(550);%TiO20=eta_Tio_2(550);
    eta_L0=eta_L0;
    eta_H0=eta_H0;
    d_L=5500/(4*eta_L0);
    d_H=5500/(4*eta_H0);

REFS1=zeros(length(theta),nth);
for i=1:40
    th=theta(i);
    %disp(th);
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
    eps_L=eta_L^2;% Permittivity of Low Index material//63 MLS sequence
    eps_L=eps_L+0.0001*1i;
    ep=[1,eps_H,eps_H,eps_H,eps_H,eps_H,eps_H,eps_H,eps_L,eps_L,eps_L,eps_L,eps_H,eps_H,eps_H,eps_L,eps_H, ...  
         eps_H,eps_H,eps_H,eps_L,eps_L,eps_H,eps_L,eps_H,eps_H,eps_L,eps_L,eps_H,eps_L,eps_L,eps_H,eps_L, ...
         eps_L,eps_L,eps_L,eps_L,eps_L,eps_H,eps_L,eps_L,eps_L,eps_H,eps_L,eps_L,eps_H,eps_H,eps_L,eps_L, ...
         eps_L,eps_H,eps_L,eps_H,eps_H,eps_H,eps_L,eps_H,eps_L,eps_H,eps_H,eps_L,eps_H,eps_H,eps_L,eps_L,...
         eps_L,eps_L,eps_L,eps_H,eps_H,eps_L,eps_L,eps_H,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_L,eps_H,...
         eps_H,eps_H,eps_L,eps_L,eps_H,eps_H,eps_H,eps_H,eps_L,eps_H,eps_H,eps_L,eps_H,eps_L,eps_L,eps_L,...
         eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_L,eps_H,eps_H,eps_H,eps_H,eps_H,eps_L,eps_H,eps_L,eps_L,... 
         eps_H,eps_L,eps_H,eps_L,eps_L,eps_L,eps_H,eps_H,eps_L,eps_H,eps_H,eps_H,eps_L,eps_L,eps_L,2.25];    
%     ep=[1,eps_H,eps_H,eps_H,eps_H,eps_H,eps_H,eps_L,eps_L,eps_L,eps_L,eps_L,eps_H,eps_L,eps_L,eps_L,eps_L, ...
%           eps_H,eps_H,eps_L,eps_L,eps_L,eps_H,eps_L,eps_H,eps_L,eps_L,eps_H,eps_H,eps_H,eps_H,eps_L,eps_H, ...
%           eps_L,eps_L,eps_L,eps_H,eps_H,eps_H,eps_L,eps_L,eps_H,eps_L,eps_L,eps_H,eps_L,eps_H,eps_H, eps_L,...
%           eps_H,eps_H,eps_H,eps_L,eps_H,eps_H,eps_L,eps_L,eps_H,eps_H, eps_L, eps_H, eps_L, eps_H, eps_L,2.25];
%    
% d_H=610;% Optical Thickeness d_H=61.08nm at 632 nm
% d_L=1082;% Optical Thickeness d_L=108.57nm at 632 nm

d=[0,d_H,d_H,d_H,d_H,d_H,d_H,d_H,d_L,d_L,d_L,d_L,d_H,d_H,d_H,d_L,d_H,...  
     d_H,d_H,d_H,d_L,d_L,d_H,d_L,d_H,d_H,d_L,d_L,d_H,d_L,d_L,d_H,d_L, ...  
     d_L,d_L,d_L,d_L,d_L,d_H,d_L,d_L,d_L,d_H,d_L,d_L,d_H,d_H,d_L,d_L, ...  
     d_L,d_H,d_L,d_H,d_H,d_H,d_L,d_H,d_L,d_H,d_H,d_L,d_H,d_H,d_L,d_L,...  
     d_L,d_L,d_L,d_H,d_H,d_L,d_L,d_H,d_H,d_L,d_H,d_L,d_H,d_L,d_L,d_H,...
     d_H,d_H,d_L,d_L,d_H,d_H,d_H,d_H,d_L,d_H,d_H,d_L,d_H,d_L,d_L,d_L,...  
     d_L,d_H,d_L,d_H,d_L,d_H,d_L,d_H,d_H,d_H,d_H,d_H,d_L,d_H,d_L,d_L,...
     d_H,d_L,d_H,d_L,d_L,d_L,d_H,d_H,d_L,d_H,d_H,d_H,d_L,d_L,d_L,0];    

	[rp(k),rs(k),tp(k),ts(k)] = ref(wl,ep,d,th);
	REFP(k) = abs(rp(k))^2;
	REFS(k) = abs(rs(k))^2;
    REFS1(i,k)=REFS(k);
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
    
    
    
%  plot(WL/10,REFS,WL/10,REFS,'LineWidth',3,'MarkerSize',8)
 RefAvgS(i)=sumS/nth;
 RefAvgP(i)=sumP/nth;
%  disp(sumS/nth);
%  disp(sumP/nth);
 
 
 
 % Read the data "REFSAltMultilayer_Ang0Deg.txt" that was saved on previous
 % run of Alternate 64 layer
% Multilayer_Ang0=load ("REFSAltMultilayer_Ang0Deg.txt");
% figure(1);
% plot(WL/10,REFP,'LineWidth',3,'MarkerSize',8)
% % plot(WL/10,REFS,'b',WL/10,Multilayer_Ang0,'r','LineWidth',3,'MarkerSize',8)
% xlabel('Wavelength [nm]')
%  ylabel('Reflection Coeff.');
%  set(findall(gcf,'type','text'),'FontSize',18);
%  set(gca,'YDir','normal');
%  axis tight;
%  hold on;
end
%  legend(' 0^{\circ}',' 30^{\circ}',' 60^{\circ}');
% legend('0^{\circ} P');
% % legend('MLS' ,'Alt.High Low');
% set(gca,'fontsize',18)

figure(2);
% REFSAltMultilayer=load ("REFSAltMultilayer.txt");
% REFPAltMultilayer=load ("REFPAltMultilayer.txt");
%  plot(theta,RefAvgS,'b--',theta,RefAvgP,'b',theta,REFSAltMultilayer,'r--',theta,REFPAltMultilayer,'r','LineWidth',3,'MarkerSize',8)
plot(theta,RefAvgS,'b',theta,RefAvgP,'r','LineWidth',3,'MarkerSize',8)
 legend('S Pol.','P Pol.')
 xlabel('Incident angle (Deg)')
 ylabel('Average Refln. Coeff.');
 ylim([0.75 1.0])
 xlim([0.0 75])
 set(findall(gcf,'type','text'),'FontSize',18);
 set(gca,'YDir','normal');
 grid on;
 %axis tight;
 ax = gca;
ax.GridColor = [0 .1 .1];
% ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;
ax.Layer = 'top';
 set(gca,'fontsize',18)


% figure(2);
% surf(REFS1);
% colorbar
% figure(2);
% surf(REFS1)
% colorbar('Direction','reverse')