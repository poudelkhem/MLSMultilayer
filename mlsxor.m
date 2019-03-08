N=input('Input the order of the MLS sequence, N (Choose from N=3,4,5,6,7,9,10,11,15):')
if N==3; x0=[-1 -1 -1]; taps=[1 2]; end
if N==4; x0=[-1 -1 -1 -1]; taps=[1 2]; end
% if N==5; x0=[-1 -1 -1 -1 -1]; taps=[1 3]; end
if N==5; x0=[-1 -1 -1 -1 -1]; taps=[1 3]; end
if N==6; x0=[-1 -1 -1 -1 -1 -1]; taps=[1 2]; end
if N==7; x0=[-1 -1 -1 -1 -1 -1 -1]; taps=[1 4]; end
if N==9; x0=[-1 -1 -1 -1 -1 -1 -1 -1 -1]; taps=[1 5]; end
if N==10; x0=[-1 -1 -1 -1 -1 -1 -1 -1 -1 -1]; taps=[1 4]; end
if N==11; x0=[-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]; taps=[1 3]; end
if N==15; x0=[-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]; taps=[1 2]; end
M=max(size(x0));
t1=taps(1);
t2=taps(2);
for k=1:((2^M)-1)
% for 1 and 0 instead of 1 and -1 use commented line
    %    xn=xor(x0(t1),x0(t2))
    xn=x0(t1)*x0(t2);
    x0=[x0 xn];
    t1=t1+1;
    t2=t2+1;
end
x02=x0;
for k=1:max(size(x02))
if(x02(k)<0)
x02(k)=0;
end
end
x02