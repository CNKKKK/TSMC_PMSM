function [sys,x0,str,ts] = lh(t,x,u,flag)
v{1,1}=[1 0 0 1 0 0];v{1,2}=[1 0 0 0 0 1];v{1,3}=[1 1 0 0 0 0];
v{2,1}=[1 0 0 0 0 1];v{2,2}=[0 0 1 0 0 1];v{2,3}=[0 0 0 0 1 1];
v{3,1}=[0 0 1 0 0 1];v{3,2}=[0 1 1 0 0 0];v{3,3}=[0 0 1 1 0 0];
v{4,1}=[0 1 1 0 0 0];v{4,2}=[0 1 0 0 1 0];v{4,3}=[1 1 0 0 0 0];
v{5,1}=[0 1 0 0 1 0];v{5,2}=[0 0 0 1 1 0];v{5,3}=[0 0 0 0 1 1];
v{6,1}=[0 0 0 1 1 0];v{6,2}=[1 0 0 1 0 0];v{6,3}=[0 0 1 1 0 0];

s{1,1}=[1 0 0 1 0 1];s{1,2}=[1 0 1 0 0 1];
s{2,2}=[1 0 1 0 0 1];s{2,3}=[0 1 1 0 0 1];
s{3,3}=[0 1 1 0 0 1];s{3,4}=[0 1 1 0 1 0];
s{4,4}=[0 1 1 0 1 0];s{4,5}=[0 1 0 1 1 0];
s{5,5}=[0 1 0 1 1 0];s{5,6}=[1 0 0 1 1 0];
s{6,6}=[1 0 0 1 1 0];s{6,7}=[1 0 0 1 0 1];
s0=[0 1 0 1 0 1];s7=[1 0 1 0 1 0];
switch flag,
   case 0
    [sys,x0,str,ts]=mdlInitializeSizes(v,s,s0,s7) ;   
   case 3
    sys=mdlOutputs(t,x,u,v,s,s0,s7);
   case {1,2,4,9}
    sys=[];  
   otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes(v,s,s0,s7)
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 20;
sizes.NumInputs      = 8;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed
sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [1/100000 0];

function sys=mdlOutputs(t,x,u,v,s,s0,s7)
global yy I_alphar I_beta k j dm dn d1 d2 d0 h
global tt A0 A1 A2 P1 sita y1 y2 y3 y4 y5 y6 tp 
tp=0.001*6;

k=u(4);%0<k<1
j=u(8);%0<j<1
dm=0;
dn=0;
d1=0;
d2=0;
d0=0;
h=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

%fix on interval
 if u(1)>0&&u(2)<0&&u(3)<=0  
   yy=1; 
   dm=-u(2)*k/u(1);
   dn=-u(3)*k/u(1);
  
elseif u(1)>0&&u(2)>=0&&u(3)<0
    yy=2; 
    dm=-u(1)*k/u(3);
    dn=-u(2)*k/u(3);
    
elseif u(2)>0&&u(1)<=0&&u(3)<0
    yy=3;
    dm=-u(3)*k/u(2);
    dn=-u(1)*k/u(2);
    
elseif u(1)<0&&u(2)>0&&u(3)>=0
    yy=4;
    dm=-u(2)*k/u(1);
    dn=-u(3)*k/u(1);
    
elseif u(3)>0&&u(2)<=0&&u(1)<0
    yy=5;
    dm=-u(1)*k/u(3);
    dn=-u(2)*k/u(3);
   
 elseif u(3)>0&&u(2)<0&&u(1)>=0
    yy=6;
    dm=-u(3)*k/u(2);
    dn=-u(1)*k/u(2);
  
 end
 
 
 
I_alphar=(2/3)* (u(5)-(1/2)*u(6)-(1/2)*u(7));  
I_beta=(2/3)*((sqrt(3)/2)*u(6)-(sqrt(3)/2)*u(7));

 if  I_beta>0     
        A0=1;
    else
        A0=0;
 end
    if  0.866*I_alphar-0.5*I_beta>0    
        A1=1;        % sin(sita+2*pi/3)
    else
        A1=0;
    end    
    if  -0.866*I_alphar-0.5*I_beta>0     
        A2=1;       % sin(sita-2*pi/3)
    else
        A2=0;
    end
    P1=4*A2+2*A1+A0;

    if  P1==1            
        tt=2;
    elseif  P1==2        
        tt=6;
    elseif  P1==3        
        tt=1;
    elseif  P1==4       
        tt=4;
    elseif  P1==5        
        tt=3;
    elseif  P1==6       
        tt=5;
    else    
    end
 if (tt==1)
    
   sita=atan(I_beta/I_alphar);
   d1=sin(sita)*j;
   d2=sin(pi/3-sita)*j;
   d0=(1-d1-d2);
    s0_=s0;s7_=s7;
elseif (tt==2&&I_alphar>0) 
    sita=atan(I_beta/I_alphar);
    d1=sin(-pi/3+sita)*j;
    d2=sin(2*pi/3-sita)*j;
     d0=(1-d1-d2);
      s0_=s7;s7_=s0;
 elseif (tt==2&&I_alphar<=0) 
    sita=atan(I_beta/I_alphar);
    sita=sita+pi;
    d1=sin(-pi/3+sita)*j;
    d2=sin(2*pi/3-sita)*j;
     d0=(1-d1-d2);
      s0_=s7;s7_=s0;
elseif (tt==3)
    sita=atan(I_beta/I_alphar);
    sita=sita+pi;
    d1=sin(sita-pi*2/3)*j;
    d2=sin(pi-sita)*j;
    d0=(1-d1-d2);
     s0_=s0;s7_=s7;
elseif (tt==4)
    sita=atan(I_beta/I_alphar);
    sita=sita+pi;
    d1=sin(sita-pi)*j;
    d2=sin(4*pi/3-sita)*j;
    d0=(1-d1-d2);
     s0_=s7;s7_=s0;
elseif (tt==5&&I_alphar<=0)
    sita=atan(I_beta/I_alphar);
    sita=sita+pi;
    d1=sin(sita-4*pi/3)*j;
    d2=sin(5*pi/3-sita)*j;
    d0=(1-d1-d2);
     s0_=s0;s7_=s7;
elseif (tt==5&&I_alphar>0)
    sita=atan(I_beta/I_alphar);
    d1=sin(sita-4*pi/3)*j;
    d2=sin(5*pi/3-sita)*j;
    d0=(1-d1-d2);
     s0_=s0;s7_=s7;
 elseif (tt==6)
    sita=atan(I_beta/I_alphar);
    d1=sin(sita-5*pi/3)*j;
    d2=sin(2*pi-sita)*j;
    d0=(1-d1-d2);
     s0_=s7;s7_=s0;
 end
 % deal in pie slice
y1=d1*dm*tp/6;
y2=d2*dm*tp/6;
y3=d2*dn*tp/6;
y4=d1*dn*tp/6;
y5=(d0/2)*dm*tp/6;
y6=(d0/2)*dn*tp/6;
if rem(t,tp/6)>=0&rem(t,tp/6)<(y5/2)
     h=[v{yy,1} s0_ yy  tt I_alphar I_beta y3 y4 y5 y6];
elseif rem(t,tp/6)>=(y5/2)&rem(t,tp/6)<(y1+y5/2)
     h=[v{yy,1} s{tt,tt} yy  tt I_alphar I_beta y3 y4 y5 y6];                
elseif rem(t,tp/6)>=(y1+y5/2)&rem(t,tp/6)<(y1+y2+y5/2)
     h=[v{yy,1} s{tt,tt+1} yy  tt I_alphar I_beta y3 y4 y5 y6]; 
elseif rem(t,tp/6)>=(y1+y2+y5/2)&rem(t,tp/6)<(y1+y2+y5)
     h=[v{yy,1} s7_ yy  tt I_alphar I_beta y3 y4 y5 y6];
elseif rem(t,tp/6)>=(y1+y2+y5)&rem(t,tp/6)<(y1+y2+y5+y6/2)
     h=[v{yy,2} s7_ yy  tt I_alphar I_beta y3 y4 y5 y6];
elseif rem(t,tp/6)>=(y1+y2+y5+y6/2)&rem(t,tp/6)<(y1+y2+y5+y6/2+y4)
     h=[v{yy,2} s{tt,tt+1} yy  tt I_alphar I_beta y3 y4 y5 y6];
elseif rem(t,tp/6)>=(y1+y2+y5+y6/2+y4)&rem(t,tp/6)<(y1+y2+y5+y6/2+y4+y3)
     h=[v{yy,2} s{tt,tt} yy  tt I_alphar I_beta y3 y4 y5 y6]; 
elseif rem(t,tp/6)>=(y1+y2+y5+y6/2+y4+y3)&rem(t,tp/6)<(y1+y2+y5+y6+y4+y3)
     h=[v{yy,2} s0_ yy  tt I_alphar I_beta y3 y4 y5 y6];
elseif rem(t,tp/6)>=(y1+y2+y5+y6+y4+y3)&rem(t,tp/6)<(tp/6)
     h=[v{yy,3} s0_ yy  tt I_alphar I_beta y3 y4 y5 y6];
end

sys=[h];  