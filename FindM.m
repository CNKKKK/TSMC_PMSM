function [sys,x0,str,ts] = FindM(t,x,u,flag)
switch flag,
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 1,
    sys=[];
  case 2,
    sys=[];
  case 3,
    sys = mdlOutputs(t,x,u);
  case 4,
    sys=[];
  case 9,
    sys=[];
   otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end
%=============================================================================
% mdlInitializeSizes
function [sys,x0,str,ts]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 4;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [1/100000 0];     % 采样时间1us,10M Hz

%=============================================================================
% mdlOutputs
function sys=mdlOutputs(t,x,u)
global Udc ULA ULB ULC ULM

    Udc=u(1);           
    ULA=u(2);           
    ULB=u(3);           
    ULC=u(4);           
    if((ULA>=ULB)&&(ULA>=ULC))
        ULM=ULA;
    elseif((ULB>=ULA)&&(ULB>=ULC))
        ULM=ULB;
    elseif ((ULC>=ULA)&&(ULC>=ULB))
        ULM=ULC;
    end
    sys(1)=ULM;
sys(2)=ULM/Udc;



