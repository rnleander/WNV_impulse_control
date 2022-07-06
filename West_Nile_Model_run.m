function [tt,x] = West_Nile_Model_run(larvicide_type)
%This function simulates a West Nile virus epidemic or mosquito population
% with no control. Disease-free or disease-free-equilibrium conditionns 
% are selected by the user prior to running on lines 56-60.

%Vector
%Es = x(1);              %eggs laid by susceptible and exposed mothers
%Ei = x(2);              %eggs laid by infected mothers
%Ls = x(3);              %susceptible larvae
%Li = x(4);              %infected larvae
%Vs = x(5);              %susceptible vectors
%Ve = x(6);              %exposed vectors
%Vi = x(7);              %infected vectors

%Host
%Hs = x(8);              %susceptible hosts
%Hi = x(9);              %infected hosts
%Hr = x(10);             %recovered hosts

%Control
%Ul = x(11);             %larvacide
%Ua = x(12);             %adultacide

%Artifical State
%int_0^t{cV*cV*(Vi(s)+Vs(s)+Ve(s)}ds=x(13)

%NH = Hs+Hi+Hr;           %total hosts


%The continuous state has 13 components. 
%These components track the value of the continous variables
%[Es,Ei,Ls,Li,Vs,Ve,Vi,Hs,Hi,Hr,Ul,Ua,Int] post treatment


%Model Parameters
%generate parameters
%Duration of simulation
Tf=150;
p = System_parametersRL(larvicide_type,Tf);

% Initial conditions for discrete/continuous state variables
%ic(2)=infected egges, %ic(4)=infected larva, %ic(7)=infected mosquitoes

%Vector parameters
rs = p(1);            %intrinsic rate of increase of uninfected mosquitoes
m_e = p(6);              %egg maturatation rate
m_l = p(7);             %larval maturation rate
muV=p(9);             %adult death rate
C = p(11);             %mosquito carrying capacity

ic_V=m_l*C/muV;
ic_E=rs*m_l*C/(m_e*muV);

NH=p(22);

% Initial conditions for discrete/continuous state variables
%healthy, summer ic. starts from DFE%
%ic = [ic_E;0;C;0;ic_V;0;0;NH;0;0;0;0;0];
%diseased
ic = [ic_E;0;C;0;ic_V;0;.01*ic_V;NH;0;0;0;0;0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



X01=ic;
%the derivative of the state with respect to its initial condition is intially one,
%so we have an identity matrix
X02=reshape(eye(length(X01)),length(X01)^2,1);
X0=[X01; X02];

f=@(t,x)West_Nile_ModelRL(t,x,p);

%solve the state equations
[tt,x]=ode45(f,[0,Tf],X0);

West_Nile_Model_plots(tt,x,0,[],[],Tf,[],[],larvicide_type);
end
