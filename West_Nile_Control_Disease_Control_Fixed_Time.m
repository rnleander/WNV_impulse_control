function [tau,ul,ua,X,J,J_comp,final_treatment_time,X0,T] = West_Nile_Control_Disease_Control_Fixed_Time(N,Tf,larvicide_type)
%This code optimizes the control level
%in order to minimize the disease burden plus control cost.
%larvicide_type: 1=long-lasting s-methorpene briquet, 2=VectoBac
%N is the maximal number of treatments
%Tf is the duration of the control
%Initial conditions are set on line 93 as the disease-free equilibrium plus  
%a low density of infectious mosquitoes.

%State variables

%Vector
%Es = x(1);               %eggs laid by susceptible and exposed mothers
%Ei = x(2);               %eggs laid by infected mothers
%Ls = x(3);              %susceptible larvae
%Li = x(4);              %infected larvae
%Vs = x(5);              %susceptible vectors
%Ve = x(6);              %exposed vectors
%Vi = x(7);              %infected vectors

%Host
%Hs = x(8);              %susceptible hosts
%Hi = x(9);              %infected hosts
%Hr = x(10);              %recovered hosts

%Control
%Ul = x(11);              %larvacide
%Ua = x(12);             %adultacide

%Artifical State
%int_0^t{cV*(Vi(s)+Hi(s))}ds=x(13)

%NH = Hs+Hi+Hr;           %total hosts

%%% Impulse equation
%%%%% lower case u is for the impulse of adulticide and larvacide
% Ul(Ti^+) = Ul(Ti^+)+ul(i);
% Ua(Ti^+) = Ua(Ti^+)+ua(i);


%The discrete state has 13 components. These components track the value of the continous variables
%[Es,Ei,Ls,Li,Vs,Ve,Vi,Hs,Hi,Hr,Ul,Ua,Int] post treatment.
%The discrete state is denoted by X(i,:)=[X(i,1), . . ., X(i,13)]

%J_comp=cv*X(13,N)+ce*(X(2,N))+cl*sum(ul.^2)+ca*sum(ua.^2)
%where X(13,N) is the value of the artifical state at the final time and 
% X(2,N) gives the density of infectious-laid eggs at the final time. 


%The Hamiltonian is -J+sum_{i=1}^{N}{<Y(i),(G(X(i-1),ul(i),ua(i))-X(i))>} + <Y(0),X0-X(0)>
% + sum_{i=1}^N{mu_a^+(i)*(u_a(i)-1)-mu_l^+(i)*(u_l(i)-1)-mu_a^-(i)*u_a(i)-mu_l^-(i)*u_l(i)}
%Here the vector valued function G determines how the discrete state at the next time is determined by
%the discrete state at the previous time. Y(i) has one component for each state component.
% <Y,G> denotes the dot product.
%The second line of the Hamiltonian contains the inequality constraints. The coefficients here are <=0.
%In addition sum_{i=1}^N{mu_a^+(i)*(u_a(i)-1)-mu_l^+(i)*(u_l(i)-1)-mu_a^-(i)*u_a(i)-mu_l^-(i)*u_l(i)}=0
                                              
%let dx/dt=g(x(t,x0)), where g is vector-valued and x depends on the initial state x0.
%G_i(x0,ul,ua,tau)=x0_i+int_0^tau{g_i(x(t,x0)) dt} for i=1, . . . 10
%G_i(x0,ul,ua,tau)=x0_i+int_0^tau{g_i(t,x0) dt} + ul(i) for i=11
%G_i(x0,ul,ua,tau)=x0_i+int_0^tau{g_i(t,x0) dt} + ua(i) for i=12
%G_i(x0,ul,ua,tau)=x0_i+int_0^tau{g_i(x(t,x0)) dt} for i=13


test=-1;

delta=10^(-8);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Model Parameters
%generate parameters
p = System_parametersRL(larvicide_type, Tf);

%Vector
rs = p(1);            %intrinsic rate of increase of uninfected mosquitoes
m_e = p(6);              %egg maturatation rate
m_l = p(7);             %larval maturation rate
muV=p(9);             %adult death rate

C = p(11);             %mosquito carrying capacity

ic_V=m_l*C/muV;
ic_E=rs*m_l*C/(m_e*muV);

NH=p(22);

% Initial conditions for discrete/continuous state variables
%disease free steady 
%ic = [ic_E;0;C;0;ic_V;0;0;1;0;0;0;0;0;0];
%start with a low density of infectious mosquitoes. 
ic = [ic_E;0;C;0;ic_V;0;.01*ic_V;NH;0;0;0;0;0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Duration of treatment
%Tf
%weight of cost of larvacide
cl = p(23);
%weight of cost of adulticide
ca = p(24); %adulticide treatment is more expensive. Cite Mosquitoes and disease Illinois dept of public health
cT=p(25);
%cost of eggs at the final time
ce=p(26);

cv=p(21);

%maximum time between controls
Maxt=p(28);
%minimum time between controls
mint=p(29);


%state variables with underlying continuous dynamics
X01=ic(1:13);
%the derivative of the state with respect to its initial condition is intially one,
%so we have an identity matrix
X02=reshape(eye(length(X01)),length(X01)^2,1);
X0=[X01; X02];


%equally-spaced waiting times with first treatment at t=0 and last
%treatment at t=T_f
tau=(Tf/(N-1))*ones(N,1);
tau(1)=0;


T=zeros(size(tau));

ua=zeros(N,1);
ul=zeros(N,1);

%vector of discrete states after control and derivatives before treatment
X=zeros(length(X0),N);
%discrete states just before control is added
XX=zeros((length(X0)),N);
%derivative of discrete state with respect to previous state value (initial value).
DX=zeros(length(X01),length(X01),N);
%discrete adjoint variables
Y=zeros(N,length(ic));


f=@(t,x)West_Nile_ModelRL_Disease(t,x,p);

J_values=[];

count = 0;

while (test<0)
    
    oldul=ul;
    oldua=ua;   %%% control variable
    oldX=X;
    
%This loop finds the discrete states just before (XX) and just after (X) each treatment
for i=1:N
        if i==1
            x0=X0;
            x=X0;
            T(i)=tau(i);
            %state just before treatment
            XX(:,i) = x(:);
        
            X(:,i) = x(:);
    
            %states after the addition of the next dose
            X(11,i)=X(11,i)+ul(i);
            X(12,i)=X(12,i)+ua(i);
        else
            T(i)=T(i-1)+tau(i);
            x01=X(1:13,i-1);
            %the derivative of the state with respect to its initial condition is intially one,
            %so we have an identity matrix
            x02=reshape(eye(length(x01)),length(x01)^2,1);
            x0=[x01; x02];
            options = odeset('RelTol',1e-8,'AbsTol',1e-10);
            sol=ode45(f,[0,tau(i)],x0,options);
            ll=length(sol.x);
            x=sol.y(:,ll);
            %state just before treatment
            XX(:,i) = x(:);
        
            X(:,i) = x(:);
    
            %states after the addition of the next dose
            X(11,i)=X(11,i)+ul(i);
            X(12,i)=X(12,i)+ua(i);
        end   
    DX(:,:,i)=reshape(x(length(X01)+1:length(X01)^2+length(X01)),length(X01),length(X01));
    %DX=[dx1/dx1(0), dx1/dx2(0),. . .dx1/dx13(0);
       % dx2/dx1(0), dx2/dx2(0), . . .dx2/dx13(0);
       %    .
       %    .
       %    .
       % dx13/dx1(0), dx13/dx2(0), . . .,dx13/dx13(0)];
end

%get the adjoint variables at steps 1,. . .,N.
Y(N,:)= [0 -ce 0 0 0 0 0 0 0 0 0 0 -cv];
for i=1:N-1
    j=N-i;
    Y(j,1:length(X01))=Y(j+1,1:length(X01))*DX(:,:,j);
end

for i=1:N
        
    ul(i) = Y(i,11)/(2*cl);
    ua(i) =Y(i,12)/(2*ca);

    ul(i)=max(0,min(1,ul(i)));
    ua(i)=max(0,min(1,ua(i)));
    
end
    
    J=cv*X(13,N)+ce*(X(2,N))+cl*sum(oldul.^2)+ca*sum(oldua.^2);
    J_values=[J_values J];

    %evaluate total relative error
    testua=delta*sum(abs(ua))-sum(abs(ua-oldua));
    testul=delta*sum(abs(ul))-sum(abs(ul-oldul));
    testX=delta*sum(abs(X(1:13,:)))-sum(abs(X(1:13,:)-oldX(1:13,:)));
    
    %update control
    ul=(.9*oldul+.1*ul);
    ua=(.9*oldua+.1*ua);
        
    test=min([testua testul testX]);
    count=count+1;
    
    final_treatment_time=sum(tau);
end

J_comp=cv*X(13,N)+ce*(X(2,N))+cl*sum(oldul.^2)+ca*sum(oldua.^2)+cT*sum(tau.^2);

M=length(J_values);
figure
plot(1:M,J_values,'*') 
ylabel('Objective functional value','FontSize', 20)
xlabel('iteration','FontSize', 20)
set(gca,'fontsize',16)

file_name=sprintf('J_fixed_times_T=%.2f_N=%.2f.eps',Tf,N);

figure_title=sprintf('N=%.2f T=%.2f',N,Tf);

title(figure_title)

%legend({'Objective functional value','Final time'},'Location','best','FontSize', 20)
exportgraphics(gcf,file_name)

hold off

%This part of the code plots the solution
x=[];
t=[];
for i=1:N
        if i==1
        x0=X0;
        else
            x0=X(:,i-1);
        end
        %if not time passes, we do not need to add any points to the
        %solution graph
    if tau(i)~=0
        %solve the state equations forward in time
        [tt,xx]=ode45(f,[0,tau(i)],x0);
        if i>1
        tt=tt+T(i-1);
        end
        t=cat(1,t,tt);
        x=cat(1,x,xx);  
    end
    
end

control_type=1;
Obj_type=2;

West_Nile_Model_plots(t,x,control_type,Obj_type,N,Tf,J,J_comp)
    
end