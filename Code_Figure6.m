% Filename: Code_Figure6.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Jan 2021
% Reference:  Y. Li, S.T. Johnston, P.R. Buenzli, P. van Heijster, M.J. Simpson (2021) 
% Dimensionality affects extinction of bistable populations.
% The script contains:
%   - one call to the function Export_columndensity1D_Dis to generate 
%     Figures 6(d), the averaged density at time t in discrete simulations.
%   - one call to the function Export_columndensity1D_Num to generate 
%     Figures 6(d), the average column density at t in numerical solutions.
%   - one call to the function Export_totaldensity_Dis to generate 
%     Figures 6(e), the total population density in discrete simulations.
%   - one call to the function Export_totaldensity_Num to generate 
%     Figures 6(e), the total population density in numerical solutions.
% This function generates the discrete simulation data with parameters
% P = 0.001, M=1, C(0)=0.16 where B=1 and w1=16, r=1 for movement event, 
% r=4 for birth/death events,
% and numerical solutions of Equation (15) with dt = 1, dx = 1

P=0.001;%Probability of attempting to grow
ini=0.16;%initial density C(0)
realisations=10; %number of realisations, use larger number for higher accuracy
AlleeParameter=0.4; %Allee threshold
r=4;%rings of the spatial template for birth/death event
neighbour=3*r*(r+1);
shape=1;%specify the initial condition is 1-D
distribution=1;%1 for B=1 and 2 for B=0.64;
MaxT=10000;%time step for producing total density 

%Obtain averaged column density in the discrete model
t=200;%time step
density_discrete=Export_columndensity1D_Dis(realisations,t,ini,P,r,neighbour,AlleeParameter,shape,distribution);

%Obtain averaged column density in the continuum model
density_continuum=Export_columndensity1D_Num(P,ini,distribution,t,AlleeParameter);

subplot(2,1,1)
plot(density_continuum(:,1),density_continuum(:,2),'b',density_discrete(:,1),density_discrete(:,2),'r--')
xlabel('x') 
ylabel('C(x,T)')
legend('continuum','discrete')


%produce total density in the discrete model
totaldensity_discrete=Export_totaldensity_Dis(realisations,MaxT,ini,P,r,neighbour,AlleeParameter,shape,distribution);

%produce total density in the continuum model
totaldensity_continuum=Export_totaldensity1D_Num(P,ini,distribution,MaxT,AlleeParameter);

subplot(2,1,2)
plot(totaldensity_continuum(:,1),totaldensity_continuum(:,2),'b',totaldensity_discrete(:,1),totaldensity_discrete(:,2),'r--')
xlabel('T') 
ylabel('C(T)') 
legend('continuum','discrete')

%Discrete model
%Export total density
function result=Export_totaldensity_Dis(realisations,MaxT,ini,Pi,r,neighbour,AlleeParameter,shape,distribution)
    tstep=MaxT/100;%obtain 100 data points
    totaldensity=Produce_totaldensity(MaxT,tstep,ini,Pi,r,neighbour,AlleeParameter,shape,distribution);
    total=zeros(realisations,length(totaldensity));
    total(1,:)=totaldensity;
    parfor i=2:realisations
        total(i,:)=Produce_totaldensity(MaxT,tstep,ini,Pi,r,neighbour,AlleeParameter,shape,distribution);
    end  
    totaldensity_final=sum(total)./realisations;
    a=(0:length(totaldensity_final)-1).*(Pi*MaxT/100);
    result=[a',totaldensity_final'];
end
function totaldensity=Produce_totaldensity(MaxT,tstep,ini,P,r,neighbour,AlleeParameter,shape,distribution)
  [N,M,A]=initial(ini,shape,distribution);
  [count,indexA]=Initialize_agentindex(A,N,M);
  totalnumber=N*M/2;
  totaldensity=ones(1,MaxT/tstep+1);
  totaldensity(1)=count;
  coor=1;
  totaldensity(1)=count;
  for i=1:MaxT
      if count>totalnumber-1
          totaldensity(coor+1:end)=totalnumber;
          break;
      elseif count<1
          totaldensity(coor+1:end)=0;
          break;
      end
      [count,indexA,A]=iteration(A,N,M,indexA,count,P,r,neighbour,AlleeParameter);
      if i>coor*tstep-1
          coor=coor+1;
          totaldensity(coor)=count;
      end
  end
  totaldensity=totaldensity./totalnumber;
end
%Export averaged column density
function result=Export_columndensity1D_Dis(realisations,t,ini,P,r,neighbour,AlleeParameter,shape,distribution)
    totaldensity=Average_columndensity(t,ini,P,r,neighbour,AlleeParameter,shape,distribution);
    total=zeros(realisations,length(totaldensity));
    total(1,:)=totaldensity;
    parfor i=2:realisations
        total(i,:)=Average_columndensity(t,ini,P,r,neighbour,AlleeParameter,shape,distribution);
    end  
    averageddensity=sum(total)./realisations;
    a=1:length(averageddensity);
    result=[a',averageddensity'];
end
function averageddensity=Average_columndensity(t,ini,P,r,neighbour,AlleeParameter,shape,distribution)
  [N,M,A]=initial(ini,shape,distribution);
  [count,indexA]=Initialize_agentindex(A,N,M);
  for i=1:t
      [count,indexA,A]=iteration(A,N,M,indexA,count,P,r,neighbour,AlleeParameter);
  end
  averageddensity=Calculate_columndensity(A,N,M);
end
function density=Calculate_columndensity(A,N,M)
    len=M/2;
    density=zeros(1,len);
    i=1;
    while i<len+1
        density(i)=sum(A(:,2*i-1)+A(:,2*i))/(N);
        i=i+1;
    end
end
function [count,indexA]=Initialize_agentindex(A,N,M)
    count=0;
    indexA=zeros(1,2);
    for i=1:N
        for j=1:M
            if A(i,j)>0
                count=count+1;
                indexA(count,1)=i;
                indexA(count,2)=j;
            end
        end
    end
end
%Initial conditions
function [N,M,A]=initial(ini,shape,distribution)
N=116;%this has to be an even number so that the periodic boundary conditions works for the hexagonal lattice
M=100*2;%this has to be an even number as well
B=0.64;%local density, will be used if 'distribution==2'.
if shape==1    %1-D
    if distribution==1%B=1
        N_up=0;
        N_down=N;
        len=ini*M/4;
        M_right=round(M/4+len)*2;
        M_left=M-M_right;
        A=zeros(N,M);
        for i=N_up+1:N_down
            for j=M_left+1:M_right
                if mod(i+j,2)>0
                    A(i,j)=1;
                end
            end
        end
    else %B!=1
        N_up=0;
        N_down=N;
        len=(M/4)*0.64;
        M_right=ceil((M/4+len))*2;
        M_left=M-M_right;
        count_totalagent=round(B*N*M/2);
        index=zeros(count_totalagent,2);
        count_index=0;
        for i=N_up+1:N_down
            for j=M_left+1:M_right
                if mod(i+j,2)>0
                    count_index=count_index+1;
                    index(count_index,1)=i;
                    index(count_index,2)=j;
                end
            end
        end
        amount=round(ini*count_totalagent);
        randomchoose=randperm(count_totalagent,amount);
        A=zeros(N,M);
        for i=1:amount
            A(index(randomchoose(i),1),index(randomchoose(i),2))=1;
        end
    end
elseif shape==2 %2-D
    if distribution==1%B=1
        len=sqrt(ini)*M/4;
        wid=len*116/100;
        N_down=round(N/2+wid);
        N_up=N-N_down;
        M_right=round((M/4+len)*2);
        M_left=M-M_right;
        A=zeros(N,M);
        for i=N_up+1:N_down
            for j=M_left+1:M_right
                if mod(i+j,2)>0
                    A(i,j)=1;
                end
            end
        end
    else %B!=1
        len=(M/4)*sqrt(0.64);
        wid=len*2/sqrt(3);
        N_down=ceil(N/2+wid);
        N_up=N-N_down;
        M_right=ceil((M/4+len))*2;
        M_left=M-M_right;
        count_totalagent=round(B*N*M/2);
        index=zeros(count_totalagent,2);
        count_index=0;
        for i=N_up+1:N_down
            for j=M_left+1:M_right
                if mod(i+j,2)>0
                    count_index=count_index+1;
                    index(count_index,1)=i;
                    index(count_index,2)=j;
                end
            end
        end
        amount=round(ini*count_totalagent);
        randomchoose=randperm(count_totalagent,amount);
        A=zeros(N,M);
        for i=1:amount
            A(index(randomchoose(i),1),index(randomchoose(i),2))=1;
        end
    end
else%0-D
    N_up=0;
    N_down=N;
    M_left=0;
    M_right=M;
    count_totalagent=N*M/2;
    index=zeros(count_totalagent,2);
    count_index=0;
    for i=N_up+1:N_down
        for j=M_left+1:M_right
            if mod(i+j,2)>0
                count_index=count_index+1;
                index(count_index,1)=i;
                index(count_index,2)=j;
            end
        end
    end
    amount=round(ini*count_totalagent);
    randomchoose=randperm(count_totalagent,amount);
    A=zeros(N,M);
    for i=1:amount
        A(index(randomchoose(i),1),index(randomchoose(i),2))=1;
    end
end
end

%Continuum model
%Export total density
function result=Export_totaldensity1D_Num(P,ini,distribution,T,A)
    tstep=100;
    a=2.5*P;
    L = 100;
    dx = 1;
    D=1/4;
    N=L/dx;
    u0 = zeros(N,1)+1e-10;
    if distribution<2 
        u_initial = 1; 
        u0 = zeros(N,1)+1e-10;
        len=ini*N/2;
        right=round(N/2+len);
        left=N/2-(right-N/2);
        for i=left+1:right
            u0(i)=u_initial;
        end
    else
        u_initial = ini/0.64; 
        len=0.64*N/2;
        right=round(N/2+len);
        left=N/2-(right-N/2);
        for i=left+1:right
            u0(i)=u_initial;
        end
    end
    u0=reshape(u0,[],1);
    tspan = 0:tstep:T;
    [t,u] = ode45(@(t,u) LineApproach_odefun_1D(t,u,N,a,A,D,dx), tspan, u0);
    u = reshape(u, [], N,1);
    result=zeros(1,length(tspan));
    for i=1:length(tspan)
        result(i)=sum(sum(u(i,:)))/(N);
    end
    a=0:length(result)-1;
    result=[(a/10)',result'];
end
%Export average density 
function result=Export_columndensity1D_Num(P,ini,distribution,T,A)%line appoach for solving 1D reaction diffusion equation to obtain total population density   
    a=2.5*P;%lambda*a/Pi
    L = 100;
    dx = 1;
    D=1/4;%D0
    N=L/dx;
    u0 = zeros(N,1)+1e-10;
    if distribution==1 
        u_initial = 1; 
        u0 = zeros(N,1)+1e-10;
        len=ini*N/2;
        right=round(N/2+len);
        left=N-right;
        for i=left+1:right
            u0(i)=u_initial;
        end
    else
        u_initial = ini; 
        len=0.64*N/2;
        right=round(N/2+len);
        left=N-right;
        for i=left+1:right
            u0(i)=u_initial;
        end
    end
    u0=reshape(u0,[],1);
    tspan = 0:T/2:T;
    [t,u] = ode45(@(t,u) LineApproach_odefun_1D(t,u,N,a,A,D,dx), tspan, u0);
    u = reshape(u, [], N,1);
    column=u(end,:);
    x=1:length(column);
    result=[x',column'];
end