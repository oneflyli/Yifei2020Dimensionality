% Filename: Code_Figure5.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Jan 2021
% Reference:  Y. Li, S.T. Johnston, P.R. Buenzli, P. van Heijster, M.J. Simpson (2021) 
% Dimensionality affects extinction of bistable populations.
% The script contains:
%   - one call to the function Export_columndensity2D_Dis to generate 
%     Figures 5(d), the averaged density along y=50 at time t in discrete simulations.
%   - one call to the function Export_columndensity2D_Num to generate 
%     Figures 5(d), the density along y=50 at t in numerical solutions.
%   - one call to the function Export_totaldensity_Dis to generate 
%     Figures 5(e), the total population density in discrete simulations.
%   - one call to the function Export_totaldensity2D_Num to generate 
%     Figures 5(e), the total population density in numerical solutions.
% This function generates the discrete simulation data with parameters
% P = 0.001, M=1, C(0)=0.16 where B=1 and w=40, r=1 for movement event, 
% r=4 for birth/death events,
% and numerical solutions of Equation (9) with dt = 1, dx = 1

P=0.001;%Probability of attempting to grow
ini=0.16;%initial density C(0)
realisations=10; %number of realisations, use larger number for higher accuracy
AlleeParameter=0.4; %Allee threshold
r=4;%rings of the spatial template for birth/death event
neighbour=3*r*(r+1);
shape=2;%specify the initial condition is 2-D
distribution=1;%1 for B=1 and 2 for B=0.64;
MaxT=10000;%time step for producing total density
t=500;%time step for producing density profile at y=50


%Obtain averaged density along y=50 in the discrete model
density_discrete=Export_columndensity2D(t,realisations,P,AlleeParameter,r,neighbour,ini,shape,distribution);

%Obtain averaged density along y=50 in the continuum model
totaldensity_continuum=Export_columndensity2D_Dis(P,ini,distribution,MaxT,AlleeParameter);

%Note that the number of realisation has to be very big in order to match
%the outcomes of the discrete and continuum models.

%produce total density in the discrete model
totaldensity_discrete=Export_totaldensity_Dis(realisations,MaxT,ini,P,r,neighbour,AlleeParameter,shape,distribution);

%produce total density in the continuum model
density_continuum=Export_totaldensity2D_Num(P,ini,distribution,t,AlleeParameter);

subplot(2,1,1)
plot(density_continuum(:,1),density_continuum(:,2),'b',density_discrete(:,1),density_discrete(:,2),'r--')
xlabel('x') 
ylabel('C(x,50,T)')
legend('continuum','discrete')

subplot(2,1,2)
plot(totaldensity_continuum(:,1),totaldensity_continuum(:,2),'b',totaldensity_discrete(:,1),totaldensity_discrete(:,2),'r--')
xlabel('T') 
ylabel('C(T)') 
legend('continuum','discrete')

%Discrete model
%Export total density
function result=Export_totaldensity_Dis(realisations,MaxT,ini,P,r,neighbour,AlleeParameter,shape,distribution)
    tstep=MaxT/100;%obtain 100 data points
    totaldensity=Produce_totaldensity(MaxT,tstep,ini,P,r,neighbour,AlleeParameter,shape,distribution);
    total=zeros(realisations,length(totaldensity));
    total(1,:)=totaldensity;
    parfor i=2:realisations
        total(i,:)=Produce_totaldensity(MaxT,tstep,ini,P,r,neighbour,AlleeParameter,shape,distribution);
    end  
    totaldensity_final=sum(total)./realisations;
    a=(0:length(totaldensity_final)-1).*(P*MaxT/100);
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
%Export averaged density along y=50
function getcutcontour=Export_columndensity2D(t1,realisations,P,AlleeParameter,r,neighbour,ini,shape,distribution)
    A_1=Export_pattern(t1,realisations,P,AlleeParameter,r,neighbour,ini,shape,distribution);
    n=length(A_1(:,1));%row
    m=length(A_1(1,:));%column
    column=zeros(1,2);
    row=round(n/2);
    for i=1:m/2
        column(i,1)=i;
        column(i,2)=((A_1(row,2*i))+(A_1(row,2*i-1)));
    end
    getcutcontour=column;
end
function [A]=Export_pattern(t1,realisations,P,AlleeParameter,r,neighbour,ini,shape,distribution)
    A_1=Calculate_pattern(t1,ini,P,r,neighbour,AlleeParameter,shape,distribution);
    A_2=zeros(length(A_1(:,1)),length(A_1(1,:)),realisations);
    A_2(:,:,1)=A_1;
    parfor i=2:realisations
        A_2(:,:,i)=Calculate_pattern(t1,ini,P,r,neighbour,AlleeParameter,shape,distribution);
    end
    for i=2:realisations
        A_1=A_1+A_2(:,:,i);
    end
    A=A_1/realisations;
end
function [A]=Calculate_pattern(t1,ini,P,r,neighbour,AlleeParameter,shape,distribution)
  [N,M,A]=initial(ini,shape,distribution);
  [count,indexA]=Initialize_agentindex(A,N,M);
  totalnumber=(N*M)/2;
  for i=1:t1
      if count>totalnumber-1||count<1
          break;
      end
      [count,indexA,A]=iteration(A,N,M,indexA,count,P,r,neighbour,AlleeParameter);
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
if shape==1    %1-D initial condition
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
elseif shape==2 %2-D initial condition
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
function total=Export_columndensity2D_Dis(P,ini,distribution,T,A)   
    tstep=100;
    a=2.5*P;
    L = 100;
    dx = 1;
    D=1/4;
    N=L/dx;
    u0 = zeros(N,N)+1e-10;
    if distribution==1
        u_initial = 1;
        len=sqrt(ini)*N/2;
        right=round(N/2+len);
        left=N/2-(right-N/2);
        for i=left+1:right
            for j=left+1:right
                u0(i,j)=u_initial;
            end
        end
    else 
        u_initial = ini/0.64;
        len=0.8*N/2;
        right=round(N/2+len);
        left=N/2-(right-N/2);
        for i=left+1:right
            for j=left+1:right
                u0(i,j)=u_initial;
            end
        end
    end
    u0=reshape(u0,[],1);
    tspan = 0:tstep:T;
    [t,u] = ode45(@(t,u) LineApproach_odefun(t,u,N,a,A,D,dx), tspan, u0);
    u = reshape(u, [], N,N);
    total=zeros(1,length(tspan));
    for i=1:length(tspan)
        total(i)=sum(sum(u(i,:,:)))/(N*N);
    end
    a=0:length(total)-1;
    total=[(a/10)',total'];
    %plot(total)
end   
%Export average density along y=50
function last=Export_totaldensity2D_Num(P,ini,distribution,T,A)
    a=2.5*P;
    L = 100;
    dx = 1;
    D=1/4;%D0
    N=L/dx;
    u0 = zeros(N,N)+1e-10;
    if distribution==1
        u_initial = 1;
        len=sqrt(ini)*N/2;
        right=round(N/2+len);
        left=N-right;
        for i=left+1:right
            for j=left+1:right
                u0(i,j)=u_initial;
            end
        end
    else 
        u_initial = ini;
        len=0.8*N/2;
        right=round(N/2+len);
        left=N/2-(right-N/2);
        for i=left+1:right
            for j=left+1:right
                u0(i,j)=u_initial;
            end
        end
    end
    u0=reshape(u0,[],1);
    tspan = 0:T/2:T;
    [t,u] = ode45(@(t,u) LineApproach_odefun(t,u,N,a,A,D,dx), tspan, u0);
    u = reshape(u, [], N,N);
    indexu=u(3,:,:);
    indexu=reshape(indexu,N,N);
    count=1;
    patternindex=zeros(1,3);
    for i=1:length(indexu(:,1))
        for j=1:length(indexu(1,:))
            patternindex(count,1)=i;
            patternindex(count,2)=j;
            patternindex(count,3)=indexu(i,j);
            count=count+1;
        end
    end
    len=length(patternindex(:,1));
    last=zeros(1,2);
    count=1;
    for i=1:len
        if patternindex(i,2)==round(L/2)
            last(count,1)=patternindex(i,1);
            last(count,2)=patternindex(i,3);
            count=count+1;
        end
    end
end