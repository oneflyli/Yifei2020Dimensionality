% Filename: Code_Figure7.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Jan 2021
% Reference:  Y. Li, S.T. Johnston, P.R. Buenzli, P. van Heijster, M.J. Simpson (2021) 
% Dimensionality affects extinction of bistable populations.
% The script contains:
%   - one call to the function Export_indexdata to generate the evolution
%     pattern with the 0-dimensional initial condition in Figures 7.
%   - one call to the function Export_totaldensity_Dis to generate 
%     Figures 7(d), the total population density in discrete simulations.
%   - one call to the function Export_totaldensity0D_Num to generate 
%     Figures 7(d), the total population density in numerical solutions.
% This function generates the discrete simulation data with parameters
% P = 0.001, M=1, C(0)=0.16 where B=0.16, 
% r=1 for movement event, r=4 for birth/death event, 
% and numerical solutions of Equation (16) with dt = 1, dx = 1.

P=0.001;%Probability of attempting to grow
ini=0.16;%initial density C(0)
realisations=10; %number of realisations, use larger number for higher accuracy
AlleeParameter=0.4; %Allee threshold
r=4;%rings of the spatial template
neighbour=3*r*(r+1);
shape=0;%specify the initial condition is 0-D
distribution=1;%1 for B=1 and 2 for B=0.64;
MaxT=10000;%time step for producing total density 

%produce snapshots of the pattern from one realisation at time t
t1=0;
index1=Export_indexdata(t1,ini,P,r,neighbour,AlleeParameter,shape,distribution);
t2=1000;
index2=Export_indexdata(t2,ini,P,r,neighbour,AlleeParameter,shape,distribution);
t3=2000;
index3=Export_indexdata(t3,ini,P,r,neighbour,AlleeParameter,shape,distribution);

figure
subplot(1,3,1)
sz = 5;
scatter(index1(:,1),index1(:,2),sz,'filled');
axis([0 100 0 100])
axis square
xlabel('x') 
ylabel('y') 
subplot(1,3,2)
scatter(index2(:,1),index2(:,2),sz,'filled');
axis([0 100 0 100])
axis square
xlabel('x') 
ylabel('y') 
subplot(1,3,3)
scatter(index3(:,1),index3(:,2),sz,'filled');
axis([0 100 0 100])
axis square
xlabel('x') 
ylabel('y') 

%Note that the function 'Export_indexdata' can also produce snapshots with
%other inital conditions by changing the parameter 'shape'

%produce total density in the discrete model
totaldensity_discrete=Export_totaldensity_Dis(realisations,MaxT,ini,P,r,neighbour,AlleeParameter,shape,distribution);

%produce total density in the continuum model
totaldensity_continuum=Export_totaldensity0D_Num(ini,MaxT*P,AlleeParameter);


figure
plot(totaldensity_continuum(:,1),totaldensity_continuum(:,2),'b',totaldensity_discrete(:,1),totaldensity_discrete(:,2),'r--')
xlabel('T') 
ylabel('C(T)') 
legend('continuum','discrete')

%Discrete model
%Export total density
function result=Export_totaldensity_Dis(realisations,MaxT,ini,Pi,r,neighbour,AlleeParameter,shape,distribution)
    tstep=MaxT/100;
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
function totaldensity=Produce_totaldensity(MaxT,tstep,ini,Pi,r,neighbour,AlleeParameter,shape,distribution)
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
      [count,indexA,A]=iteration(A,N,M,indexA,count,Pi,r,neighbour,AlleeParameter);
      if i>coor*tstep-1
          coor=coor+1;
          totaldensity(coor)=count;
      end
  end
  totaldensity=totaldensity./totalnumber;
end
%%Export an index indicating the position of all agents
function [index]=Export_indexdata(t,ini,Pi,r,neighbour,AlleeParameter,shape,distribution)
  [N,M,A]=initial(ini,shape,distribution);
  [count,indexA]=Initialize_agentindex(A,N,M);
  totalnumber=(N*M)/2;
  i=0;
  while i<t
      if count>totalnumber-1||count<1
          break;
      end
      [count,indexA,A]=iteration(A,N,M,indexA,count,Pi,r,neighbour,AlleeParameter);
      i=i+1;
  end
  index=[indexA(:,2)/2,indexA(:,1).*100./N];
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
function result=Export_totaldensity0D_Num(ini,T,A)
tspan=[0 T];
u0=ini;
[t,u] = ode45(@(t,u) 2.5*u*(1-u)*(u-A), tspan, u0);
result=[t,u];
end