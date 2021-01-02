% Filename: PhaseDiagram_Discrete.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Jan 2021
% Reference:  Y. Li, S.T. Johnston, P.R. Buenzli, P. van Heijster, M.J. Simpson (2021) 
% Dimensionality affects extinction of bistable populations.
% The script contains:
%   - one call to the function 'drawphaseplane' to generate the phase
%   diagram of discrete simulations with the 0-dimensional initial condition. 
%   While this function can also generate other phase diagrams via changing 
%   the parameters 'shape' and 'distribution'.

%Produce phase diagrams with P/M and C_0 in the discrete model
realisations=10;
AlleeParameter=0.4;
r=4;
neighbour=3*r*(r+1);
shape=0;%0 for 0-D; 1 for 1-D initial condition; 2 for 2-D;
distribution=1;%1 for B=1 and 2 for B=0.64;
%In summary, shape=0, distribution=1 refers to Figure 8 (0D Case), 
%shape=1,distribution=1 refers to Figure 9 (1D with B=1), 
%shape=1,distribution=2 refers to Figure 10 (2D with B=1),
%shape=1,distribution=2 refers to Figure S1 (1D with B!=1) in the supplementary materials
%and shape=2,distribution=2 refers to Figure S2 (2D with B!=1) in the supplementary materials

drawphaseplane(realisations,r,neighbour,AlleeParameter,shape,distribution)

%produce phase diagrams in the discrete model. Please note that the practical 
%parameter regimes are time-consuming. Parallel compution with 20 
%cpus will need around 12 hours to obtain the result for one type of
%initial condition. To shorten the operation time, here we first provide a toy
%regime.
function drawphaseplane(realisations,r,neighbour,AlleeParameter,shape,distribution)
    if distribution==1
        %a toy regime
        Pi=[0.01,0.02];
        ini=[0.3,0.4];
        %the practical regime
        %Pi=0.001:0.001:0.04;
        %ini=0.1:0.01:0.6;
    else
        %a toy regime
        Pi=[0.01,0.02];
        ini=[0.3,0.4];
        %the practical regime
        %Pi=0.001:0.001:0.04;
        %ini=0.3:0.02:1;
    end
    n=length(Pi);
    m=length(ini);
    result=zeros(n*m,3);
    count=1;
    for i=1:n
        for j=1:m
            MaxT=max(30/Pi(i),10000);%iterate a sufficient long time
            if distribution==1
                result(count,1)=ini(j);
            else
                result(count,1)=ini(j)*0.64;
            end
            result(count,2)=Pi(i);
            result(count,3)=Export_survivalrate(realisations,MaxT,ini(j),Pi(i),r,neighbour,AlleeParameter,shape,distribution);
            count=count+1;
        end
    end
    sz=150;
    colorMap = [linspace(1,0,256)',linspace(1,0,256)',ones(256,1)];
    colormap(colorMap);
    scatter(result(:,1),result(:,2),sz,result(:,3),'filled');
    xlabel('C_0') 
    ylabel('P/M') 
end
%Export survival rate.
function surviverate=Export_survivalrate(realisations,MaxT,ini,Pi,r,neighbour,AlleeParameter,shape,distribution)
    total=zeros(1,realisations);
    parfor i=1:realisations
        total(i)=Calculate_endstate(MaxT,ini,Pi,r,neighbour,AlleeParameter,shape,distribution);
    end  
    surviverate=sum(total)./realisations;
end
function density=Calculate_endstate(t,ini,Pi,r,neighbour,AlleeParameter,shape,distribution)
  [N,M,A]=initial(ini,shape,distribution);
  [count,indexA]=Initialize_agentindex(A,N,M);
  totalnumber=N*M/2;
  threshold=0.7*totalnumber;%stop when the total density is higher than 0.7, this is for speeding up the simulation
  for i=1:t
      if count>threshold||count<1
          break
      end
      [count,indexA,A]=iteration(A,N,M,indexA,count,Pi,r,neighbour,AlleeParameter);
  end
  if count/totalnumber>AlleeParameter
      density=1;
  else
      density=0;
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
    else %B=0.64
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
        wid=len*2/sqrt(3);
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
    else %B=0.64
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