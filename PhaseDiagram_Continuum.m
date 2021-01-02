% Filename: PhaseDiagram_Continuum.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Jan 2021
% Reference:  Y. Li, S.T. Johnston, P.R. Buenzli, P. van Heijster, M.J. Simpson (2021) 
% Dimensionality affects extinction of bistable populations.
% The script contains:
%   - one call to the function 'Obtain_toy' to generate the phase
%   diagram of numerical solutions with the 1-dimensional initial condition. 
%   While this function can also generate other phase diagrams via changing 
%   the parameters 'shape' and 'distribution'.

%Produce phase diagrams with P/M and C_0 in the continuum model

type=1;%type1:1D, B=1; type 2: 1D, B!=1; type 3: 2D, B=1; type 4: 2D, B!=1
AlleeParameter=0.4;
result=Obtain_toy(type,AlleeParameter);
sz=150;
colorMap = [linspace(0.95,0,256)',linspace(0.95,0,256)',ones(256,1)];
colormap(colorMap);
scatter(result(:,1),result(:,2),sz,result(:,3),'filled');
xlabel('C_0') 
ylabel('P/M') 

%As producing one phase diagram is time-consuming, especially for solving 
%the 2-D RDE, this code uses the function 'Obtain_toy' to draw a phase diagram 
%with less data. You can replace it with the function 'Obtain' to 
%obtain the data related to Figures 9--11 and Figures S1, S2.

function result=Obtain_toy(type,A)
    if type==1
        a=0.005:0.005:0.04;%P/M
        ini=0.1:0.1:0.6;%C(0)
        distribution=1;
        n=length(a);
        m=length(ini);
        survive=zeros(n,m);
        parfor i=1:n
            a=0.005:0.005:0.04;
            ini=0.1:0.1:0.6;
            for j=1:m
                survive(i,j)=issurvive1D(a(i),ini(j),distribution,A);
            end
        end
        a=0.005:0.005:0.04;
        ini=0.1:0.1:0.6;
    elseif type==2
        a=0.005:0.005:0.04;
        ini=0.3:0.1:1;
        distribution=2;
        n=length(a);
        m=length(ini);
        survive=zeros(n,m);
        parfor i=1:n
            a=0.005:0.005:0.04;
            ini=0.3:0.1:1;
            for j=1:m
                survive(i,j)=issurvive1D(a(i),ini(j),distribution,A);
            end
        end
        a=0.005:0.005:0.04;
        ini=0.3:0.1:1;
        ini=ini.*0.64;
    elseif type==3
        a=0.005:0.005:0.04;
        ini=0.1:0.1:0.6;
        distribution=1;
        n=length(a);
        m=length(ini);
        survive=zeros(n,m);
        parfor i=1:n
            a=0.005:0.005:0.04;
            ini=0.1:0.1:0.6;
            for j=1:m
                survive(i,j)=issurvive2D(a(i),ini(j),distribution,A);
            end
        end
        a=0.005:0.005:0.04;
        ini=0.1:0.1:0.6;
    else
        a=0.005:0.005:0.04;
        ini=0.3:0.1:1;
        distribution=2;
        n=length(a);
        m=length(ini);
        survive=zeros(n,m);
        parfor i=1:n
            a=0.005:0.005:0.04;
            ini=0.3:0.1:1;
            for j=1:m
                survive(i,j)=issurvive2D(a(i),ini(j),distribution,A);
            end
        end
        a=0.005:0.005:0.04;
        ini=0.3:0.1:1;
        ini=ini.*0.64;
    end
    result=zeros(n*m,3);
    count=1;
    for i=1:n
        for j=1:m
            result(count,1)=ini(j);
            result(count,2)=a(i);
            result(count,3)=survive(i,j);
            count=count+1;
        end
    end
end
function result=Obtain(type,A)
    if type==1
        a=0.001:0.002:0.1;
        ini=0.1:0.01:0.6;
        distribution=1;
        n=length(a);
        m=length(ini);
        survive=zeros(n,m);
        parfor i=1:n
            a=0.001:0.002:0.1;
            ini=0.1:0.01:0.6;
            for j=1:m
                survive(i,j)=issurvive1D(a(i),ini(j),distribution,A);
            end
        end
        a=0.001:0.002:0.1;
        ini=0.1:0.01:0.6;
    elseif type==2
        a=0.001:0.002:0.08;
        ini=0.3:0.02:1;
        distribution=2;
        n=length(a);
        m=length(ini);
        survive=zeros(n,m);
        parfor i=1:n
            a=0.001:0.002:0.08;
            ini=0.3:0.02:1;
            for j=1:m
                survive(i,j)=issurvive1D(a(i),ini(j),distribution,A);
            end
        end
        a=0.001:0.002:0.08;
        ini=0.3:0.02:1;
        ini=ini.*0.64;
    elseif type==3
        a=0.001:0.002:0.1;
        ini=0.1:0.01:0.6;
        distribution=1;
        n=length(a);
        m=length(ini);
        survive=zeros(n,m);
        parfor i=1:n
            a=0.001:0.002:0.1;
            ini=0.1:0.01:0.6;
            for j=1:m
                survive(i,j)=issurvive2D(a(i),ini(j),distribution,A);
            end
        end
        a=0.001:0.002:0.1;
        ini=0.1:0.01:0.6;
    else
        a=0.001:0.002:0.08;
        ini=0.3:0.02:1;
        distribution=2;
        n=length(a);
        m=length(ini);
        survive=zeros(n,m);
        parfor i=1:n
            a=0.001:0.002:0.08;
            ini=0.3:0.02:1;
            for j=1:m
                survive(i,j)=issurvive2D(a(i),ini(j),distribution,A);
            end
        end
        a=0.001:0.002:0.08;
        ini=0.3:0.02:1;
        ini=ini.*0.64;
    end
    result=zeros(n*m,3);
    count=1;
    for i=1:n
        for j=1:m
            result(count,1)=ini(j);
            result(count,2)=a(i);
            result(count,3)=survive(i,j);
            count=count+1;
        end
    end
end
function total=issurvive2D(Pi,ini,distribution,A)%line appoach for solving 1D reaction diffusion equation in parameter space
    T=max(50/Pi,10000);
    a=2.5*Pi;
    L = 100;
    dx = 1;
    N=L/dx;
    D=1/4;
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
    total=sum(sum(u(3,:,:)))/(N*N);
    if total>A
        total=1;
    else
        total=0;
    end
end
function total=issurvive1D(Pi,ini,distribution,A)%line appoach for solving 1D reaction diffusion equation %line appoach for solving 1D reaction diffusion equation    
    T=max(50/Pi,10000);
    a=2.5*Pi;
    L = 100;
    dx = 1;
    N=L/dx; 
    D=1/4;%D0
    u0 = zeros(N,1)+1e-10;
    if distribution==1 
        u_initial = 1; 
        u0 = zeros(N,1)+1e-10;
        len=ini*N/2;
        right=round(N/2+len);
        left=N/2-(right-N/2);
        for i=left+1:right
            u0(i)=u_initial;
        end
    else
        u_initial = ini; 
        len=0.64*N/2;
        right=round(N/2+len);
        left=N/2-(right-N/2);
        for i=left+1:right
            u0(i)=u_initial;
        end
    end
    u0=reshape(u0,[],1);
    tspan = 0:T/2:T;
    [t,u] = ode45(@(t,u) LineApproach_odefun_1D(t,u,N,a,A,D,dx), tspan, u0);
    u = reshape(u, [], N,1);
    total=sum(u(3,:))/(N);
    if total>A
        total=1;
    else
        total=0;
    end
end