% Filename: LineApproach_odefun.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Jan 2021
% Reference:  Y. Li, S.T. Johnston, P.R. Buenzli, P. van Heijster, M.J. Simpson (2021) 
% Dimensionality affects extinction of bistable populations.
% The script contains the function 'LineApproach_odefun', which applies
% the method of lines for solving a 2-dimensional RDE

function dudt = LineApproach_odefun(t,u,N,a,A,D,dx)
u=reshape(u,N,N);
dudt1=zeros(N,N);

for i=2:N-1
    for j=2:N-1
        dudt1(i,j)=D*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))/(dx^2)...
                   +a*u(i,j)*(1-u(i,j))*(u(i,j)-A);
    end
end

i=1;
for j=2:N-1
    dudt1(i,j)=D*(u(i+1,j)+u(N,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))/(dx^2)...
                   +a*u(i,j)*(1-u(i,j))*(u(i,j)-A);
end

i=N;
for j=2:N-1
    dudt1(i,j)=D*(u(1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))/(dx^2)...
                   +a*u(i,j)*(1-u(i,j))*(u(i,j)-A);
end

j=1;
for i=2:N-1
    dudt1(i,j)=D*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,N)-4*u(i,j))/(dx^2)...
                   +a*u(i,j)*(1-u(i,j))*(u(i,j)-A);
end

j=N;
for i=2:N-1
    dudt1(i,j)=D*(u(i+1,j)+u(i-1,j)+u(i,1)+u(i,j-1)-4*u(i,j))/(dx^2)...
                   +a*u(i,j)*(1-u(i,j))*(u(i,j)-A);
end

i=1;
j=1;
    dudt1(i,j)=D*(u(i+1,j)+u(N,j)+u(i,j+1)+u(i,N)-4*u(i,j))/(dx^2)...
                   +a*u(i,j)*(1-u(i,j))*(u(i,j)-A);
i=N;
j=1;
    dudt1(i,j)=D*(u(1,j)+u(i-1,j)+u(i,j+1)+u(i,N)-4*u(i,j))/(dx^2)...
                   +a*u(i,j)*(1-u(i,j))*(u(i,j)-A);
i=1;
j=N;
    dudt1(i,j)=D*(u(i+1,j)+u(N,j)+u(i,1)+u(i,j-1)-4*u(i,j))/(dx^2)...
                   +a*u(i,j)*(1-u(i,j))*(u(i,j)-A);
i=N;
j=N;
    dudt1(i,j)=D*(u(1,j)+u(i-1,j)+u(i,1)+u(i,j-1)-4*u(i,j))/(dx^2)...
                   +a*u(i,j)*(1-u(i,j))*(u(i,j)-A);
dudt=reshape(dudt1,[],1);
