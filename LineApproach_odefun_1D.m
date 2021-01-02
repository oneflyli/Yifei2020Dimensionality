% Filename: LineApproach_odefun.m
% Author: Yifei Li
% Queensland University of Technology, Brisbane, Australia, Jan 2021
% Reference:  Y. Li, S.T. Johnston, P.R. Buenzli, P. van Heijster, M.J. Simpson (2021) 
% Dimensionality affects extinction of bistable populations.
% The script contains the function 'LineApproach_odefun', which applies
% the method of lines for solving a 1-dimensional RDE

function dudt = LineApproach_odefun_1D(t,u,N,a,A,D,dx)
u=reshape(u,N,1);
dudt1=zeros(N,1);
for i=2:N-1
        %logistic
        %dudt1(i)=D*(u(i+1)+u(i-1)-2*u(i))/(dx^2)...
        %           +a*u(i)*(1-u(i));
        %Allee
        dudt1(i)=D*(u(i+1)+u(i-1)-2*u(i))/(dx^2)...
                   +a*u(i)*(1-u(i))*(u(i)-A);
end

i=1;
%logistic
%dudt1(i)=D*(u(i+1)-u(i))/(dx^2)...
%                   +a*u(i)*(1-u(i));
%Allee
dudt1(i)=D*(u(i+1)+u(N)-2*u(i))/(dx^2)...
                   +a*u(i)*(1-u(i))*(u(i)-A);
i=N;
%logistic
%dudt1(i)=D*(u(i-1)-u(i))/(dx^2)...
%                   +a*u(i)*(1-u(i));
%Allee
dudt1(i)=D*(u(1)+u(i-1)-2*u(i))/(dx^2)...
                   +a*u(i)*(1-u(i))*(u(i)-A);
dudt=dudt1;
end