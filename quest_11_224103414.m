%224103414
clc
clear all
syms dia1 len1 dia2 len2 G w Ip z t

%Transfer Matrix Method


J1=pi*dia1^4/32; % (m^4)
J2=pi*dia2^4/32; % (m^4)


kt1=G*J1/len1; % (N/m)
kt2=G*J2/len2; % (N/m)


P=[1 0;-w^2*Ip 1];  % point mat
F1=[1 1/kt1;0 1];    % field mat
F2=[1 1/kt2;0 1];

%  Boundary Condition:
s=[z;t];
A=P*F1*F2*s;
A1=subs(A,{z},{0});
A1=subs(A,{dia1,len1,dia2,len2,G,Ip,z},{0.01,0.6,0.03,0.4,8*10^10,0.01,0});


wnf= solve(A1(2),w);
disp('Transfer Matrix Method')
fprintf('\n')
disp('The Frequency Using Transfer Matrix Method :')
wnf=vpa(wnf(2))




disp('Using FEM')


Mass1=[0 0;0 0];
k1=(G*J2/len2)*[1 -1;-1 1];                         % Element 1
k1=subs(k1,{G,dia2,len2},{8*10^10,0.03,0.4});


k2=(G*J1/len1)*[1 -1;-1 1];
Mass2=[0 0;0 Ip];                        % Element 2
Mass2=subs(Mass2,Ip,0.01);
k2=subs(k2,{G,dia1,len1},{0.8*10^11,0.01,0.6});


Mass=zeros(3,3);
Mass(1:2,1:2)=Mass1;                       % Global Mass matrix
Mass(2:3,2:3)=Mass(2:3,2:3)+Mass2

K=zeros(3,3);
K(1:2,1:2)=k1;              % Global Stiffness matrix
K(2:3,2:3)=K(2:3,2:3)+k2


Mass=vpa(-w^2*Mass);
A2=vpa(Mass+K);                 % Boundary condition


A2=A2(2:3,2:3);              % Reduced Mass and Stiffness mastrix



A2=det(A2);
wnf2=solve(A2,w);
disp('The Frequency using FEM:')            %solving for frequency
wnf2=wnf2(2)










