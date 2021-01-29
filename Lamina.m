%LAMINA SOLVER FUNCTION

function [T_theta,Q_ij,Q_bar]=Lamina(theta,v12,E1,E2,G12)
%Dummy Givens

v21=v12*E2/E1;

m=cosd(theta);
n=sind(theta);

%Transformation Matrix
T_theta=zeros(3);
T_theta(1,1)=m^2;
T_theta(1,2)=n^2;
T_theta(1,3)=2*m*n;
T_theta(2,1)=T_theta(1,2);
T_theta(2,2)=T_theta(1,1);
T_theta(2,3)=-T_theta(1,3);
T_theta(3,1)=-m*n;
T_theta(3,2)=m*n;
T_theta(3,3)=m^2-n^2;

%Stiffness Matrix Lamina Axis
Q_ij=zeros(3);
den= 1-v12*v21;
Q_ij(1,1) = E1/den;
Q_ij(1,2) = v21*E1/den;
Q_ij(2,1) = Q_ij(1,2);
Q_ij(2,2) = E2/den;
Q_ij(3,3) = G12;

%Stiffness Matrix General Axis
u1 = (3*Q_ij(1,1)+3*Q_ij(2,2)+2*Q_ij(1,2)+4*Q_ij(3,3))/8;
u2 = 0.5*(Q_ij(1,1) - Q_ij(2,2));
u3 = (Q_ij(1,1)+Q_ij(2,2)-2*Q_ij(1,2)-4*Q_ij(3,3))/8 ;   
u4 = (Q_ij(1,1)+Q_ij(2,2)+6*Q_ij(1,2)-4*Q_ij(3,3))/8;
u5 = (Q_ij(1,1)+Q_ij(2,2)-2*Q_ij(1,2)+4*Q_ij(3,3))/8;

m1 = cosd(2*theta);     
n1 = sind(2*theta);
m2 = cosd(4*theta); 
n2 = sind(4*theta);

Q_bar=zeros(3);
Q_bar(1,1) = u1+u2*m1+u3*m2;
Q_bar(1,2) = u4-u3*m2;
Q_bar(2,1) = Q_bar(1,2);
Q_bar(2,2) = u1-u2*m1+u3*m2;
Q_bar(1,3) = 0.5*u2*n1+u3*n2;
Q_bar(3,1) = Q_bar(1,3);
Q_bar(2,3) = 0.5*u2*n1-u3*n2;
Q_bar(3,2) = Q_bar(2,3);
Q_bar(3,3) = u5-u3*m2;

