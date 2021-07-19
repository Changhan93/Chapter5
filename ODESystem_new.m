% This code describes Model II with nutrient preference. 

% Written by Changhan He in 2021.

function dy=ODESystem(t,y)

global kr0 kr1 dr 
global kg J 
global c1 c2 c 
global a   
global gama q1 q2 r1 r2 K 

AarC=y(1);
N=y(2);
P=y(3);

gA = 1/(AarC/J+1)+c2 ;
fA = kr0+kr1*(AarC^2)/(AarC^2+c1);
W1 = 0.5*(sign(P-q1)+1);
W2 = 0.5*(sign(q1-P)+1);
G1 = W1*(r1*(1-q1/P));
G2 = W2*(r2*(1-q2/P));
GR = kg* gA * ( G1 + G2 );
F1 = 1/(K*GR+1);
F2 = (a*GR^0.5+1)/(c*GR^2+1);

dN=GR*N;
dAarC=fA*( W1*F1 + W2*F2 )-(dr+GR)*AarC; %AarC
dP=-gama*N*(P-q2)^1.2;

dy=[dAarC; dN; dP];
 
end


