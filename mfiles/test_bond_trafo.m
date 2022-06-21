
% test of Bond transformation VTI--> TTI

C11=34.3e9; C33=22.7e9; C55=5.4e9; C13=10.7e9;

CVTI=[C11, C13,0; C13, C33, 0; 0, 0, C55];

theta=20.0;

t=theta*pi/180.0;

l1=cos(t); 
l2=sin(t); 
l12=l1*l1; 
l22=l2*l2; 
l14=l12*l12; 
l24=l22*l22; 
l13=l1*l12; 
l23=l2*l22;

D=[l12 l22 2.0*l1*l2;
    l22 l12 -2.0*l1*l2;
    -l1*l2 l1*l2 (l12-l22)];


CTTI=D*CVTI*D';

