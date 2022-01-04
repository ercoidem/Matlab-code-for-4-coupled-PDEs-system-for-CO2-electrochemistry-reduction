%See https://www.mathworks.com/help/matlab/math/solve-system-of-pdes.html
%for official guidelines for solving PDEs

function pdex4
m = 0; %Refer to https://www.mathworks.com/help/matlab/math/partial-differential-equations.html for definitions of f,m,c,s,etc. here and following codes
x = 0:0.00001:0.0001; %set range of x
t = 0:1:20;%set range of t

%Set up pdepe solver
sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t); 
u1 = sol(:,:,1);
u2 = sol(:,:,2);
u3 = sol(:,:,3);
u4 = sol(:,:,4);

%Plotting
figure
surf(fliplr(x),t,u1)
title('CO2(x,t)')
xlabel('Distance x')
ylabel('Time t')

figure
surf(fliplr(x),t,u2)
title('HCO3-(x,t)')
xlabel('Distance x')
ylabel('Time t')

figure
surf(fliplr(x),t,u3)
title('CO32-(x,t)')
xlabel('Distance x')
ylabel('Time t')

figure
surf(fliplr(x),t,u4)
title('OH-(x,t)')
xlabel('Distance x')
ylabel('Time t')
end
% --------------------------------------------------------------
%Setting up PDEs' coefficients
function [c,f,s] = pdex4pde(x,t,u,DuDx)
DCO2 = 1.91e-9;
DHCO3 = 9.23e-10;
DCO3 = 1.19e-9;
DOH = 5.27e-9;
k1f = 5.93e3;
k1r = 1.34e-4;
k2f = 1e8;
k2r = 2.15e4;

c=[1;1;1;1];    
s=[-u(1)*u(4)*k1f+u(2)*k1r;u(1)*u(4)*k1f-u(2)*k1r-u(2)*u(4)*k2f+u(3)*k2r;
    u(2)*u(4)*k2f-u(3)*k2r;-u(1)*u(4)*k1f+u(2)*k1r-u(2)*u(4)*k2f+u(3)*k2r];
f=[DCO2;DHCO3;DCO3;DOH].*DuDx;   
end
% --------------------------------------------------------------
%Initial Conditions
function u0 = pdex4ic(x);
u0=[0.0342;0.499;7.6e-4;3.3e-7]; 
end
% --------------------------------------------------------------
%Boundary Conditions
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
pl = [(-1.45e-7);0;0;(1.45e-7)];                               
ql = [1;1;1;1];                                  
pr = [ur(1)-0.0342;ur(2)-0.499;ur(3)-(7.6e-4);ur(4)-(3.3e-7)];                            
qr = [0;0;0;0];          
end