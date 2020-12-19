%%

%This program takes the coupled Lorenz ODEs which describe the behavior of
%weather and use them to show how deterministic systems become chaotic or
%appear random in the phenomenon called the butterfly effect. First I solve
%the ODEs using Runge Kutta 2nd order from t=0 to t=50 with nonchaotic paramters 
%and plot it.  Then I change rho to 28 which makes the system chaotic which
%I demonstrate with a small change in initial conditions we have big
%changes in values of y later on. I then finally make graphs of the famous
%Lorenz attractor for y vs. z in 2d and then x y and z in 3d.

clear all
sigma =10;
rho= 13;   %non-chaotic parameter
beta=8/3;
x(1)=0;
y(1)=1;
z(1)=0;
dt= 0.0001;    %extremely small time step in order to minimize numerical error
nSteps= 500000;
t(1)=0;
for n= 1:nSteps    %Runge Kutta 2nd order used
    ymid= y(n) +(rho*x(n) -y(n)-x(n)*z(n))*dt/2;  %calculate all the mid values
    xmid= x(n) +sigma*(y(n)-x(n))*dt/2;
    zmid= z(n) + (x(n)*y(n)-beta*z(n))*dt/2;
    x(n+1)= x(n)+ sigma*(ymid-x(n))*dt;   %Now plug the mids into Eulers
    y(n+1)= y(n) + (rho*xmid -y(n)-xmid*zmid)*dt;  
    z(n+1)= z(n) + (xmid*ymid-beta*z(n))*dt;
    t(n+1)=n*dt;
end
figure(1)
plot(t, y, 'k');  %First figure shows how y changes with respect to t
xlabel('Time');
ylabel('Y for Lorenz Coupled ODEs');
title('Non-Chaotic Plot for Y as a Function of t');
%The Lorenz System exhibits non chaotic behavior so far. It seems
%deterministic and initial conditions for t=0 to t=50
%don't really seem to change the course of the function y(t).

%%
clear all
sigma =10;
rho= 28;   %Chaotic parameter
beta=8/3;
x1(1)=0;
y1(1)=1;
z1(1)=0;
dt= 0.0001;
nSteps= 500000;
t(1)=0;
for n= 1:nSteps   %Solve the ODEs again!
    y1mid= y1(n) +(rho*x1(n) -y1(n)-x1(n)*z1(n))*dt/2;
    x1mid= x1(n) +sigma*(y1(n)-x1(n))*dt/2;
    z1mid= z1(n) + (x1(n)*y1(n)-beta*z1(n))*dt/2;
    x1(n+1)= x1(n)+ sigma*(y1mid-x1(n))*dt;
    y1(n+1)= y1(n) + (rho*x1mid -y1(n)-x1mid*z1mid)*dt;
    z1(n+1)= z1(n) + (x1mid*y1mid-beta*z1(n))*dt;
    t(n+1)=n*dt;
end
figure(2)
subplot(2, 1, 2)
plot(t, y1, 'k')
ylabel('Chaotic Y value')
xlabel('Time')
title('Chaotic Y vs t for Lorenz ODEs')

x(1)=0;
y(1)=1+10^-9;   % very small change in initial condition y
z(1)=0;
dt= 0.0001;
nSteps= 500000;
t(1)=0;
for n= 1:nSteps
    ymid= y(n) +(rho*x(n) -y(n)-x(n)*z(n))*dt/2;
    xmid= x(n) +sigma*(y(n)-x(n))*dt/2;
    zmid= z(n) + (x(n)*y(n)-beta*z(n))*dt/2;
    x(n+1)= x(n)+ sigma*(ymid-x(n))*dt;
    y(n+1)= y(n) + (rho*xmid -y(n)-xmid*zmid)*dt;
    z(n+1)= z(n) + (xmid*ymid-beta*z(n))*dt;
    t(n+1)=n*dt;
end

subplot(2, 1, 1)
% plots y and y1 on subplot with only extremely small change in inital
% condition to see how different they become by t =50
plot(t, y, 'b')
ylabel('Chaotic Y value')
xlabel('Time')
title('Chaotic Y vs t with 10^{-9} Change in IC')

%%
figure(3)
plot(y1, z1); % plots 2D version of Lorenz attractor with y and z
xlabel('Y value')
ylabel('Z value')
title('2D Lorenz Attractor')
figure(4)
scatter3(x1, y1, z1, 1, 'b'); % plots all three x, y and z together on a 3D plot to show full Lorenz attractor
title('3D Lorenz Attractor')
