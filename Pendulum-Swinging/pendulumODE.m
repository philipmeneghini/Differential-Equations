%%
%This code takes the differential eqution that models a pendulum swinging
%with no damping and since it is a differential equation of order two I
%break it into a system of two differential equations. I then solve both
%differential equations by first using Euler's method and then Runge-Kutta
%of order two and three. This is then used to compare the numerical
%difference between the three methods and see how the error changes.

clear all
g = 9.81; %in meters per second squared
l = input('What length do you want the pendulum? (in meters)');
%user inputs the length of the pendulum
theta(1) = input('What theta do you want the pendulum to be pushed back to? (in radians)');
%user inputs initial theta value
w(1) = 0; %this is the w value array that represents the change in theta
dt=0.001; %very small dt to minimize numerical error
nSteps = 50000; %number of dt time steps, equivalent to 50 seconds of the pendulum moving
for k = 1:nSteps
    t(k)=dt*(k-1); %sets up all times in time t array
end

for n = 1:nSteps-1
    theta(n+1)=theta(n)+ w(n)*dt; %Euler's Method of numerically approximating the ODE
    w(n+1)= w(n) -(g/l)*sin(theta(n))*dt;
end

figure(1)
plot(t, theta, 'k');
xlabel('Time in Seconds');
ylabel('Angle of Pendulum in Radians');
title('Numerical Approximation of Pendulum with Eulers');

%Runge-Kutta method of order 2
thetar(1)=theta(1); %initialize new arrays for next numerical method
wr(1)=0;
for n = 1:nSteps-1
    thetamid = thetar(n) + wr(n)*(dt/2); %We calculate half steps to reduce numerical error
    wmid = wr(n) -(g/l)*sin(thetar(n))*(dt/2);
    thetar(n+1)=thetar(n)+ wmid*dt;
    wr(n+1)= wr(n) -(g/l)*sin(thetamid)*dt;
end

figure(2) %plot to show numerical aproximation with Runge-Kutta
plot(t, thetar, 'k');
xlabel('Time in Seconds');
ylabel('Angle of Pendulum in Radians');
title('Numerical Approximation of Pendulum with Runge-Kutta 2nd order');

%Finally Runge Kutta 3rd order
fw = @(t,thet, wi) -(g/l)*sin(thet);
ftheta = @(t, thet, wi) wi;
theta3(1)=theta(1); %initialize new arrays for next numerical method
w3(1)=0;

for n = 1:nSteps-1
     f1 = dt*ftheta(t(n), theta3(n), w3(n));
     k1 = dt*fw(t(n), theta3(n), w3(n));
     f2 = dt*ftheta(t(n)+ (dt/2), theta3(n) + (f1/2), w3(n) + (k1/2));
     k2 = dt*fw(t(n)+ (dt/2), theta3(n) + (f1/2), w3(n) + (k1/2));
     f3 = dt*ftheta(t(n)+dt,  theta3(n)-f1 + 2*f2, w3(n)- k1 + 2*k2);
     k3 = dt*fw(t(n)+dt,  theta3(n)-f1 + 2*f2, w3(n)- k1 + 2*k2);
     w3(n+1) = w3(n) + (1/6)*(k1 +4*k2+k3);
     theta3(n+1)=theta3(n)+ (1/6)*(f1 + 4*f2+f3);
end

figure(3) %plot to show numerical aproximation with Runge-Kutta 4th order
plot(t, theta3, 'k');
xlabel('Time in Seconds');
ylabel('Angle of Pendulum in Radians');
title('Numerical Approximation of Pendulum with Runge-Kutta 3rd order');

errore = abs(theta3-theta); %Error terms for Euler and 2nd Order when looking at 3rd order
error2 = abs(theta3-thetar);

figure(4)
subplot(2, 1, 1);
plot(t, errore, 'k');
ylabel('Error in Euler')
xlabel('Time in Seconds')
title('Difference between Euler and Runge Kutta 3rd order')
subplot(2, 1, 2)
plot(t, error2, 'k');
ylabel('Error in Runge Kutta 2nd Order')
xlabel('Time in Seconds')
title('Difference between Runge Kutta 2nd and 3rd order')

%%
% Initially going from Euler's method to Runge Kutta second order there 
% is a very considerable accuracy jump which becomes transparent 
% when looking at the errors. Euler's method had an accuracy difference of 
% 0.7381 radians towards the end of the pendulum swinging for 50 seconds which
% is considerably large. You can see it in the graph towards the end as the
% oscillating motion gets bigger and bigger amplitudes. This is all
% numerical error. Since there is no damping it should be a very constant
% oscillating function with the same amplitude for all times. When looking
% at Runge Kutta's second order approach we see a much smaller error of
% around 0.0002739 radians. Runge Kutta even of the second order is a much
% more accurate way to approximate ODEs as opposed to the traditional Euler's
% method.
