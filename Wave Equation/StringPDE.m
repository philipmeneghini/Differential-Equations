%%
% This code uses finite differencing to solve the wave PDE. For
% simplicity's sake I am simulating only a 1D wave through time. The wave
% has an initial displacement of a triangular wave right in the middle of
% the string and I will show how the wave propogates across the string with
% initially both ends of the string fixed for 500 seconds. The user inputs
% an alpha value which gives me a dt for the approximation and then I put
% all the dislacement values in a matrix with the columns representing the
% string at different points through time (so each column represents the
% string at a distinct time state). I then graph this using a mesh plot and
% observe what happens with different alpha and dt values.

clear all

L = 100; %length of string in cm
v=0.5; %wave velocity of string (cm/s)

dx=0.5;  % spatial grid spacing of string (in cm)

alpha= input('what alpha value do you want');
%user inputs alpha
dt = sqrt((dx*dx*alpha)/(v^2))% time step (in sec)
%calculated because alpha= (v*dt/dx)^2 rearrange that to find dt
tfinal = 500;   % total time for simulation (in sec)

X = 0:dx:L;       % space grid for string
T = 0:dt:tfinal;   % time grid

M = length(X);         % number of grid points in length
N = length(T);         % number of steps in time


u=zeros(M,N);   % intialize the displacement array, setting its size

%IC
%set first two time states with triangular wave as specified above
%since this is the initial displacement. I need to do this for
%first two time steps because my finite differencing iteration refers to
%time state two behind one I am trying to calculate.
for k=1:M
    if k*dx >46 && k*dx<50.5
        u(k,2)=u(k-1,2)+(1/8); %set triangular wave from 46 cm to 50 increasing by (1/8)
         u(k,1)=u(k,2);
    else if k*dx >= 50.5 && k*dx<=54
        u(k, 2)= u(k-1, 2)-(1/8); %finish the triangular wave by going back down to zero displacement
         u(k,1)=u(k,2);
    else
        u(k, 2)=0;
        u(k,1)=0;
        end
    end
end

%Dirichlet BC:
u(1, :)= 0; %For both ends of string
u(M, :)=0;

figure(1)
%Shows initial triangular pulse with a plot of displacement vs. x (length
%of string)
plot(X,u(:,1));
title('Initial Displacement')
xlabel('Length of String (cm)')
ylabel('Displacement of String (cm)')

%These boundary conditions cause
%the wave to invert over the displacement axis. So the wave points up and
%then when it reaches the end of the string the wave flips over and points
%down with a negative displacement. 

for n= 2: N-1      %time loop 
     for j= 2: M-1  %space loop
         %finite differencing iteration uses u(j, n-1) as said above in why
         %we need two initial condition states
         u(j,n+1) = 2*(1-alpha)*u(j,n)-u(j, n-1)+alpha*(u(j+1, n)+u(j-1, n));
     end       
end
 


 figure(2)
 mesh(T, X, u); %mesh plot showing triangular wave moving through time and across the 1d string
 title('Numerical Approximation of 1D String Through Time')
 xlabel('Time (sec)')
 ylabel('String Length (cm)')
 zlabel('Displacement of String (cm)')
 colorbar
 
%%

% When alpha equals 0.25 dt equals 0.5 and the wave behaves normally with
% a very nice mesh plot showing the wave moving exactly the way it is
% suppose to through time. Nothing in our numerical approximation of the
% PDE is acting weird. When I change alpha to 1.0 dt changes to 1.0 and I
% can see in the mesh plot how the bigger dt makes the numerical solution
% less accurate and more sparse. (for 0.25 the wave looked continuous going
% through time while when alpha and dt =1 I can see the discrete time steps
% we have taken). Other than that the solution is pretty good and their is
% no weird behavior from my numerical approximation. As soon as I move the
% alpha value a little above one to 1.0025 all credence to this numerical
% approximation is destroyed. The dt because lower than alpha which is what
% changes when going above one and in this case of alpha (1.0025) the data
% all spikes up at the beginning time step to make an oval shape and then
% as it goes through time the string has no displacement (the string just
% stays on zero). So this is obviously terribly wrong and a terrible
% representation of the wave PDE because the wave doesnt move through time
% and the initial wave insn't even trigular like it is suppose to be it is
% a weird circle shape.
% The values of alpha above zero obviously and at or below one give good
% approximations of the wave but an alpha value anywhere above one will give
% wildly wrong approximations for the PDE. It essentially breaks our
% numerical model. When alpha equals 1.0 is a boundary for our numerical 
% solution because it is where the time step gets twice as big as dx 
% (the way we laid outour spatial grid) of 0.5. This makes the dt too big 
% and leads to nonsensical answers.
