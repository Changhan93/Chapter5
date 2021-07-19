% This code is used to estimate parameter kr1, a, gama, q1, K, r1, r2.

% Written by Changhan He and Juan Melendez-Alvarez in 2021.

clc
clear
global kr0 kr1 dr 
global kg J 
global c1 c2 c
global a  
global gama q1 q2 r1 r2 K

% Fixed parameters
J=2.123;
kr0=0.0514;  
kg=1.2540;
c1=1; 
c2=0; 
c=6.7;
q2=0.25;
dr=0.4813;

% Initial guess for the fiting parameters
% Here we select the best fitted parameters as the initial guess
Para0=[1.76784, 10.6323, 2.2911, 0.7234, 0.0862, 4.67, 1.633];

% Load and scale experimental data
rol=importdata('M9LB.mat');
Time_EXP=rol.time; 
rol.GFP_Lara(:,6) = rol.GFP_Lara(:,6)*0.001;
rol.OD_Lara(:,6) = rol.OD_Lara(:,6)*1;
GFP_EXP20=(rol.GFP_Lara(:,6))/(rol.GFP_Lara(1,6)); 
OD_EXP20=rol.OD_Lara(:,6);

Score_Best=1000000;
Para_best=Para0;
Score0=Score_Best;

%--------------------------------------------------------------------------
% Iteration time = 5000
for step=1:5000
  
    step
    pause(.000000000000001)
    
    Para1=Para0.*(.95+.1*rand(size(Para0)));
 
    temp=num2cell(Para1);
    [kr1, a, gama, q1, K, r1, r2]=deal(temp{:});
     
    % Solve the system in two steps
    % Step 1: solving the system with random AraC initial concentration so
    % that the circuit expression will reach its steady state
    % Step 2: Setting the steady state of circuit expression as the initial
    % condition
    sol20=ode23s(@ODESystem_new,[0 20],[0.1  OD_EXP20(1) 1]);    
    sol20=ode23s(@ODESystem_new,[0 20],[sol20.y(1,end)  OD_EXP20(1) 1]);
    
    OD_Sim=deval(sol20,Time_EXP,2);
    GFP_Sim=deval(sol20,Time_EXP,1);
    GFP_Sim=GFP_Sim/GFP_Sim(1);
    
    % Compute fitting error
    Score1=5*sum(((OD_Sim'-OD_EXP20)./OD_EXP20).^2+1*((GFP_Sim'-GFP_EXP20)./GFP_EXP20).^2);
 
    if Score1<Score0 || rand < .01 
        Para0=Para1;
        Score0=Score1;
    end
    
   
    if Score1<Score_Best
        
        Score_Best=Score1;
        Para_best=Para1;
        
        C20=sol20.y(1,:);
        N20=sol20.y(2,:);
        P20=sol20.y(3,:);
    
        % Plotting figures
        subplot(3,2,3)
        plot(sol20.x,N20,'linewidth',2)
        hold on
        plot(Time_EXP, OD_EXP20,'o')
        ylabel('N')
        xlabel('Time')
        hold off
        
        
        subplot(3,2,4)
        plot(sol20.x,C20/C20(1),'linewidth',2)
        hold on
        plot(Time_EXP, GFP_EXP20,'o')
        hold off
        ylabel('AraC')
        xlabel('Time')
        
        subplot(3,2,5)
        plot(sol20.x,P20)
        ylabel('Nutrient')
        xlabel('Time')


    end
    
end