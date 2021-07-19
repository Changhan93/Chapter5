clc
clear
global kr0 kr1 dr 
global kg J 
global c1 c2 c
global a  
global gama q1 q2 r1 r2 K

J=2.123;
kr0=0.0514;  
kg=1.2540;
c1=1; 
c2=0; 
c=6.7;

% r1=2.75; 
% r2=1.75; 
% q1=0.40; 
q2=0.25;
% gama=2.7167;

a=7.6334;
dr=0.4813;
kr1=0.7630;

Para0=[[1.5842, 5.3552, 1.3008, 0.7247, 0.09181, 4.6606, 1.0499]];
% Para0=[1.76784, 10.6323, 2.2911, 0.6234, 0.0862, 3.67, 1.633]

rol=importdata('M9LB.mat');
Time_EXP=rol.time; 

rol.GFP_Lara(:,7) = rol.GFP_Lara(:,7)*0.001;
rol.OD_Lara(:,7) = rol.OD_Lara(:,7)*1;

GFP_EXP40=(rol.GFP_Lara(:,7))/(rol.GFP_Lara(1,7)); 
OD_EXP40=rol.OD_Lara(:,7);

Score_Best=1000000;
Para_best=Para0;
Score0=Score_Best;
%%
%for step=1:5000
for step=1:1
    
    step
    pause(.000000000000001)
    
%    Para1=Para0.*(.95+.1*rand(size(Para0)));
    Para1=Para0;
   
    temp=num2cell(Para1);
    [kr1, a, gama, q1, K, r1, r2]=deal(temp{:});
     
%     kr0=0.1*R_Timescale; dr=R_Timescale;
%     
%     syms Rm
%     eqn(Rm) = kr0+kr1*kg*Rm^n/(Rm^n+K^n)-kg*Rm-dr*Rm== 0;
%     sol=vpasolve(eqn, Rm,[0 Inf]);
%     Rmax=max(double(sol)); % the maximum ribosome number at growth rate=kg
%       
%     
%     Rmin=kr0/dr; % the ribosome number at growth rate=0
%     R0=Rmin/2; % at growth rate=0, half of the ribosome can be used for cell maintenance
%         
%     k0=0.1/(Rmin-R0)*100*C_Timescale;k1=2/(Rmin-R0)*100*C_Timescale;dc=1*100*C_Timescale;

    
    %%
 
    sol40=ode23s(@ODESystem_new,[0 20],[0.1  OD_EXP40(1) 1]);    
    sol40=ode23s(@ODESystem_new,[0 20],[sol40.y(1,end)  OD_EXP40(1) 1]);
    
    OD_Sim=deval(sol40,Time_EXP,2);
    GFP_Sim=deval(sol40,Time_EXP,1);
    GFP_Sim=GFP_Sim/GFP_Sim(1);

    Score1=2*sum(((OD_Sim'-OD_EXP40)./OD_EXP40).^2 +1*((GFP_Sim'-GFP_EXP40)./GFP_EXP40).^2);
    
    Score1
    
    figure(1)
    if Score1<Score0 || rand < .01 %exp(-(Score1-Score0)*1)
        Para0=Para1;
        Score0=Score1;
%         subplot(3,1,1)
%         plot(step,Score1,'o')
%      %   hold on
%      hold off
    end
    
   
    if Score1<Score_Best
        
        Score_Best=Score1;
        Para_best=Para1;
        
        C40=sol40.y(1,:);
        N40=sol40.y(2,:);
        P40=sol40.y(3,:);
    
%         subplot(3,2,3)
%       %  hold off
%         plot(sol40.x,N40,'linewidth',2)
%         hold on
%         plot(Time_EXP, OD_EXP40,'o')
%         ylabel('N')
%         xlabel('Time')
%         hold off
%         
%         
%         subplot(3,2,4)
%      %   hold off
%         plot(sol40.x,C40/C40(1),'linewidth',2)
%         hold on
%         plot(Time_EXP, GFP_EXP40,'o')
%         hold off
%         ylabel('AraC')
%         xlabel('Time')
%         
%         subplot(3,2,5)
%      %   hold off
%         plot(sol40.x,P40)
%         ylabel('Nutrient')
%         xlabel('Time')
        
        
        figure(1)
yyaxis left
plot(sol40.x,N40,'linewidth',2)
hold on
plot(Time_EXP, OD_EXP40,'o','Markersize',5)
hold off

yyaxis right
plot(sol40.x,C40/C40(1),'linewidth',2)
hold on
plot(Time_EXP, GFP_EXP40,'o','Markersize',5)
hold off

yyaxis left
xlabel('Time (Hour)')
ylabel('OD')

yyaxis right
ylabel('GFP')

legend

    end
    
end