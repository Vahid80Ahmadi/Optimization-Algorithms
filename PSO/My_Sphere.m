clc;
clear all;
close all;
%%
%Problem Parameters

cost=@(x) costfunction(x);

nvar=5;
varsize=[1,nvar];
varmax=[2 2 2 2 2];
varmin=[-1 -1 -1 -1 -1];

%%
%PSO Parameters


npop=1000;
iter=10000;
% 
% w=0.7298;
% wp=0.99;
% c1=1.4192;
% c2=2;
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));
w=chi;          % Inertia Weight
wp=1;        % Inertia Weight Damping Ratio
c1=chi*phi1;    % Personal Learning Coefficient
c2=chi*phi2;    % Global Learning Coefficient
velmax=0.1*(varmax-varmin);
velmin=-velmax;
%%
%Setup PSO
particle.position=[];
particle.object=[];
particle.velocity=[];
particle.best.position=[];
particle.best.object=[];

globalbest.position=[];
globalbest.object=inf;

p=repmat(particle,npop,1);

for i=1:npop
    
    p(i).position=unifrnd(varmin,varmax);
    p(i).velocity=zeros(varsize);
    p(i).object=cost(p(i).position);
    
    p(i).best.position=p(i).position;
    p(i).best.object=p(i).object;
    
    if p(i).best.object<globalbest.object
       globalbest= p(i).best; 
    end
    
end
%%
%PSO Main Loop
bestcost=zeros(iter,1);

for it=1:iter

    for i=1:npop
    
        %Update Velocity
        p(i).velocity=w*p(i).velocity+c1*rand(varsize).*(p(i).best.position-p(i).position)...
            +c2*rand(varsize).*(globalbest.position-(p(i).position));
        
        %Apply Velocity Limits
        p(i).velocity=min(max(p(i).velocity,velmin),velmax);
        
        %Update Position
        p(i).position=p(i).position+p(i).velocity;
        
        %Apply Position Limits
        p(i).position=min(max(p(i).position,varmin),varmax);
        
        %Mirroring
        IsOutside=(p(i).position<varmin | p(i).position>varmax);
        p(i).velocity(IsOutside)=-p(i).velocity(IsOutside);
        
        %Evaluation
        p(i).object=cost(p(i).position);
        
        
        if p(i).object<p(i).best.object
            
            p(i).best.position=p(i).position;
            p(i).best.object=p(i).object;
            
            if p(i).best.object<globalbest.object
                
                globalbest=p(i).best;
            end
        end   
    end
    bestcost(it)=globalbest.object;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(bestcost(it))]);
    
    w=w*wp;
    
end
%%
%Plot Results
semilogy(bestcost);

%%
%Function
function z=costfunction(x)
     z=sum(x.^2);
end