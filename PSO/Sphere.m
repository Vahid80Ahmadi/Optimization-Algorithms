clc;
clear all;
close all;
%% 
%Problem definition

Costfunction=@(x) sphere(x);
%Costfunction=@sphere;

nvar=5; %Number of decition variables
varsize=[1 nvar];

varmin=[1 1 1 1 1];
varmax=[2 2 2 2 2];
%%
%PSO parameters

maxit=1000;

npop=20;


w=1;
wdamp=0.99;
c1=2;
c2=2;

velmax=0.1*(varmax-varmin);
velmin=-velmax;

%%
%Initialization

empty_particle.position=[];
empty_particle.cost=[];
empty_particle.velocity=[];
empty_particle.best.position=[];
empty_particle.best.cost=[];

p=repmat(empty_particle,npop,1);

globalbest.cost=inf;

for i=1:npop
    
   p(i).position=unifrnd(varmin,varmax);
   
   p(i).velocity=zeros(varsize);
   
   p(i).cost=Costfunction(p(i).position);
   
   p(i).best.position=p(i).position;
   p(i).best.cost=p(i).cost;
   
   if p(i).best.cost<globalbest.cost
       %globalbest.cost=p(i).best.cost;
       %globalbest.position=p(i).best.cost;
       globalbest=p(i).best;
   end
   
end

bestcost=zeros(maxit,1);

%%
%PSO main loop

for it=1:maxit
    for i=1:npop
        %update velocity
        p(i).velocity=w.*(p(i).velocity)...
            +c1.*rand(varsize).*(p(i).best.position-p(i).position)...
            +c2.*rand(varsize).*(globalbest.position-p(i).position);
        
        %Apply velocity limits
        
        p(i).velocity=min(max(p(i).velocity,velmin),velmax);
        
        %update position
        p(i).position=p(i).position+p(i).velocity;
        
        %Velocity mirror effect
        
        isoutside=(p(i).position<varmin | p(i).position>varmax);
        p(i).position(isoutside)=-p(i).position(isoutside);
        
        %Apply position limits
        p(i).position=min(max(p(i).position,varmin),varmax);
        
        
        
        %evaluation
        p(i).cost=Costfunction(p(i).position);
        
        %update personal best
        if p(i).cost<p(i).best.cost
            p(i).best.cost=p(i).cost;
            p(i).best.position=p(i).position;
        
        
            %update global best
            if p(i).best.cost<globalbest.cost
                globalbest=p(i).best;
            end
        end
        
            
    end
    bestcost(it)=globalbest.cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(bestcost(it))]);
    
    w=w*wdamp;
end
%%
%Result

figure;
semilogy(bestcost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');

%%
%Func
function z=sphere(x)
    z=sum(x.^2);
end