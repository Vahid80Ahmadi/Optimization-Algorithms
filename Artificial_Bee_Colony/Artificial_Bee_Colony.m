clc;
clear all; 
close all;
%% Problem Parameters

cost=@(x) Sphere(x);

nvar=10;
varsize=[1 nvar];

varmin=-10;
varmax= 10;

%% ABC Parameters

npop=300;
maxit=1000;


L=round(npop*nvar);

a=1;            % Acceleration Coeffitiont Upper Bound
damp=1.001;

%% Setup ABC

empty.position=[];
empty.cost=[];
empty.L=0;
empty.fitness=0;

bestsol.cost=inf;


pop=repmat(empty,1,npop);

for i=1:npop
    
   pop(i).position=unifrnd(varmin,varmax,varsize);
   pop(i).cost=cost(pop(i).position);
    
   if pop(i).cost<=bestsol.cost
       bestsol=pop(i);
   end
   
end
   
bestcost=zeros(maxit,1);

%% Main Loop

for it=1:maxit
    
    % Recruited Bee
    for i=1:npop
    
        newbee=pop(i);
        
        K=[1:i-1 i+1:npop];
        k=K(randi([1,numel(K)]));
        
        phi=unifrnd(-a,a,varsize);
        
        newbee.position=pop(i).position+...
            phi.*(pop(i).position-pop(k).position);
        
        newbee.position=min(max(newbee.position,varmin),varmax);
        
        newbee.cost=cost(newbee.position);
    
        if newbee.cost<pop(i).cost
            pop(i)=newbee;
        else
            pop(i).L=pop(i).L+1;
        end
        
    end
    
    % Fitness Value
    for i=1:npop
        if pop(i).cost>=0
            pop(i).fitness=1/(1+pop(i).cost);
        else
            pop(i).fitness=1+abs(pop(i).cost);
        end
    end
    
    F=[pop.fitness];
    
    p=F/sum(F);
    
    % Onlooker Bee
    for b=1:npop
        
        i=Roulett(p);
        
        K=[1:i-1 i+1:npop];
        k=K(randi([1,numel(K)]));
        
        phi=unifrnd(-a,a,varsize);
        
        newbee.position=pop(i).position+...
            phi.*(pop(i).position-pop(k).position);
        
        newbee.position=min(max(newbee.position,varmin),varmax);
        newbee.cost=cost(newbee.position);
    
        if newbee.cost<pop(i).cost
            pop(i)=newbee;
        else
            pop(i).L=pop(i).L+1;
        end
        
        
    end
    
    % Scout Bee
    for i=1:npop     
        
        if pop(i).L>=L
            pop(i).position=unifrnd(varmin,varmax,varsize);
            pop(i).cost=cost(pop(i).position);
            pop(i).L=0;
        end
        if pop(i).cost<bestsol.cost
            bestsol=pop(i);
        end
    end
    
    bestcost(it)=bestsol.cost;
     a=a*damp;
    disp(['iteration :' num2str(it) ' Cost :' num2str(bestcost(it))]);
    
end

%% Plot

loglog(bestcost);

%%  Functions
function z=Sphere(x)
    z=sum(x.^2);
end

function R=Roulett(x)
    r=rand;
    C=cumsum(x);
    
    R=find(r<=C,1,'first');
end