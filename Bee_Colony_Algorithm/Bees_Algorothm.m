clc;
clear all;
close all;
%% Problem Parameters

cost=@(x) Sphere(x);

nvar=5;
varsize=[1 nvar];

varmin=-10;
varmax= 10;

%% Bees Parameter

maxiter=10000;

nScoutBee=30;

nSelectedSite=round(nScoutBee/2);

nEliteSite=round(0.4*nSelectedSite);

nSelectedSiteBee=round(0.5*nScoutBee);

nEliteSiteBee=2*nSelectedSiteBee;

r=0.2*(varmax-varmin);

r_damp=0.99;
%% Seup Bees

emptyBee.position=[];
emptyBee.cost=[];

pop=repmat(emptyBee,nScoutBee,1);

for i=1:nScoutBee
    
    pop(i).position=unifrnd(varmin,varmax,varsize);
    pop(i).cost=cost(pop(i).position);
    
end

[~,so]=sort([pop.cost]);
pop=pop(so);

bestsol=pop(1);

bestcost=zeros(maxiter,1);

%% Main Loop

for it=1:maxiter
    
    % Elit Sites
    for i=1:nEliteSite
        for j=1:nEliteSiteBee
           
          newbee=pop(i);
          
          newbee.position=UniformDance(newbee.position,r);
          newbee.position=min(max(newbee.position,varmin),varmax);
          
          newbee.cost=cost(newbee.position);
          
          if newbee.cost<pop(i).cost
              pop(i)=newbee;
          end
            
        end
    end
        
    % Selected Non-Elit Sites
    for i=nEliteSite+1:nSelectedSite
                   
        newbee=pop(i);

        newbee.position=UniformDance(newbee.position,r);
        newbee.position=min(max(newbee.position,varmin),varmax);

        newbee.cost=cost(newbee.position);

        if newbee.cost<pop(i).cost
            pop(i)=newbee;
        end
        
    end
    
    % Non-Selected Sites
    for i=nSelectedSite+1:nScoutBee
        
        pop(i).position=UniformDance(pop(i).position,r);
        pop(i).position=min(max(pop(i).position,varmin),varmax);
        pop(i).cost=cost(pop(i).position);
        
    end
    [~,so]=sort([pop.cost]);
    pop=pop(so);

    bestsol=pop(1);
    
    bestcost(it)=bestsol.cost;
    
    r=r*r_damp;
    
    disp(['Iteration: ' num2str(it) ' Cost: ' num2str(bestcost(it))]);
    
end

%% Plot

semilogy(bestcost);

%% Functions 
function z=Sphere(x)
    z=sum(x.^2);
end

function y=UniformDance(x,r)
    
    nvar=numel(x);
    
    k=randi([1 nvar]);
    
    y=x;
    y(k)=x(k)+unifrnd(-r,r);
end