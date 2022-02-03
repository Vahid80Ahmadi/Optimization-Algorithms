clc;
clear all; 
close all;
%% Problem Parameters

model=CreateModel();
n=model.n;

cost=@(x) Cost(x,model);

%% ACO Parameters

maxiter=300;

npop=100;

alpha=1;
beta=0.1;

rho=0.05;

tau0=1;

Q=1;
%% Setup ACO

N=[0 1];


eta=[model.w./model.v
     model.v./model.w];

tau=tau0*ones(2,n);

ant.path=[];
ant.cost=[];
ant.x=[];
ant.sol=[];

pop=repmat(ant,npop,1);

bestsol.path=[];
bestsol.cost=inf;

bestcost=zeros(maxiter,1);
%% Main Loop

for it=1:maxiter

    for k=1:npop
        
        pop(k).path=[];
        for l=1:n
            
            p=(tau(:,l).^alpha).*(eta(:,l).^beta);
            
            p=p/sum(p);
            
            j=Roulett(p);
            
            pop(k).path=[pop(k).path j];
            
        end
        
        pop(k).x=N(pop(k).path);
        
        [pop(k).sol,pop(k).cost]=cost(pop(k).x);
        
        
        if pop(k).cost<bestsol.cost
            bestsol=pop(k);
        end

    end
        
    %Update Phromones

    for k=1:n
        
        tour=pop(k).path;

        for l=1:n
            
            tau(tour(l),l)=tau(tour(l),l)+Q/pop(k).cost;
            
        end

    end    
    
    % Evaporation
    
    tau=(1-rho)*tau;
    
    bestcost(it)=bestsol.cost;
    
    disp(['iteration: ' num2str(it) ' Cost: ' num2str(bestcost(it))]);
    
    

    
end


%% Functions
function model=CreateModel()
    
    v=[426 463 150 466 353 139 211 319 483 486 163 489 483 294 420 156 269 467 417 484];
    w=[53 21 63 67 54 58 57 40 53 28 56 21 34 22 24 61 55 36 68 21];
    
    W=500;
    n=numel(w);
    
    model.n=n; 
    model.W=W;
    model.w=w;
    model.v=v;
    
end

function [sol,z]=Cost(x,model)
    
    v=model.v;
    w=model.w;
    W=model.W;
    
    V1=sum(v.*x);
    W1=sum(w.*x);
    V0=sum(v.*(1-x));
    W0=sum(w.*(1-x));
    
    Vio=max(W1/W-1,0);
    
    beta=100000;
    
    z=V0+beta*Vio;
    
    sol.V1=V1;
    sol.W1=W1;
    sol.V0=V0;
    sol.W0=W0;
    sol.IsFeasible=(Vio==0);
    sol.W=W;
    sol.Vio=Vio;
end

function i=Roulett(p)
    
    r=rand;
    c=cumsum(p);
    i=find(r<=c,1,'first');
end
