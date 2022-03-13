clc;
clear all; 
close all;
%%
%Peroblem parameters

model=CreateModel();
cost=@(p) Cost(p,model);


nvar=model.N;
varmax=model.pmax;
varmin=model.pmin;
varsize=[1,nvar];

%%
%PSO Parameters
w=1;
wp=0.99;
c1=2;
c2=2;

velmax=0.1*(varmax-varmin);
velmin=-velmax;

npop=200;
iter=1000;
%%
%Setup PSO
particle.position=[];
particle.velocity=[];
particle.cost=[];
particle.best.position=[];
particle.best.cost=[];

globalbest.position=[];
globalbest.cost=inf;

p=repmat(particle,nvar,1);


for i=1:npop
    
    p(i).position=CreateRandomSolution(model);
    p(i).velocity=zeros(varsize);
    p(i).cost=cost(p(i).position);
    p(i).best.position=p(i).position;
    p(i).best.cost=p(i).cost;
    
    if p(i).best.cost<globalbest.cost
        globalbest=p(i).best;
    end
end
%%
%PSO Main Loop
bestcost=zeros(iter,1);

for it=1:iter

    for i=1:npop
        
        p(i).velocity=w*p(i).velocity...
            +c1*rand(varsize).*(p(i).best.position-p(i).position)...
            +c2*rand(varsize).*(globalbest.position-p(i).position);

        p(i).velocity=min(max(p(i).velocity,velmin),velmax);
        
        p(i).position=p(i).position+p(i).velocity;
        
        IsOutside=(p(i).position<varmin | p(i).position>varmax);
        p(i).velocity(IsOutside)=-p(i).velocity(IsOutside);
        
        p(i).position=min(max(p(i).position,velmin),velmax);
        
        
        p(i).cost=cost(p(i).position);
        
        if p(i).cost<p(i).best.cost
            p(i).best.position=p(i).position;
            p(i).best.cost=p(i).cost;
           if p(i).best.cost<globalbest.cost
               globalbest=p(i).best;
           end
            
        end
        
    end
    bestcost(it)=globalbest.cost;
    
    disp(['Iteration:' num2str(it) '  Best Cost:' num2str(bestcost(it))]);
    w=w*wp;
    
    
end
%%
%results
semilogy(bestcost);


%%
%Functions
function model=CreateModel()
    
    pmin=[808 853 463 857 716 448 539 673];
    pmax=[2166 2172 1526 2177 2166 1788 2041 1513];
    a0=[7067 9488 8882 9702 8213 5175 9161 9577];
    a1=[9 9 9 7 8 6 9 5];
    a2=[-0.1092 -0.1194 -0.2647 -0.2390 -0.1634 -0.2900 -0.1069 -0.1877]*1.0e-4;
    
    N=numel(pmin);
    
    pL=10000;
    
    model.N=N;
    model.pmin=pmin;
    model.pmax=pmax;
    model.a0=a0;
    model.a1=a1;
    model.a2=a2;
    model.pL=pL;
 
end
function p=CreateRandomSolution(model)

    pmin=model.pmin;
    pmax=model.pmax;
    
    p=unifrnd(pmin,pmax);
    
end
function  z=Cost(p,model)
 
    z=sum(model.a0+model.a1.*(p)+model.a2.*(p.^2));
    
    %violation
    
    v=abs((sum(p)/model.pL)-1);
    
    beta=2; 
    
    z=z*(1+beta*v);
    
end