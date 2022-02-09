clc;
clear all;
close all;
%% Problem Parameters

cost=@(x) Cost(x);

nvar=5;

varsize=[1 5];

varmin=-10;
varmax= 10;

%% FireFly Parameters

maxit=1000;

npop=100;

gamma=1;
alpha=0.2;
delta=0.5*(varmax-varmin);
m=1.2;
beta=2;

B=@(x1,x2,m,gamma,beta) beta*exp(-gamma*(norm(x1-x2))^m);

%% FireFly Setup

firefly.position=zeros(varsize);
firefly.cost=zeros(varsize);

pop=repmat(firefly,npop,1);

for i=1:npop
    
    pop(i).position=unifrnd(varmin,varmax,varsize);
    
    pop(i).cost=cost(pop(i).position);

end

Costs=[pop.cost];
[~,so]=sort(Costs);
pop=pop(so);

bestcost=zeros(maxit);

%% Main Loop

for it=1:maxit

    for i=1:npop
        for j=1:npop
               if pop(i).cost<pop(j).cost 
                    pop(i).position=pop(i).position...
                        +B(pop(j).position,pop(i).position,m,gamma,beta)...
                        *(pop(j).position-pop(i).position)...
                        +alpha*delta*randn(varsize);
                    
                    pop(i).position=min(max(pop(i).position,varmin),varmax);
                    
                    pop(i).cost=cost(pop(i).position);
               end
        end
    end

    Costs=[pop.cost];
    [~,so]=sort(Costs);
    pop=pop(so);
    
    alpha=alpha*(1-exp(it/maxit));
    
    bestcost(it)=pop(1).cost;
    
    disp(['Iteration : ' num2str(it)  ' Cost : ' num2str(bestcost(it))]);
    
end

loglog(bestcost);
%% Functions

function z=Cost(x)
    z=sum(x.^2);
end