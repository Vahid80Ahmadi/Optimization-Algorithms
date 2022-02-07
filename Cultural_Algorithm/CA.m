clc;
clear;
close all;
%% Problem Parameters

cost = @(x) Sphere(x);

nvar=5;
varsize=[1 nvar];

varmin=-10;
varmax= 10;

%% Cultural Algorithm Parameters

maxit=100;
npop=20;

pAccept=0.33;
nAccept=round(pAccept*npop);

alpha=0.25;
beta=1;

%% Setup Cultural Algorithm

% Empty culture
culture.situational.cost=inf;
culture.normative.size=[];
culture.normative.max=-inf(varsize);
culture.normative.min= inf(varsize);
culture.normative.L=inf(varsize);
culture.normative.U=inf(varsize);

% Empty individual Strcture
empty.position=[];
empty.cost=[];

pop=repmat(empty,npop,1);

for i=1:npop
    
    pop(i).position=unifrnd(varmin,varmax,varsize); 
    pop(i).cost=cost(pop(i).position);

end

pop=Sorting(pop);

spop=pop(1:nAccept);
culture=AdjustCulture(culture,spop);


bestsol=culture.situational;

bestcosts=zeros(maxit,1);
%% Main Loop
for it=1:maxit
    
    % Influence of Culture
    for i=1:npop
        
        sigma=alpha*culture.normative.size;
        pop(i).position=pop(i).position+sigma.*randn(varsize);
        
        pop(i).cost=cost(pop(i).position);
        
    end
    pop=Sorting(pop);

    spop=pop(1:nAccept);
    culture=AdjustCulture(culture,spop);

    bestsol=culture.situational;
    
    bestcosts(it)=bestsol.cost;
    
    disp(['Iteration : ' num2str(it) ' Cost : ' num2str(bestcosts(it))]);
    
end
%% Functions
function z=Sphere(x)
    z=sum(x.^2);
end

function pop=Sorting(pop)
    [~,SortOrder]=sort([pop.cost]);
    pop=pop(SortOrder);
end

function culture=AdjustCulture(culture,pop)
    
    n=numel(pop);
    nvar=numel(pop(1).position);
    
    for i=1:n
        
        if pop(i).cost < culture.situational.cost
            culture.situational=pop(i);
        end
        
        for j=1:nvar
            
            if pop(i).position(j) < culture.normative.min(j) ...
                  ||  pop(1).cost < culture.normative.L(j)
                
              culture.normative.min(j) = pop(i).position(j);
              culture.normative.L(j) = pop(i).cost;
            end
            
            if pop(i).position(j) > culture.normative.max(j) ...
                  ||  pop(1).cost < culture.normative.U(j)
                
              culture.normative.max(j) = pop(i).position(j);
              culture.normative.U(j) = pop(i).cost;
            end
  
        end
        
    end
    
    culture.normative.size=culture.normative.max-culture.normative.min;
    
end