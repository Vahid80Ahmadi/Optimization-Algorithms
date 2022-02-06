clc;
clear all;
close all;
%% Problem Parameters

cost=@(x) costFUNC(x);

nvar=5;

varsize=[1 nvar];

varmin=-10;
varmax= 10;

%% BBO Parameters

maxiter=2000;

npop=1000;

keeprate=0.3;  %keep rate

nkeep=round(keeprate*npop); % number of kept habitat

nnew=npop-nkeep;

mu=linspace(1,0,npop); % emmigration rates
lambda=1-mu;           % immigration rates

alpha=0.9;

pm=0.2;

sigma=0.02*(varmax-varmin);

%% Setup BBO

habitat.position=zeros(varsize);
habitat.cost=0;

pop=repmat(habitat,npop,1);

for i=1:npop

    pop(i).position=unifrnd(varmin,varmax,varsize);
        
    pop(i).cost=cost(pop(i).position);
    
end

pop=sorting(pop);

bestcosts=zeros(maxiter,1);
%% Main Loop

for it=1:maxiter

    % Migration
    % Copy
    newpop=pop;
    for i=1:npop
        for k=1:nvar
            if rand <= lambda(i)
                Ep=mu;
                Ep(i)=0;
                Ep=Ep/sum(Ep);
                
                j=roulet(Ep);
                
                newpop(i).position(k)=pop(i).position(k)...
                    +alpha*(pop(j).position(k)-pop(i).position(k));
                
            end
            % Mutation

            if rand <= pm
                    newpop(i).position(k)=pop(i).position(k)...
                        +sigma*randn;
            end
                
        end
        newpop(i).position=min(max(newpop(i).position,varmin),varmax);
        newpop(i).cost=cost(pop(i).position);
    end
    
    newpop=sorting(newpop);
    
    pop=[pop(1:nkeep)
         newpop(1:nnew)];
     
    pop=sorting(pop);
    
    bestcosts(it)=pop(1).cost;
    
    disp(['iteration: ' num2str(it) ' Cost: ' num2str(pop(1).cost)]);
    
    
end

loglog(bestcosts);
%% Functions
function z=costFUNC(x)
    z=sum(x.^2);
end

function pop=sorting(pop)

    Costs=[pop.cost];
    [~,so]=sort(Costs);
    pop=pop(so);

end

function j=roulet(p)

    i=rand;
    c=cumsum(p);
    j=find(i<=c,1,'first');
    
end