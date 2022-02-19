clc;
clear all;
close all;
%% Problem Parameter

cost=@(x) Sphere(x);

nvar=5;
varsize=[1 nvar];

varmin=-10*ones(1,nvar);
varmax= 10*ones(1,nvar);

%% KH Parameters

maxit=100;
npop=20;

vf=0.02;
Dmax=0.005;
Dt=mean(abs(varmax-varmin))/2;
Nmax=0.01;
sr=0;
sf=zeros(varsize);
%% Setup KH

krill.position=[];
krill.N=0;
krill.F=0;
krill.d=0;
krill.dx=0;
krill.cost=[];

pop=repmat(krill,npop,1);


for i=1:npop
    
    pop(i).position=unifrnd(varmin,varmax);
    pop(i).cost=cost(pop(i).position);
    
end

pop=Sorting(pop);

bestsol=pop(1);
worstsol=pop(end);

Pfood_position=zeros(varsize);
Pfood_cost=inf;

bestcosts=zeros(maxit,1);
%%  Main Loop

for it=1:maxit

    BWcostD=bestsol.cost-worstsol.cost;
    %%% Virtual Food %%%
    
    positions=[pop.position];
    positions=reshape(positions,npop,nvar);
    Costs=[pop.cost];
    
    for j=1:nvar
        sf(j)=sum(positions(:,j)/Costs(j));     
    end
    
    food_position=sf/(sum(Costs));
    food_position=min(max(food_position,varmin),varmax);
    food_cost=cost(food_position);
    
    if food_cost<Pfood_cost
        Pfood_cost=food_cost;
        Pfood_position=food_position;
    else
        food_cost=Pfood_cost;
        food_position=Pfood_position;
    end
    
    w=(0.1+0.8*(1-it/maxit));
    
    for i=1:npop
        
        %%% Distances %%%
        
        food_dis=sqrt(sum((food_position-pop(i).position).^2));
        Bdis=sqrt(sum((bestsol.position-pop(i).position).^2));
        Distance=zeros(npop,npop);
        
        for m=1:npop
            for n=m:npop
                
                Distance(m,n)=sqrt(sum((pop(m).position-pop(n).position).^2));
                Distance(n,m)=Distance(m,n);
                
            end
        end
    
        %%% Best Krill Effect %%%
        
        if bestsol.cost < pop(i).cost
            alpha_b= -2*(1+rand*(it/maxit))*(bestsol.cost-pop(i).cost)/BWcostD;
        else
            alpha_b=0;
        end
    
        %%% Neihgbors Kill Effect %%%
        
        NeiNumber=0;
        ds = mean(Distance(i,:))/5;
        alpha_n=0;
        
        for l=1:npop
            NeiNumber=NeiNumber+1;
            if NeiNumber<=4 && Distance(i,l)<ds && i~=l
                alpha_n = alpha_n-(pop(l).cost-pop(i).cost)/(BWcostD*Distance(i,l));
            end        
        end
    
        %%% Movement Induced %%%
        pop(i).N=w*pop(i).N+Nmax*(alpha_b+alpha_n);
        
        %%% FOOD Attraction %%%
        
        if food_cost<pop(i).cost
            beta_f= -2*(1-it/maxit)*(food_cost-pop(i).cost)/(BWcostD*food_dis);
        else
            beta_f=0;
        end
        
        pop(i).F= w*pop(i).F+vf*(beta_f);
        
        %%% Physical Diffusion %%% 
        pop(i).d = Dmax*(1-it/maxit)*floor(rand+(pop(i).cost-bestsol.cost)/BWcostD)*(2*rand);
        
        pop(i).dx = pop(i).N + pop(i).F + pop(i).dx;
        
        pop(i).position = pop(i).position + (varmax-varmin) * pop(i).dx;
        pop(i).position = min(max(pop(i).position,varmin),varmax);
        
        pop(i).cost = cost(pop(i).position);
        
    end
    
    pop=Sorting(pop);

    bestsol=pop(1);
    worstsol=pop(end); 
    
    bestcosts(it)=bestsol.cost;
    
    disp(['Iteration: ' num2str(it) ' Cost: ' num2str(bestcosts(it))]);
    
end


%% Functions

function z=Sphere(x)
    z=sum(x.^2);
end

function pop=Sorting(pop)
    [~,so]=sort([pop.cost]);
    pop=pop(so);
end