clc;
clear all;
close all;
%% Problem Parameters

cost=@(x) Sphere(x);

nvar=5;
varsize=[1 nvar];

varmax= 10;
varmin=-10;
%% ACOR Parameters

npop=50;
maxit=10000;

q=0.5;
zeta=1;
nSample=50;
%% Setup ACOR

ant.position=[];
ant.cost=[];

pop=repmat(ant,npop,1);

for i=1:npop
   
    pop(i).position=unifrnd(varmin,varmax,varsize);
    pop(i).cost=cost(pop(i).position);
    
end

[~,so]=sort([pop.cost]);

pop=pop(so);

bestsol=pop(1);

bestcost=zeros(1,maxit);

%Weights
l=1:npop;
w=1/(sqrt(2*pi)*q*npop)*exp(-0.5*((l-1)/(q*npop)).^2);

%Probabilities
p=w/sum(w);

%% Main Loop

for it=1:maxit
    
    s=zeros(npop,nvar);
    
    for i=1:npop
       s(i,:)=pop(i).position;      
    end
    
    sigma=zeros(npop,nvar);
    
    for l=1:npop
        for i=1:nvar
            D=0;
            for r=1:npop
               D=D+abs(s(l,i)-s(r,i));
            end    
           sigma(l,i)=zeta*D/(npop-1); 
        end
    end
    
    newpop=repmat(ant,nSample,1);
    
    for t=1:nSample
        
        newpop(t).position=zeros(varsize);
        
        for i=1:nvar
            
            j=Roulett(p);
            
            newpop(t).position(i)=s(j,i)+sigma(j,i)*randn;
            
        end
        
       
       newpop(t).cost=cost(newpop(t).position);
       
    end
    
    pop=[pop;newpop];

    [~,so]=sort([pop.cost]);
    pop=pop(so);

    pop=pop(1:npop);
    
    bestsol=pop(1);
    
    bestcost(it)=bestsol.cost;
    
    disp(['iteration: ' num2str(it) ' Cost: ' num2str(bestcost(it))]);
    
end
%% Functions

function z=Sphere(x)
    z=sum(x.^2);
end

function i=Roulett(x)
    
    r=rand;
    c=cumsum(x);
    
    i=find(r<=c,1,'first');
end
