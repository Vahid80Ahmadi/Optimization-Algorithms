clc;
clear all;
close all;
%%
%Problem Paarameters

cost=@(x) sphere(x);
nvar=100;
varsize=[1 nvar];

varmin=-100;
varmax= 100;
%%
%GWO Parameters

maxit=100;
npop=50;


%%
%Setup GWO

empty.position=[];
empty.cost=inf;

alpha=empty;
beta=empty;
delta=empty;

pop=repmat(empty,npop,1);

for i=1:npop
    
    pop(i).position=unifrnd(varmin,varmax,varsize);
    pop(i).cost=cost(pop(i).position);
    
    if pop(i).cost<alpha.cost
        alpha=pop(i);
        if pop(i).cost<beta.cost
            beta=pop(i);
            if pop(i).cost<delta.cost
                delta=pop(i);
            end
        end
    end
    
end

%%
% Main Loop

for it=1:maxit
    
   a=2-it*(2/maxit); 
   for i=1:npop
       
       %Alpha Part
       r1=rand(varsize);
       r2=rand(varsize);
    
       A=a*(2*r1-1);
    
       c=2*r2;
        
       d=abs(c.*alpha.position-pop(i).position);
       x1=alpha.position-A.*d;
       
       
       %Beta Part
       r1=rand(varsize);
       r2=rand(varsize);
    
       A=a*(2*r1-1);
    
       c=2*r2;
        
       d=abs(c.*beta.position-pop(i).position);
       x2=beta.position-A.*d;
       
       
       %Delta Part
       r1=rand(varsize);
       r2=rand(varsize);
    
       A=a*(2*r1-1);
    
       c=2*r2;
        
       d=abs(c.*delta.position-pop(i).position);
       x3=delta.position-A.*d;
       
       pop(i).position=(x1+x2+x3)/3;  
       pop(i).position=min(max(pop(i).position,varmin),varmax);
%        FlagUb = pop(i).position>varmax;
%        FlagLb = pop(i).position<varmin;
%         
%        pop(i).position = (pop(i).position.*(~(FlagUb+FlagLb))) ...
%             + varmax.*FlagUb + varmin.*FlagLb;
       % Commented Way Was Wonderful
       pop(i).cost=cost(pop(i).position);
       

    if pop(i).cost<alpha.cost
        alpha=pop(i);
        if pop(i).cost<beta.cost
            beta=pop(i);
            if pop(i).cost<delta.cost
                delta=pop(i);
            end
        end
    end
      
   end
   disp(['Iteration ' num2str(it) ' Cost : ' num2str(alpha.cost)]);
   
end

%%
%Functions
function z=sphere(x)
    z=sum(x.^2);
end