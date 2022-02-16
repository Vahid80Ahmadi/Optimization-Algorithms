clc;
clear all;
close all;
%%
%Problem Parameters

cost=@(x) Sphere(x);

nvar=5;
varsize=[1 nvar];

varmin=-10;
varmax= 10;


%%
%Harmony Search Parameters

maxit=4000;
HMS=60;      %Harmony Memory Size

nNew=200;    %Number Of New Harmonies 


PAR=0.1;
HMCR=0.5;
FW=0.02*(varmax-varmin); %Fret Width
damp=0.99;
%%
%Setup Harmony Search

empty.position=zeros(varsize);
empty.cost=[];

pop=repmat(empty,HMS,1);

for i=1:HMS

    pop(i).position=unifrnd(varmin,varmax,varsize);
    pop(i).cost=cost(pop(i).position);

end

Costs=[pop.cost];
[~,so]=sort(Costs);

pop=pop(so);

bestsol=pop(1);

bestcost=zeros(maxit);

%%
%Main Loop

for it=1:maxit
    
    new=repmat(empty,nNew,1);
    
    for k=1:nNew
        
        %Create New Harmony Search
        new(k).position=unifrnd(varmin,varmax,varsize);
        for j=1:nvar            
            if rand<=HMCR                
                i=randi([1 HMS]);
                new(k).position(j)=pop(i).position(j);
            end    
           % Pitch Adjustment
           if rand<=PAR
              delta=unifrnd(-FW,FW);
              new(k).position(j)=new(k).position(j)+delta;
              
           end
            
        end
        
        new(k).position=min(max(new(k).position,varmin),varmax);
        new(k).cost=cost(new(k).position);

    end
    
    
    
    %Merging
    
    pop=[pop
         new];
    
    [~,so]=sort([pop.cost]);
    pop=pop(so);
    
    pop=pop(1:HMS);
    
    bestsol=pop(1);
    bestcost(it)=bestsol.cost;
    
    disp(['Iteration: ' num2str(it) ' Best Cost: ' num2str(bestcost(it))]);

    FW=FW*damp;
    
end

%%
%Functions
function z=Sphere(x)
    z=sum(x.^2);
end


