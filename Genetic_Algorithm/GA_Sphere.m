clc
clear all;
close all;
%%%Problem Parameters
global NFE;
NFE=0;

cost=@(x) Cost(x);

nvar=5;

varmax= 10;
varmin=-10;
varsize=[1 nvar];
%%
%GA Parameters

maxit=500;
npop=200;

pc=0.8;
nc=2*round((pc*npop)/2);

pm=0.3;
nm=round(pm*npop);

beta=5.5;

TournamentSize=3;

mu=0.05;
%%
%Initialization
chromosome.position=[];
chromosome.cost=[];

c=repmat(chromosome,npop,1);

for i=1:npop
    c(i).position=unifrnd(varmin,varmax,1,nvar);
    c(i).cost=cost(c(i).position);
end

Costs=[c.cost];
[Costs,SortOrder]=sort(Costs);
c=c(SortOrder);

bestsol=c(1);

bestcost=zeros(maxit,1);
worstcost=c(end).cost;
nfe=zeros(maxit,1);
%%
%GA Main Loop
for it=1:maxit
    %Calculate Selction Probabities
    P=exp(-beta*(Costs/worstcost));
    P=P/sum(P);
    
    %Crossover
    popc=repmat(chromosome,nc/2,2);
    for k=1:nc/2
         
        i1=RWS(P);
        i2=RWS(P);
        %i1=TournamentSelection(c,TournamentSize);
        %i2=TournamentSelection(c,TournamentSize);
        
        
        %Copy
        c1=c(i1);
        c2=c(i2);
        
        [popc(k,1).position,popc(k,2).position]=Crossover(c1.position,c2.position,varmin,varmax);
        
        popc(k,1).cost=cost(popc(k,1).position);
        popc(k,2).cost=cost(popc(k,2).position);
       
    end
    popc=popc(:);
    
    %Mutation
    popm=repmat(chromosome,nm,1);
    for k=1:nm
        
        i=randi([1 npop]);
        
        popm(k).position=Mutation(c(i).position,mu,varmin,varmax);
        popm(k).cost=cost(popm(k).position);
    end
    
    %Merge
    c=[c
       popc
       popm];
    
    %Sort
    Costs=[c.cost];
    [Costs,SortOrder]=sort(Costs);
    c=c(SortOrder);
    
    %Truncate
    c=c(1:npop);
    Costs=Costs(1:npop);
    
    %Stor Best Solution
    bestsol=c(1);
    
    %Stor Best & Worst Cost
    bestcost(it)=bestsol.cost;
    worstcost=max(c(end).cost,worstcost);
    nfe(it)=NFE;
   
    
    disp(['Iteratin ' num2str(it) ' Best Cost ' num2str(bestcost(it)) ' NFE ' num2str(nfe(it))]);
end




%%
%Results
figure;
semilogy(nfe,bestcost,'LineWidth',2);
xlabel('NFE');
ylabel('Best Cost');

%%
%Function
function z=Cost(x)
    
    global NFE;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;
    z=sum(x.^2);
end

function [y1,y2]=Crossover(c1,c2,varmin,varmax)
    
    gama=0.1;
    alfa=unifrnd(-gama,1+gama,size(c1));
   
    y1=alfa.*c1+(1-alfa).*c2;
    y2=(1-alfa).*c1+alfa.*c2;
    
    y1=min(max(y1,varmin),varmax);
    y2=min(max(y2,varmin),varmax);
end

function y=Mutation(x,mu,varmin,varmax)
    
    
    n=numel(x);
    
    SelectedN=ceil(n*mu);
    
    List=randsample(n,SelectedN);
    sigma=0.1*(varmax-varmin);
    y=x;
    y(List)=x(List)+sigma*randn(size(List));
    y=min(max(y,varmin),varmax);
    
end

function i=RWS(P) 
    
    r=rand;
    c=cumsum(P);
    i=min(find(r<c,1,'first'));
    
end

function i=TournamentSelection(c,m)
    
    n=numel(c);
    S=randsample(n,m);
    
    SPop=c(S);
    
    sc=[SPop.cost];
    
    [~,j]=min(sc);
    i=S(j);
end