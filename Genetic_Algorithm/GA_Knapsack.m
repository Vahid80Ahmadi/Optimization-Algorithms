clc;
clear all;
close all;
%%
%Problem Parameters
global NFE;
NFE=0;

model=CreateModel();
cost=@(x) Cost(x,model);

nvar=model.N;

varmax=1;
varmin=0;
varsize=[1 nvar];
%%
%GA Parameters

maxit=200;
npop=100;

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
chromosome.sol=[];

c=repmat(chromosome,npop,1);

for i=1:npop
    c(i).position=CreatRandomSolution(model);
    [c(i).cost,c(i).sol]=cost(c(i).position);
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
        
        
        c1=c(i1);
        c2=c(i2);
        
        [popc(k,1).position,popc(k,2).position]=Crossover(c1.position,c2.position);
        
        [popc(k,1).cost,popc(k,1).sol]=cost(popc(k,1).position);
        [popc(k,2).cost,popc(k,2).sol]=cost(popc(k,2).position);
       
    end
    popc=popc(:);
    
    %Mutation
    popm=repmat(chromosome,nm,1);
    for k=1:nm
        
        i=randi([1 npop]);
        
        popm(k).position=Mutation(c(i).position,mu);
        [popm(k).cost,popm(k).sol]=cost(popm(k).position);
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
plot(nfe,bestcost,'LineWidth',2);
xlabel('NFE');
ylabel('Best Cost');
grid on;

%%
%Functions
function model=CreateModel()
    
    N=10;
    vi=[4 20 20 11 17 4 10 19 17 20];
    w=[857 408 537 324 236 535 853 169 500 456];
    W=2100;
    
    model.N=N;
    model.vi=vi;
    model.w=w;
    model.W=W;
end

function x=CreatRandomSolution(model)
    
    N=model.N;
    x=randi([0 1],[1 N]);
    
end

function [z,sol]=Cost(x,model)
    
    global NFE;
    if isempty(NFE)
        NFE=0;
    end
    
    NFE=NFE+1;
    
    vi=model.vi;
    W=model.W;
    w=model.w;
    
    Violation=max((sum(x.*w)/W)-1,0);
    alpha=1000000;
    
    z=sum(vi.*(1-x))+alpha*Violation;
    
    sol.Violation=Violation;
    sol.Weight=sum(x.*w);
    sol.LostWeight=sum((1-x).*w);
end

function [y1,y2]=Crossover(c1,c2)
    
    Selection=randi([1 3]);
    
    switch Selection
        case 1
            [y1,y2]=SinglePointCrossover(c1,c2);
        case 2
            [y1,y2]=DoublePointCrossover(c1,c2);
        case 3
            [y1,y2]=UniformCrossover(c1,c2);
    end
end

function [y1,y2]=SinglePointCrossover(c1,c2)
    
    n=numel(c1);
    cutpoint=randi([1,n-1]);
    y1=[c1(1:cutpoint) c2(cutpoint+1:end)];
    y2=[c2(1:cutpoint) c1(cutpoint+1:end)];
    
end

function [y1,y2]=DoublePointCrossover(c1,c2)
    
    n=numel(c1);
    cutpoint=randsample(n-1,2);
    cc1=min(cutpoint);
    cc2=max(cutpoint);
    
    y1=[c1(1:cc1) c2(cc1+1:cc2) c1(cc2+1:end)];
    y2=[c2(1:cc1) c1(cc1+1:cc2) c2(cc2+1:end)];
    
end

function [y1,y2]=UniformCrossover(c1,c2)
    
    n=numel(c1);
    
    alfa=randi([0 1],1,n);
    
   y1=alfa.*c1+(1-alfa).*c2;
   y2=(1-alfa).*c1+alfa.*c2;

end

function y=Mutation(x,mu)
    
    pim=mu;
    n=numel(x);
    
    SelectedN=ceil(n*pim);
    
    List=randsample(n,SelectedN);
   
    x(List)=1-x(List);
    y=x;
end

function i=RWS(P) 
    
    r=rand;
    c=cumsum(P);
    i=min(find(r<c,1,'first'));
    
end
