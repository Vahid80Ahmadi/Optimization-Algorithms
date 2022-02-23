clc;
clear;
close all;
%%
%Probelm Parameters

cost=@(x) MOP2(x);

nvar=3;
varsize=[1 nvar];
varmax= 4;
varmin=-4;

%%
%MOPSO Parameters

maxit=400;
npop=200;
nrep=50;    %Archive Size

w=1;
wdamp=0.99;
c1=1;
c2=2;

ngrid=5;    %Number Of Grids Pre Dimension
alpha=0.2;  %Infilatin Rate   
beta=5;     %Leader Selection Pressure
gamma=5;    %Deletion Selection Pressure
mu=0.02;    %Mutation Rate
%%
%Setup MOPSO
particle.position=[];
particle.velocity=[];
particle.cost=[];
particle.best.position=[];
particle.best.cost=[];
particle.isdominated=[];
particle.gridindex=[];
particle.gridsubindex=[];

p=repmat(particle,npop,1);

for i=1:npop
    
    p(i).position=unifrnd(varmin,varmax,varsize);
    p(i).velocity=zeros(varsize);
    p(i).cost=cost(p(i).position);
        
    p(i).best.position=p(i).position;
    p(i).best.cost=p(i).cost;
    
end

%Determine Domination
p=DominateFlag(p);

%Repository
rep=p(~[p.isdominated]);

%Making Grid
grid=CreateGrid(rep,ngrid,alpha);

for i=1:numel(rep)
   rep(i)= FindGridIndex(rep(i),grid);
end
%%
%MOPSO Main Loop

for it=1:maxit
       
    for i=1:npop
        
        %Select Leader
        leader=SelectLeader(rep,beta);
        
        %Update Velocity
        p(i).velocity=w*p(i).velocity...
            +c1*rand(varsize).*(p(i).best.position-p(i).position)...
            +c2*rand(varsize).*(leader.position-p(i).position);
        
        %Update Position
        p(i).position=p(i).position+p(i).velocity;
        
        %Evaluate
        p(i).cost=cost(p(i).position);
        
        %Apply Mutation
        pm=(1-((it-1)/(maxit-1)))^(1/mu);
        newsol.position=Mutate(p(i).position,pm,varmin,varmax);
        newsol.cost=cost(newsol.position);
        
        if Dominate(p(i).cost,newsol.cost)
            newsol.cost=p(i).cost;
            newsol.position=p(i).position;
        
        elseif Dominate(newsol.cost,p(i).cost)
            %Do Nothing
        else
            if rand<0.5
                newsol.cost=p(i).cost;
                newsol.position=p(i).position;  
            end
        end
        
        %Update Personal Best
        if Dominate(p(i).cost,p(i).best.cost)
            p(i).best.cost=p(i).cost;
            p(i).best.position=p(i).position;
        
        elseif Dominate(p(i).best.cost,p(i).cost)
            %Do Nothing
        else
            if rand<0.5
                p(i).best.cost=p(i).cost;
                p(i).best.position=p(i).position;  
            end
        end
        
    end
    
    
    %Repository
    
    %Add Non-Dominated Particles to REPOSITORY    
    rep=[rep; p(~[p.isdominated])];

    %Determine Domination of New Resository Members
    rep=DominateFlag(rep);
    
    %Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.isdominated]);

    %Update Grid
    grid=CreateGrid(rep,ngrid,alpha);

    %Update Grid Indices
    for i=1:numel(rep)
        rep(i)= FindGridIndex(rep(i),grid);
    end
    
    %Check Repository Is Full
    if numel(rep)>nrep
       Extra=numel(rep)-nrep;
       
       for r=1:Extra
          rep=DeleteRepository(rep,gamma); 
       end
    end
    
    w=w*wdamp;
    
    %Show
    disp(['Iteration ' num2str(it)...
        ' Number Of Rep Member ' num2str(numel(rep))]);
    figure(1);
    MyPlot(p,rep);
    
    
end
%%
%Functions

function z=MOP2(x)

    n=numel(x);
    
    z1=1-exp(-sum((x-1/sqrt(n)).^2));
    
    z2=1-exp(-sum((x+1/sqrt(n)).^2));
    
    z=[z1 z2]';

end

function b=Dominate(x,y)

    b=(all(x<=y) && any(x<y));

end

function pop=DominateFlag(pop)
    
    npop=numel(pop);
    
    for i=1:npop
        pop(i).isdominated=false;
    end
    
    for i=1:npop-1
       for j=i+1:npop
           
           if Dominate(pop(i).cost,pop(j).cost)
                pop(j).isdominated=true;
           end
           if Dominate(pop(j).cost,pop(i).cost)
                pop(i).isdominated=true;
           end
           
       end
    end
    
    
end

function MyPlot(pop,rep)
    
    pop_costs=[pop.cost];
    rep_costs=[rep.cost];
  
    plot(pop_costs(1,:),pop_costs(2,:),'ko');
    hold on;
    plot(rep_costs(1,:),rep_costs(2,:),'r*');
    
    xlabel('1st Objective');
    ylabel('2nd Objective');
    %grid on;
    hold off;
    
end

function grid=CreateGrid(pop,ngrid,alpha)

    c=[pop.cost];
    
    cmin=min(c,[],2);
    cmax=max(c,[],2);
    
    dc=cmax-cmin;
    cmin=cmin-alpha*dc;
    cmax=cmax+alpha*dc;
    
    nobj=size(c,1);
    
    empty_grid.lb=[];
    empty_grid.ub=[];
    
    grid=repmat(empty_grid,nobj,1);
    
    for j=1:nobj
        
       cj=linspace(cmin(j),cmax(j),ngrid+1);
       
       grid(j).lb=[-inf cj];
       grid(j).ub=[cj +inf];
       
    end
    
    
end

function particle=FindGridIndex(particle,grid)
    
    nobj=numel(particle.cost);
    ngrid=numel(grid(1).ub);
    
    particle.gridsubindex=zeros(1,nobj);
    
    for i=1:nobj
        
        particle.gridsubindex(i)...
            =find(particle.cost(i)<grid(i).ub,1,'first');

    end
    
    particle.gridindex=particle.gridsubindex(1);
    for j=2:nobj
        particle.gridindex=particle.gridindex-1;
        particle.gridindex=ngrid*particle.gridindex;
        particle.gridindex=particle.gridindex+particle.gridsubindex(j);
    end
    
end

function leader=SelectLeader(rep,beta)
    
    % Grid Index Of All Repository Members
    GI=[rep.gridindex];
    
    % Occupied Cells
    OC=unique(GI);
    
    %Number of Particle in Occupied Cells
    N=zeros(size(OC));
    for k=1:numel(OC)
        N(k)=numel(find(GI==OC(k)));
    end
    
    %Selection Probabilities
    P=exp(-beta*N);
    P=P/sum(P);
    
    %Selected Cell Index
    sci=RWS(P);
    
    %Selected Cell
    sc=OC(sci);
    
    %Selected Cell Members
    SCM=find(GI==sc);
    
    %Selected Member Index
    smi=randi([1 numel(SCM)]);
    sm=SCM(smi);
    
    %Select Leader
    leader=rep(sm);
    
end

function i=RWS(p)
    
    r=rand;
    c=cumsum(p);
    
    i=find(r<=c,1,'first');
end

function rep=DeleteRepository(rep,gamma)
  
    % Grid Index Of All Repository Members
    GI=[rep.gridindex];
    
    % Occupied Cells
    OC=unique(GI);
    
    %Number of Particle in Occupied Cells
    N=zeros(size(OC));
    for k=1:numel(OC)
        N(k)=numel(find(GI==OC(k)));
    end
    
    %Selection Probabilities
    P=exp(gamma*N);
    P=P/sum(P);
    
    %Selected Cell Index
    sci=RWS(P);
    
    %Selected Cell
    sc=OC(sci);
    
    %Selected Cell Members
    SCM=find(GI==sc);
    
    %Selected Member Index
    smi=randi([1 numel(SCM)]);
    sm=SCM(smi);
    
    %Delete Selected Member
    rep(sm)=[];
    
end

function xnew=Mutate(x,pm,varmin,varmax)
    
    nvar=numel(x);
    j=randi([1 nvar]);
    
    dx=pm*(varmax-varmin);
    
    lb=x(j)-dx;
    if lb<varmin
        lb=varmin;
    end
    
    ub=x(j)+dx;
    if ub>varmax
        ub=varmax;
    end
    
    xnew=x;
    xnew(j)=unifrnd(lb,ub);  
        
end