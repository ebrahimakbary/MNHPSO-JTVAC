% % Based on the paper:
% % Mojtaba Ghasemi, Ebrahim Akbari, Mohammad Zand, Morteza Hadipour,
% % Sahand Ghavidel, and Li Li. "An Efficient Modified HPSO-TVAC-Based
% % Dynamic Economic Dispatch of Generating Units." Electric Power
% % Components and Systems 47, no. 19-20 (2019): 1826-1840.

clear
tic
clc
disp('MNHPSO-JTVAC');

CostFunction=@(x) ObjFun(x);

nVar=30;            % Number of decision variables

VarMin=-100;        % Minimum value of decision variables
VarMax=-VarMin;     % Maximum value of decision variables

nPop=100;           % Number of algorithm population

npop_0=nPop;
npop_f=20;

ci=0.5;
cf=0.0;

dx=VarMax-VarMin;
vmax(1:nVar)=0.5*dx;

MaxNFE=100e3;           % Maximum Number of Function Evaluations
MaxIt=MaxNFE/npop_f;

velocity=zeros(nPop,nVar);
position=velocity;
cost=zeros(1,nVar);
gbestcost=inf;
pbest=position;
pbestcost=cost;
gcost=nan(1,MaxIt);
NFES=zeros(1,MaxIt);
NFE=0;
it=0;

while NFE<=MaxNFE
    it=it+1;
    
    if NFE==0
        for i=1:nPop
            velocity(i,1:nVar)=vmax;
            position(i,:)=VarMin+rand(1,nVar).*(VarMax-VarMin);
            cost(i)=CostFunction(position(i,:));
            NFE=NFE+1;
            pbest(i,:)=position(i,:);
            pbestcost(i)=cost(i);
            
            if pbestcost(i)<gbestcost
                gbest=pbest(i,:);
                gbestcost=pbestcost(i);
            end
        end
        gcost(1)=gbestcost;
    else
        c_it=((cf-ci)*(NFE/MaxNFE))+ci;
        
        % % Population Reduction Policy
        nPop=npop_0-(npop_0-npop_f)*(NFE/MaxNFE);
        nPop=round(nPop);
        nPop=max(nPop,npop_f);
        nPop=min(npop_0,nPop);
        
        for i=1:nPop
            
            % % Select a random particle other than the i-th one
            a0=randi(nPop-1,1);
            if a0>=i
                a0=a0+1;
            end
            
            w=randn;
            c1_it=(abs(w))^(c_it*w);
            c2_it=(abs(1-w))^(c_it*(1/(1-w)));
            %%%%%%%%%%%%% 
            velocity(i,:)=(c1_it*(rand(1,nVar)).*((pbest(i,:))-position(i,:)))...
                +(c2_it*(rand(1,nVar)).*((gbest(1,:)+pbest(a0,:))-2*position(i,:)));
            
            for gg=1:nVar
                if velocity(i,gg)==0
                    if rand<0.5
                        velocity(i,gg)=rand*vmax(1,gg);
                    else
                        velocity(i,gg)=-rand*vmax(1,gg);
                    end
                end
                velocity(i,gg)=sign(velocity(i,gg)).*min(abs(velocity(i,gg)),abs(vmax(1,gg)));
            end
            
            %%%%%%%%%%%%%%%%
            
            velocity(i,:)=min(max(velocity(i,:),-vmax),vmax);
            position(i,:)=position(i,:)+velocity(i,:);
            position(i,:)=min(max(position(i,:),VarMin),VarMax);
            cost(i)=CostFunction(position(i,:));
            NFE=NFE+1;
            
            if cost(i)<pbestcost(i)
                pbest(i,:)=position(i,:);
                pbestcost(i)=cost(i);
                if pbestcost(i)<gbestcost
                    gbest=pbest(i,:);
                    gbestcost=pbestcost(i);
                end
            end
        end
    end
    NFES(it)=NFE;
    gcost(it)=gbestcost;
    fprintf('Iteration %3.0f,  NFE %6.0f,   Best Cost %g\n',it,NFE,gbestcost)
end
gcost(it+1:end)=[];
NFES(it+1:end)=[];
plot(NFES,gcost,'b','linewidth',2);

BestSol.Position=gbest;
BestSol.Cost=gbestcost;
