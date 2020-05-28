% % Based on the paper:
% % M. Ghasemi, J. Aghaei, and M. Hadipour, "New self-organising  hierarchical PSO with jumping time-varying acceleration
% % coefficients," Electron. Lett., vol. 53, no. 20, pp. 1360–1362, 2017. DOI: 10.1049/el.2017.2112.

clear
clc
disp('HPSO-TVAC');

CostFunction=@(x) ObjFun(x);

nVar=40;            % Number of decision variables

VarMin=-100;        % Minimum value of decision variables
VarMax=-VarMin;     % Maximum value of decision variables

nPop=40;            % Number of algorithm population

ci=0.5;
cf=0.0;

MaxNFE=100e3;           % Maximum Number of Function Evaluations
MaxIt=MaxNFE/nPop;

dx=VarMax-VarMin;
vmax(1:nVar)=0.5*dx;

velocity=zeros(nPop,nVar);
position=velocity;
cost=zeros(1,nVar);
gbestcost=inf;
pbest=position;
pbestcost=cost;
gcost=zeros(1,MaxIt);
NFES=zeros(1,MaxIt);
NFE=0;
it=0;

while NFE<MaxNFE
    it=it+1;
    if it==1    % Initialization

        for i=1:nPop
            velocity(i,1:nVar)=vmax;
            
            position(i,:)=rand(1,nVar).*VarMin+(VarMax-VarMin);
            cost(i)=CostFunction(position(i,:));
            NFE=NFE+1;
            
            pbest(i,:)=position(i,:);
            pbestcost(i)=cost(i);
            
            if pbestcost(i)<gbestcost(it)
                
                gbest=pbest(i,:);
                gbestcost=pbestcost(i);
            end
        end
    else
        c_it=((cf-ci)*(it/MaxIt))+ci;
        
        for i=1:nPop
            A=randperm(nPop);
            A(A==i)=[];
            a0=A(1);
            
            w=randn;
            c1_it=(abs(w))^((c_it)*(w));
            c2_it=(abs((1-w)))^((c_it)*(1/(1-w)));
            %%%%%%%%%%%%%HPSO
            velocity(i,:)=(c1_it*(rand(1,nVar)).*((pbest(i,:))-position(i,:)))...
                +(c2_it*(rand(1,nVar)).*((gbest+pbest(a0,:))-2*position(i,:)));
            
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


