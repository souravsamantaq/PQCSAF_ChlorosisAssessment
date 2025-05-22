%---------------------------------------------------------------------------
% The multics_feature_selection() function is designed for feature selection 
% using a multi-swarm cuckoo search optimization approach.
%===========================================================================
function multics_feature_selection()
Max_iteration=10;
dim=408;
SearchAgents_no=20;
mos=301;
nos=550;
fobj=@obj_fs;

%---------------------------------
GF=0.04;
%---------------------------------
ds=load('dataset.mat');
data=ds.data;
datalabel=ds.label;
orgdataset=data;
%-----------------------------------
data1=horzcat(data(:,1:102),datalabel);
data2=horzcat(data(:,103:204),datalabel);
data3=horzcat(data(:,205:306),datalabel);
data4=horzcat(data(:,307:408),datalabel);

data_org=data;

dim_1=102;
dim_2=102;
dim_3=102;
dim_4=102;
dim=102; %data=data4;

%----------------------------------------------
[bswarm]=initialization_cs_init(dim);
dataswram1.data=data1;
dataswram1.bswarm=bswarm;
%--------------------------------------------
[bswarm]=initialization_cs_init(dim);
dataswram2.data=data2;
dataswram2.bswarm=bswarm;
%--------------------------------------------
[bswarm]=initialization_cs_init(dim);
dataswram3.data=data3;
dataswram3.bswarm=bswarm;
%--------------------------------------------
[bswarm]=initialization_cs_init(dim);
dataswram4.data=data4;
dataswram4.bswarm=bswarm;
%----------------------------------------------
main_dataswarmList{1}=dataswram1;
main_dataswarmList{2}=dataswram2;
main_dataswarmList{3}=dataswram3;
main_dataswarmList{4}=dataswram4;
%----------------------------------------------


funList = {@BCS_SUB_B,@BCS_SUB_B,@BCS_SUB_B,@BCS_SUB_B};

tic

t=1; cs_gen=1;
%===========================main multi cs===========================================
%-----------------------------------------
%parpool local
%tic
while(t<40)

for i=1:4
outList{i}=funList{i}(main_dataswarmList{i});
end

%toc/60
%------------------------------------------

SF1=outList{1}.bgfitsol;
SF2=outList{2}.bgfitsol;
SF3=outList{3}.bgfitsol;
SF4=outList{4}.bgfitsol;

SF=[SF1 SF2 SF3 SF4];

[data,k]=Optimum_Feature_Selection(SF,data_org);
dataM=horzcat(data(1:250,:),datalabel);

[bcsm]=BCS_MAIN(dataM,k);

MF=bcsm.gbest;
MF_SOL=bcsm.gfitsol;
BMF_SOL=bcsm.bgfitsol;
  %-------------------------------------
     if(MF<GF)
         GF=MF;
         MF_SOL=MF_SOL;
         GMF_SOL=BMF_SOL;
         GSF_SOL=SF;
     end
     MF;
     Perf(cs_gen,1)=GF;
     Perf(cs_gen,2)=MF;
     Perf(cs_gen,3)=nnz(BMF_SOL);
     cs_gen=cs_gen+1;
  %-------------------------------------
  main_dataswarmList{1}.bswarm(1,:)=rand(1,102);
  main_dataswarmList{1}.bswarm(5,:)=rand(1,102);
  
  main_dataswarmList{2}.bswarm(6,:)=rand(1,102);
  main_dataswarmList{2}.bswarm(2,:)=rand(1,102);
  
  main_dataswarmList{3}.bswarm(8,:)=rand(1,102);
  main_dataswarmList{3}.bswarm(4,:)=rand(1,102);
  
  main_dataswarmList{4}.bswarm(4,:)=rand(1,102);
  main_dataswarmList{4}.bswarm(5,:)=rand(1,102);
  %--------------------------------------
   t=t+10;
end   
toc/60;
%==============================================================================
%delete(gcp('nocreate'))


SF=GSF_SOL;
SZ=GMF_SOL;
nnz(SZ);
%--------------------------------------
nos=250;
%orgdataset=data(1:nos,1:408);
[orgdata1,k]=Optimum_Feature_Selection(SF,orgdataset);
[orgdata2,k]=Optimum_Feature_Selection(SZ,orgdata1);
%---------------------------------------
% dataL=dataset(1:nos,409);
% orgdataFLSVM=horzcat(orgdata2,dataL);
dataLN=datalabel(1:nos,:);
orgdataFLNN=horzcat(orgdata2,dataLN);

%----------------------------------------------------------------
end

function [bcs]=BCS_MAIN(data,dim)
global qps ps  datalabel

Max_iteration=5;
%dim=102;
SearchAgents_no=20;
%nos=200;
fobj=@obj_fs;

%%-------------------------

bcs.gbest={};
bcs.gfitsol={};
bcs.bgfitsol={};
bcs.acc={};
bcs.fs={};
bcs.st={};
%-------------------------

ps=SearchAgents_no;

%-----------------Cuckoo search---------------------

[bswarm]=initialization_cs(dim);
[gbest,bgfitsol,gfitsol,fitarr]=objective_cs(bswarm,fobj,data);

%---------------------------------------
t=1;
tic
while(t<=Max_iteration)
%     
      [bswarm]=get_cuckoos(bswarm,gbest,gfitsol,fitarr);
     [xgbest,xbgfitsol,xgfitsol,xfitarr]=objective_cs(bswarm,fobj,data);
%       xgfs=nnz(xgfitsol);
%        %if((xgbest<=gbest)&& (xgfs<gfs))
       if((xgbest<=gbest))
           gbest=xgbest;
           bgfitsol=xbgfitsol;
           gfitsol=xgfitsol;
           %gfs=xgfs;
       end
%        
%        
%        
        gfitstack(t)=gbest;
%        
%       [xgbest,xgfitsol,xfitarr,bswarm]=empty_nest(bswarm,fobj,xfitarr);
        [bswarm]=empty_nest(bswarm,fobj,xfitarr);
%        
%        
% %         if((xgbest<=gbest)&& (xgfs<gfs))
         %if((xgbest<=gbest))
         %  gbest=xgbest;
         %  bgfitsol=xbgfitsol;
         %  gfitsol=xgfitsol;
%           gfs=xgfs;
        %end
%        xgfs=nnz(xgfitsol)
%        (1-gbest)*100
%     itr=itr+1
    listfit3(t)=gbest;
    t=t+1;
% 
end
t_2=toc/60;

    acc=1-gbest;
    
    bcs.gbest=gbest;
    bcs.gfitsol=gfitsol;
    bcs.bgfitsol=bgfitsol;
    bcs.acc=acc;
    bcs.fs=nnz(gfitsol);
    bcs.st=gfitstack;
%    save('qbcs_testing_mdataset.mat','qbcs')
%    
%    plot(qbcs.st)
%-------------------------------------------------------------------

CS_Best_score=gbest;
CS_Best_ant=bgfitsol;
CS_cg_curve=listfit3;
% CS_cg_curve=vertcat(listfit2,listfit3');

CS_CSDATA.Best_score=gbest;
CS_CSDATA.Best_chrosome=bgfitsol;
CS_CSDATA.cg_curve=CS_cg_curve;


%save('CS_CSDATA_SEV.mat','CS_CSDATA')
%-------------------------------------------------------------------
end

function [bcs]=BCS_SUB_B(dataswram)
global ps
bswarm=dataswram.bswarm;
data=dataswram.data;

Max_iteration=5;
dim=102;
SearchAgents_no=10;

fobj=@obj_fs;

%-------------------------

bcs.gbest={};
bcs.gfitsol={};
bcs.bgfitsol={};
bcs.acc={};
bcs.fs={};
bcs.st={};
%-------------------------

ps=SearchAgents_no;

%-----------------Cuckoo search---------------------

% [bswarm]=initialization_cs(dim);
[gbest,bgfitsol,gfitsol,fitarr]=objective_cs(bswarm,fobj,data);

%---------------------------------------
t=1;
tic
while(t<=Max_iteration)
%     
     [bswarm]=get_cuckoos(bswarm,gbest,gfitsol,fitarr);
     [xgbest,xbgfitsol,xgfitsol,xfitarr]=objective_cs(bswarm,fobj,data);
%       xgfs=nnz(xgfitsol);
%        %if((xgbest<=gbest)&& (xgfs<gfs))
       if((xgbest<=gbest))
           gbest=xgbest;
           bgfitsol=xbgfitsol;
           gfitsol=xgfitsol;
           %gfs=xgfs;
       end
%        
%        
%        
        gfitstack(t)=gbest;
%        
%       [xgbest,xgfitsol,xfitarr,bswarm]=empty_nest(bswarm,fobj,xfitarr);
        [bswarm]=empty_nest(bswarm,fobj,xfitarr);
%        
%        
% %         if((xgbest<=gbest)&& (xgfs<gfs))
         %if((xgbest<=gbest))
         %  gbest=xgbest;
         %  bgfitsol=xbgfitsol;
         %  gfitsol=xgfitsol;
%           gfs=xgfs;
        %end
%        xgfs=nnz(xgfitsol)
%        (1-gbest)*100
%     itr=itr+1
    listfit3(t)=gbest;
    t=t+1;
% 
end
t_2=toc/60;

    acc=1-gbest;
    
    bcs.gbest=gbest;
    bcs.gfitsol=gfitsol;
    bcs.bgfitsol=bgfitsol;
    bcs.acc=acc;
    bcs.fs=nnz(gfitsol);
    bcs.st=gfitstack;
%    save('qbcs_testing_mdataset.mat','qbcs')
%    
%    plot(qbcs.st)
%-------------------------------------------------------------------

CS_Best_score=gbest;
CS_Best_ant=bgfitsol;
CS_cg_curve=listfit3;
% CS_cg_curve=vertcat(listfit2,listfit3');

CS_CSDATA.Best_score=gbest;
CS_CSDATA.Best_chrosome=bgfitsol;
CS_CSDATA.cg_curve=CS_cg_curve;

%save('BGACS_GACSDATA_DIS_2.mat','BGACS_GACSDATA')
%save('CS_CSDATA_SEV.mat','CS_CSDATA')
%-------------------------------------------------------------------
end


%----------------------------get cuckoo-----------------------------
function [bswarm]=get_cuckoos(bswarm,gbest,gfitsol,fitarr)
global ps 
% Levy flights
nest=bswarm;
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
%V_shaped_transfer_function=abs((2/pi).*atan((pi/2).*qswarm(i,:)));
m=1;
a=ones(1,size(nest,2));
for j=1:1:ps
    s=bswarm(j,:);
    best=gfitsol;
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
    stepsize=1*step.*(s-best);
    as=s+stepsize.*randn(size(s));

%     tf=abs((2/pi).*atan((pi/2).*as));
%     tr=rand(1,dim);
%     nest(m,:)=(tr>tf);
%    as=a./(a+exp(-s));
    nest(m,:)=as;
    m=m+1;

end

bswarm=nest;
end



%-------------------------------------------------------------------

function [bswarm]=empty_nest(bswarm,fobj,fitarr)
%xgbest,xgfitsol,xfitarr,bswarm
global dim data
[ival,indx]=sort(fitarr);

fit_c=size(fitarr,2);
fit_cm=round(fit_c*0.8);  %0.6 to 0.8


for i=fit_c:-1:fit_cm
    idx=indx(i);
    
    as=rand(1,dim);
    bswarm(i,:)=as;

end

end


%----------------------------------------------------------------
function [fitval,bgfisol,gfitsol,fitarr]=objective_cs(bswarm,fobj,data)
%global data 
global ps dim 
dim=size(bswarm,2);
rpop=zeros(1,dim);
fitarr=zeros(1,ps);

tswarm=zeros(ps,dim);

for p=1:ps
    rpop(1,:)=bswarm(p,:);
    %---------------------
%     tf=abs((2/pi).*atan((pi/2).*rpop));
%     tr=rand(1,dim);
%     tswarm(p,:)=(tr>tf);
%     rpop=tswarm(p,:);
%-----------------------
     tf=sigmoid(rpop);
     %tf=abs((2/pi).*atan((pi/2).*rpop));
     tr=rand(1,dim);
     tswarm(p,:)=(tr>tf);
     rpop=tswarm(p,:);
     fitarr(p)=obj_fs(rpop,data);
end
    [fitval,ind]=min(fitarr);
    gfitsol=bswarm(ind,:);
    bgfisol=tswarm(ind,:);
  
end   
%----------------------------------------------------------------------
function [data,k]=Optimum_Feature_Selection(SF,dataset)
ld=size(SF,2);
[h w]=size(dataset);
k=1;
% for i=1:ld
%     if(SF(i)==1)
%         k=k+1;
%     end
% end
%----------------------------------
%data=zeros(h,k+15);
%----------------------------------
for i=1:ld
     if(SF(i)==1)
       data(:,k)=dataset(:,i);
       k=k+1;
     end
end
size(data);
size(dataset);
%data=horzcat(data,dataset(:,31:45));
%data=horzcat(data,dataset(:,171:174));
k=k-1;
end 
%---------------------------------------------------------------------

function [bswarm]=initialization_cs(dim)
global ps;
%global dim;

%qswarm=zeros(ps,dim);

bswarm=round(rand(ps,dim));
end

function [bswarm]=initialization_cs_init(dim)
%global ps;
ps=20;
%global dim;

%qswarm=zeros(ps,dim);

bswarm=round(rand(ps,dim));
end


function y = sigmoid(x)
    y = 1 ./ (1 + exp(-x));
end


function [fval]=obj_fs(s,data)

[v]=size(data,2);
    % Read Data Elements
%     x=data.x;
%     t=data.t;
     x=data(:,1:v-4)';
     t=data(:,v-3:v)';
%    x=data(:,1:v-1)';
%    t=data(:,v)';

%    t=data(:,171)';
%       x=X1';
%       t=Y1';

    % Selected Features
    S=find(s~=0);

    % Number of Selected Features
    nf=numel(S);
    
    % Ratio of Selected Features
    rf=nf/numel(s);
    
    % Selecting Features
    xs=x(S,:);
    
    % Weights of Train and Test Errors
    wTrain=0.8;
    wTest=1-wTrain;

    % Number of Runs
    nRun=2;
    EE=zeros(1,nRun);
    for r=1:nRun
        % Create and Train ANN
        results=CreateAndTrainANN(xs,t);

        % Calculate Overall Error
        EE(r) = wTrain*results.TrainData.E + wTest*results.TestData.E;
    end
    
    E=mean(EE);
    %if isinf(E)
    %    E=1e10;
    %end
    
    % Calculate Final Cost
    beta=0.5;
    z=E*(1+beta*rf);

    % Set Outputs
    out.S=S;
    out.nf=nf;
    out.rf=rf;
    out.E=E;
    out.z=z;
    out.net=results.net;
    out.net;
    ganet=out.net;
    %save ganet
    %view(out.net)
    %out.Data=results.Data;
    %out.TrainData=results.TrainData;
    %out.TestData=results.TestData;
    %x=rand(1,20);
    %predict(results.net,x)
    fval=z;
end



function results=CreateAndTrainANN(x,t)

    if ~isempty(x)
        
        % Choose a Training Function
        % For a list of all training functions type: help nntrain
        % 'trainlm' is usually fastest.
        % 'trainbr' takes longer but may be better for challenging problems.
        % 'trainscg' uses less memory. NFTOOL falls back to this in low memory situations.
        trainFcn = 'trainlm';  % Levenberg-Marquardt

        % Create a Fitting Network
        hiddenLayerSize = 10;
        net = fitnet(hiddenLayerSize,trainFcn);

        % Choose Input and Output Pre/Post-Processing Functions
        % For a list of all processing functions type: help nnprocess
        net.input.processFcns = {'removeconstantrows','mapminmax'};
        net.output.processFcns = {'removeconstantrows','mapminmax'};

        % Setup Division of Data for Training, Validation, Testing
        % For a list of all data division functions type: help nndivide
        net.divideFcn = 'dividerand';  % Divide data randomly
        net.divideMode = 'sample';  % Divide up every sample
        net.divideParam.trainRatio = 70/100;
        net.divideParam.valRatio = 15/100;
        net.divideParam.testRatio = 15/100;

        % Choose a Performance Function
        % For a list of all performance functions type: help nnperformance
        net.performFcn = 'mse';  % Mean squared error

        % Choose Plot Functions
        % For a list of all plot functions type: help nnplot
        net.plotFcns = {};
        % net.plotFcns = {'plotperform','plottrainstate','ploterrhist', 'plotregression', 'plotfit'};

        net.trainParam.showWindow=false;
        net.trainParam.epochs=50;
        
        % Train the Network
        [net,tr] = train(net,x,t);

        % Test the Network
        y = net(x);
        e = gsubtract(t,y);
        E = perform(net,t,y);
        
    else        
        
        y=inf(size(t));
        e=inf(size(t));
        E=inf;
        
        tr.trainInd=[];
        tr.valInd=[];
        tr.testInd=[];
        
    end

    % All Data
    Data.x=x;
    Data.t=t;
    Data.y=y;
    Data.e=e;
    Data.E=E;
    
    % Train Data
    TrainData.x=x(:,tr.trainInd);
    TrainData.t=t(:,tr.trainInd);
    TrainData.y=y(:,tr.trainInd);
    TrainData.e=e(:,tr.trainInd);
    if ~isempty(x)
        TrainData.E=perform(net,TrainData.t,TrainData.y);
    else
        TrainData.E=inf;
    end
    
    % Validation and Test Data
    TestData.x=x(:,[tr.testInd tr.valInd]);
    TestData.t=t(:,[tr.testInd tr.valInd]);
    TestData.y=y(:,[tr.testInd tr.valInd]);
    TestData.e=e(:,[tr.testInd tr.valInd]);
    if ~isempty(x)
        TestData.E=perform(net,TestData.t,TestData.y);
    else
        TestData.E=inf;
    end
    
    % Export Results
    if ~isempty(x)
        results.net=net;
    else
        results.net=[];
    end
    results.Data=Data;
    results.TrainData=TrainData;
    % results.ValidationData=ValidationData;
    results.TestData=TestData;
    
end