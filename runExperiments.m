% Author: Matthew Blaschko - matthew.blaschko@esat.kuleuven.be
% code for arXiv:1606.05918

more off
addpath LPoptimizer

saveDir = '/Users/mblaschk/work/ConvexExtensions/Experiments/octaveTmp/';

useSavedMatFile = true;
fname1 = [saveDir 'baselineResults.mat'];
if(exist(fname1,'file')~=2 || useSavedMatFile==false)
    load ~/work/ConvexExtensions/Experiments/cifar-10-batches-mat/test_batch.mat
    clear batch_label labels
    dataSubsample = 400;
    dataInd = [ones(dataSubsample,1);zeros(size(data,1)-dataSubsample,1)];
    dataInd = find(dataInd(randperm(length(dataInd))));
 
    [func,D] = coverageFunction(data(dataInd,:));
    %dataInd = dataInd(coverageInd);

    numRandRepititions = 50;
    budget = round(dataSubsample*.4);

    for i=1:budget
      for j=1:numRandRepititions
        valRand(i,j) = func(randomSelection(func,dataSubsample,i));
      end
    end

    ind = coverageNemhauserWolsey(D,budget);
    valNemhauser = [];
    for i=1:length(ind)
      valNemhauser(i) = func(ind(1:i));
    end


    myGamma = findMaxGammaMarginRescaling(func,size(D,1));

    save(fname1);
else
    load(fname1);
end

useSavedMargin = true;
sparms = [];
xs = [];
valMargin = [];
valMarginLP = [];
budgetsLP = 1:budget;
for i=1:budget
  fname = [saveDir sprintf('marginResult%d.mat',i)];
  if(exist(fname,'file')~=2 || useSavedMargin==false)
    disp(sprintf(' running margin rescaling with budget %d %s',i,datestr(now)));
    clear sparm
    sparm.formulationType = 'margin';
    sparm.budget = i;
    sparm.p = size(D,1);
    sparm.func = @(A)(myGamma*func(A));
    if(i==1)
      [x,sparm,state] = LPcuttingPlane(sparm);
    else
      state.budget = i;
      state.hardb = -i;
      [x,sparm,state] = LPcuttingPlane(sparm,state);
    end
    save(fname,'x','sparm','state');
  else
    load(fname);
  end
  [s,xSortInd] = sort(x,'descend');
  valMargin(i) = sparm.func(xSortInd(1:i));
  valMarginLP(i) = sparm.bestobjective;
  sparms(i).sparm = sparm;
  xs(i).x = x;
end

sparms = [];
xs = [];
valSlackmargin = [];
valSlackmarginLP = [];
budgetsLP = 1:budget;
for i=1:budget
  fname = [saveDir sprintf('slackmarginResult%d.mat',i)];
  if(exist(fname,'file')~=2 || useSavedMargin==false)
    disp(sprintf(' running slack-margin rescaling with budget %d %s',i,datestr(now)));
    clear sparm
    sparm.formulationType = 'slackmargin';
    sparm.budget = i;
    sparm.p = size(D,1);
    sparm.func = @(A)(myGamma*func(A));
    if(i==1)
      [x,sparm,state] = LPcuttingPlane(sparm);
    else
      state.budget = i;
      state.hardb = -i;
      [x,sparm,state] = LPcuttingPlane(sparm,state);
    end
    save(fname,'x','sparm','state');
  else
    load(fname);
  end
  [s,xSortInd] = sort(x,'descend');
  valSlackmargin(i) = sparm.func(xSortInd(1:i));
  valSlackmarginLP(i) = sparm.bestobjective;
  sparms(i).sparm = sparm;
  xs(i).x = x;
end

figure()
hold on
  plot(valMarginLP(1:4)./myGamma,'color','cyan');
plot(valMargin./myGamma,'color','black');
plot(valSlackmarginLP./myGamma,'c-.');
plot(valSlackmargin./myGamma,'k-.');
plot(valNemhauser,'color','green');
plot(valNemhauser./(1-1/exp(1)),'g-.');
plot(mean(valRand'),'color','blue');
legend('margin LP bound', 'margin LP discretized', 'LP bound','LP discretized','Nemhauser Wolsey greedy','Greedy offline bound','random')
patch( [1:budget, budget:-1:1],[mean(valRand')-std(valRand')./sqrt(numRandRepititions), mean(valRand(end:-1:1,:)')+std(valRand(end:-1:1,:)'./sqrt(numRandRepititions))],[1,0.5,0.5])
plot(mean(valRand'),'color','blue');
hold off

