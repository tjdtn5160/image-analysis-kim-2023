function [] = LearnAndEvaluate(dataMatrix,labels,labelDist,plotFigures,Color)

% learns a model from data and performs cross validation to evaluate the

% model

% dataMatrix is the set of data we learn from. Each row is a sample

% lables is

% label dist is the distribution type of the labels. e,g, labelDist='binomial';



%Example data for run

% dataMatrix = randn(1000,2);

% noise = randn(1000,1);

% y = 0.7 .*dataMatrix(:,1) + 0.5 .* dataMatrix(:,2) + noise;

% labels = y > 0;

% labelDist='binomial';

% plotFigures=1;

if nargin<5
    Color=[0 0 0]
end


[numSamples,~]=size(dataMatrix);

randOrder = randperm(numSamples);

numDivisions=10;

divisionSize = floor(numSamples/numDivisions);



for divisionIndex=1:numDivisions
    
    if (divisionIndex==1)
        
        firstIndexTrain=1;
        
        lastIndexTrain=divisionSize;
        
        firstIndexTest=divisionSize+1;
        
        lastIndexTest=numSamples;
        
        
        
        trainIndices = randOrder(firstIndexTrain:lastIndexTrain);
        
        testIndices = randOrder(firstIndexTest:lastIndexTest);
        
        
        
    elseif (divisionIndex==numDivisions)
        
        firstIndexTrain=numSamples-divisionSize+1;
        
        lastIndexTrain=numSamples;
        
        firstIndexTest=1;
        
        lastIndexTest=numSamples-divisionSize;
        
        
        
        trainIndices = randOrder(firstIndexTrain:lastIndexTrain);
        
        testIndices = randOrder(firstIndexTest:lastIndexTest);
        
        
        
    else
        
        firstIndexTrain = ((divisionIndex-1)*divisionSize) +1;
        
        lastIndexTrain = (divisionIndex)*divisionSize;
        
        
        
        trainIndices = randOrder(firstIndexTrain:lastIndexTrain);
        
        testIndices = randOrder([ 1:(firstIndexTrain-1), (lastIndexTrain+1):numSamples]);
        
        
        
    end
    
    
    
    %train
    
    resultModel{divisionIndex}=glmfit(dataMatrix(trainIndices,:),labels(trainIndices),labelDist);
    
    %test
    
    modelTrainPrediction{divisionIndex} = glmval(resultModel{divisionIndex}, dataMatrix(trainIndices,:), 'logit');
    
    modelTestPrediction{divisionIndex} = glmval(resultModel{divisionIndex}, dataMatrix(testIndices,:), 'logit');
    
    realLabelsTrain{divisionIndex} = labels(trainIndices);
    
    realLabelsTest{divisionIndex} = labels(testIndices);
    
end

thresholds=[1:-0.01:0];

numThresholds=length(thresholds);

tprateTrain=nan(1,numThresholds);

fprateTrain=nan(1,numThresholds);

tprateTest=nan(1,numThresholds);

fprateTest=nan(1,numThresholds);



for threshInd=1:numThresholds
    
    thresh = thresholds(threshInd);
    
    tprateTrainDiv=nan(1,numDivisions);
    
    tprateTestDiv=nan(1,numDivisions);
    
    fprateTrainDiv=nan(1,numDivisions);
    
    fprateTestDiv=nan(1,numDivisions);
    
    
    
    for divisionIndex=1:numDivisions
        
        [tprateTrainDiv(divisionIndex),fprateTrainDiv(divisionIndex)]=calcROCVals(realLabelsTrain{divisionIndex},modelTrainPrediction{divisionIndex},thresh);
        
        [tprateTestDiv(divisionIndex),fprateTestDiv(divisionIndex)]=calcROCVals(realLabelsTest{divisionIndex},modelTestPrediction{divisionIndex},thresh);
        
    end
    
    tprateTrain(threshInd)=mean(tprateTrainDiv);
    
    fprateTrain(threshInd)=mean(fprateTrainDiv);
    
    tprateTest(threshInd)=mean(tprateTestDiv);
    
    fprateTest(threshInd)=mean(fprateTestDiv);
    
end



if (plotFigures~=0)
    
    %figure;
    
    %line([0 1],[0 1],'color','k');
    
    hold on;
    
    %plot(fprateTrain,tprateTrain,'LineWidth',2);
        
    plot(fprateTest,tprateTest,'color',Color,'LineWidth',2);
        
    %legend('control line','Train','Test');
    
    %title('ROC Curve');
    
    xlabel('False Positive Rate');
    
    ylabel('True Positive Rate');
    
    
    
    set(gca, 'FontName', 'Arial');
    
    set(gca, 'FontSize', 18);
    
end



end



function [tprate,fprate]=calcROCVals(trueLabels,modelPredictions,threshold)

numSamples = length(trueLabels);

pos = sum(trueLabels);

neg = numSamples-pos;

posPredictions = (modelPredictions>threshold);

trueP = sum(trueLabels.*posPredictions);

falseP = sum(posPredictions)-trueP;

%TP rate = TP/P (sensitivity)

tprate=trueP/pos;

%FP rate = FP/N (1-specificity)

fprate=(falseP/neg);

end