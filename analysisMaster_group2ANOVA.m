function analysisMaster_group2ANOVA(which,exptSettings,ROIs,factNames,textThis,tails)
% extract the values our ANOVA needs from the either GLM or PSC group
% results
% 
% clear all; close all;
% which = 'GLM';
% exptSettings = 'LGNfigureAttn';
% ROIs = {'LGN' 'V1'};
% factNames = {'eye cond' 'orientation'};
% textThis = 1;
% 
load(['exptParams/' exptSettings '.mat']);

%tails = 1; % are our ttests one or two-tailed

%
% saveThis = 0;
% textThis = 0;

resultsDir = [fMRIdir '/' expt '/groupResults/'];
if ~exist(resultsDir) mkdir(resultsDir); end

% save output to a text file
if textThis == 1
    diaryName =  [resultsDir '/' expt '_' which '_ANOVA_group.txt'];
    if exist(diaryName)>0 delete(diaryName); end
    eval(['diary ' diaryName]);
    
    % to put at the top of our text file
    exptSubjs
    runTime = datestr(now)
end


ANOVA = struct;
ANOVA.subjects = exptSubjs;
for r = 1:length(ROIs)
    ANOVA(r).ROI = ROIs{r};
    groupData = [];
    for n = 1:length(exptSubjs)
        eval(['load ' resultsDir 'group_' ROIs{r} '_' which '.mat']);
        
        groupData = [groupData.subjMeans]; % mean for each condition for each subject
    end
    anovaData = reshape(groupData(:,1:length(condNames)),1,length(condNames)*length(exptSubjs))';
    %%%% ANOVA TIME!
    % Parameters:
    %    Y          dependent variable (numeric) in a column vector
    %    S          grouping variable for SUBJECT
    %    F1         grouping variable for factor #1
    %    F2         grouping variable for factor #2
    %    FACTNAMES  a cell array w/ two char arrays: {'factor1', 'factor2'}
    %
    %    Y should be a 1-d column vector with all of your data (numeric).
    %    The grouping variables should also be 1-d numeric, each with same
    %    length as Y. Each entry in each of the grouping vectors indicates the
    %    level # (or subject #) of the corresponding entry in Y.
    %    stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
    
    subjInd= repmat([1:length(exptSubjs)],1,length(condNames))';
    fact1 = kron([1:2],ones(1,length(exptSubjs)*2))';
    fact2 = repmat(kron([1:2],ones(1,length(exptSubjs))),1,2)';
   
    ANOVA(r).results = rm_anova2(anovaData,subjInd,fact1,fact2,factNames);
    %[anovaData,subjInd,attention,salience]
    fprintf(ROIs{r})
    mean(groupData)
    ANOVA(r).results
 
 %areas
%diary OFF

% some additional ttests, relevant for the LGNfigure experiments
for c = 1:length(compNames)
fprintf([compNames{c}  ', ' ANOVA(r).ROI '\n']);
if tails == 2
    fprintf(['two-tailed tests:\n']);
    [~,P,~,STATS] = ttest(groupData(:,compNums(c,2))-groupData(:,compNums(c,1)),0,'tail','both');
else if tails == 1
    fprintf(['one-tailed tests:\n']);
    [~,P,~,STATS] = ttest(groupData(:,compNums(c,2))-groupData(:,compNums(c,1)),0,'tail','right');
end
fprintf(['p = ' num2str(P) ', tstat = ' num2str(STATS.tstat) ', df = ' num2str(STATS.df) '\n\n']);
end

%diary ON
end
end
