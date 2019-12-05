% loads & plots scatter plot inv vs upr par difference vs. behavior

clear all; close all;

expt = 'fixPRF';

behavFile = [raid '/behavFIE/analysis/N9-prfRec2.mat'];
load(behavFile);
behText = 'FIE at Lower Left';
%behText = 'FIE at Center';
%behText = 'FIE Change from Center to Lower Left';

% set behav metric
switch behText
    case 'FIE at Lower Left'
behav = (grp(4).subjPC - grp(1).subjPC); % FIE at lower left
    case 'FIE at Center'
behav = (grp(5).subjPC - grp(2).subjPC); % FIE at center
    case 'FIE Change from Center to Lower Left'
        behav = (grp(4).subjPC - grp(1).subjPC)...
     - (grp(5).subjPC - grp(2).subjPC); % FIE across locations: 4-1, 5-2, 6-3
end

subjs = an.subjs;

saveFig = 1;

compPar = 'Y'; %'size'%; 1:'Y'    2:'X'    3:'sd'    4:'gain'    5:'exp'

minR2 = 'r2-20';['perc-50'];%;          % cutoff for vox selection
ROIs= standardROIs('face');%'mFus_faces';

whichStim = 'outline';%'photo';%'internal';%
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
hems = {'rh' 'lh'};
fitSuffix = '';
trimGains = 5; % if we're looking at gains?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix,''));
subjNum = cellNum(subjs,info.subjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: distribution of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.1 .1 .9 .9]; end
niceFig(figSize,fontSize); sp = [2, ceil(length(ROIs)/2)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        titleText = [ whichModel ' ' compPar ' vs. Behavior, Subj: '];
        titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];

        for r = 1:length(ROIs)
             ROI= ROIs{r};ROInum = cellNum(ROI,info.ROIs);
            subplot(sp(1),sp(2),r); title(ROI);
        % get param values for this condition
        for ss = 1:length(subjs)
            s = subjNum(ss);
            parNum = cellNum(compPar,subj(s).roi(1).fits(1).parNames);
            if ~isempty(parNum)
                pars1 = vertcat(subj(s).roi(ROInum).fits(1).vox.params);
                pars2 = vertcat(subj(s).roi(ROInum).fits(2).vox.params);
                sPars{ss} = pars2(:,parNum)'-pars1(:,parNum)';
            else
                eval(['sPars{ss} = [[subj(s).roi(ROInum).fits(2).vox.' compPar ']-[subj(s).roi(ROInum).fits(1).vox.' compPar ']];']);  end
            
            if containsTxt(compPar,'gain') && trimGains>0 % for the time being, only look at reasonable-ish gains
                z = sPars{ss};
                z(find(z>trimGains))=NaN;
                sPars{ss}=z;
            end
            medPars(ss) = median(sPars{ss})/roi(1).fits(1).ppd;
        end
        
        hold on;
        s = scatter(behav,medPars,30,condColors(randi(10)));
        text(behav+.01,medPars+.01,subjs);
        [rho,pval] = corr(behav',medPars')
        xlabel(['Behavior: ' behText]); ylabel(['Delta ' compPar]); title([ROI ' (R^2 = ' num2str(rho*rho) ', p= ' num2str(pval) ')']);
        hold on; l =lsline; l=fliplr(l);
        if strcmp(compPar,'r2') xlim([minR2 100]);end
        end
 
    superTitle(titleText,titleSize,.025);
       
        
% if saveFig == 1
%     if length(subjs) == 1
%         txt = ['subj' subjs{1}];
%     else txt = ['groupN' num2str(length(subjs))]; end
%     if exist('task','var')
%         txt = [compPars{1} 'V' compPars{2} '_' task 'Task_' txt ];
%     else txt = [compPars{1} 'V' compPars{2} '_' txt ];end
%     if length(hems) == 1
%         txt = [txt '_' hems{1}]; end
%     
%     txt = [whichModel '_' whichStim '_' txt];
%     if length(plotConds)==1
%         txt = [roi(1).fits(plotConds).cond '_' txt]; end
%     
%     niceSave([dirOf(pwd) 'figures/' expt '/paramsComp/'],txt); % just save pngs, since these can be generated pretty quickly
% end

if onLaptop playSound; end