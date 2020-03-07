% loads & plots scatter plot of one par vs the other (e.g. gain vs r2)
% within condition rather than across conditions

clear all; close all;

subjs = prfSubjs;
task = '';
expt = 'fixPRF';

saveFig = 1;

compPars = {'gain' 'size'}; %; 1:'Y'    2:'X'    3:'sd'    4:'gain'    5:'exp'
trimGains = 5;
plotConds = [1 2]; % allows us to just vis one condition at a time

minR2 = 50;          % cutoff for vox selection
ROIs= {'V1' 'hV4' 'pFus_faces' 'mFus_faces'};%'V1' 'V2' 'V3' 'hV4'

whichStim = 'outline';%'photo';%'internal';%
whichModel = 'kayCSS';%tempCSSn';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
hems = {'rh' };
fitSuffix = '';%'_orig';%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix,task));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

if length(subjNum) == 1 roi = subj(subjNum).roi; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: distribution of parameters for this ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFig && onLaptop figSize = [.1 .1 .6 .6]; else figSize = [.1 .1 .9 .9]; end
niceFig(figSize,fontSize);
numPlots = [1 length(ROIs)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(ROIs)
    for c = plotConds
        
        titleText = [compPars{1} ' vs. ' compPars{2}];
        titleText = [titleText ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
        
        
        subplot(numPlots(1),numPlots(2),r)
        % get param values for this condition
        
        for p = 1:2
            parNum = cellNum(compPars{p},roi(1).fits(1).parNames);
            if ~isempty(parNum)
                pars = vertcat(roi(ROInum(r)).fits(c).vox.params);
                plPars{p} = pars(:,parNum)';
            else
                eval(['plPars{p} = [roi(ROInum(r)).fits(c).vox.' compPars{p} '];']);  end
            
            if containsTxt(compPars{p},'gain') && trimGains>0 % for the time being, only look at reasonable-ish gains
                z = plPars{p};
                z(find(z>trimGains))=NaN;
                plPars{p}=z;
            end
        end
        
        hold on;
        pl(c) = scatterCent(plPars{1},plPars{2},condColors(c),...
            compPars{1},compPars{2},[ROIs{r} ' (' num2str(length(roi(ROInum(r)).fits(1).vox)) ' vox)'],fontSize,0,0);
        alpha(.2);
        hold on; l =lsline; l=fliplr(l);
        if strcmp(compPars{1},'r2') xlim([minR2 100]);end
        
    end
    g = legend(pl(plotConds),{roi(1).fits(:).cond});
    for c = plotConds set(l(c),'Color', condColors(c),'LineWidth',2); end
    superTitle(titleText,titleSize,.025);
end
if saveFig == 1
    
    txt = [compPars{1} '-vs-' compPars{2} '_r2-' num2str(minR2) '_' ];
    if length(hems) == 1
        txt = [txt '_' hems{1}]; end
    
    txt = [txt whichModel '_' whichStim];
    if length(plotConds)==1
        txt = [txt '_'  roi(1).fits(plotConds).cond]; end
    
    niceSave([dirOf(pwd) 'figures/' expt '/paramComp/'],txt,[],[],{'svg'}); % just save pngs, since these can be generated pretty quickly
end

if onLaptop playSound; end