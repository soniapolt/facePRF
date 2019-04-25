% loads & plots scatter plot of one par vs the other (e.g. gain vs r2)
% within condition rather than across conditions

clear all; close all;

subjs = {'SP' 'MG' 'JG' 'TH' 'EM' 'DF' };%
task = '';
expt = 'fixPRF';
noCenters = 0;

saveFig = 1;

compPars = {'gain' 'r2'}; %; 1:'Y'    2:'X'    3:'sd'    4:'gain'    5:'exp'
trimGains = 5;

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs;%{'IOG_faces' 'pFus_faces' 'mFus_faces'};%'V1' 'V2' 'V3' 'hV4'

whichStim = 'photo';%'internal';%
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
hems = {'rh' 'lh'};
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
    for c = 1:length(roi(1).fits)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.2 .1 .8 .8]; end
        niceFig(figSize,fontSize);
        numPlots = [2 ceil(length(ROIs)/2)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for r = 1:length(ROIs)
        titleText = [roi(1).fits(c).cond ' ' whichModel ' ' compPars{1} 'vs. ' compPars{2} ', Subj: '];
        titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
        
        
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
            
            scatterCent(plPars{1},plPars{2},condColors(4),...
                compPars{1},compPars{2},[ROIs{r} ' (' num2str(length(roi(r).fits(1).vox)) ' vox)'],fontSize,0,1);
            
        superTitle(titleText,titleSize,.025);
        end
        if saveFig == 1
            if length(subjs) == 1 
                txt = ['subj' subjs{1}]; 
            else txt = ['groupN' num2str(length(subjs))]; end
            if exist('task','var')
                txt = [compPars{1} 'V' compPars{2} '_' task 'Task_' txt ];
            else txt = [compPars{1} 'V' compPars{2} '_' txt ];end
            if length(hems) == 1
                txt = [txt '_' hems{1}]; end
            
            txt = ['scatter_' roi(1).fits(c).cond '_' whichModel '_' txt];
            niceSave([dirOf(pwd) 'figures/' expt '/params/'],txt); % just save pngs, since these can be generated pretty quickly
        end
        
    end
if onLaptop playSound; end