% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = {'SP' 'MG' 'JG' 'TH' 'EM' 'DF' };%
task = '';
expt = 'fixPRF';
noCenters = 0;

saveFig = 0;
convertDVA = 1; % convert X and Y measurements into DVA

minR2 = 20;          % cutoff for vox selection
ROIs= {'hV4' 'mFus_faces'};%standardROIs;%'V1' 'V2' 'V3' 'hV4'

whichStim = 'internal';%'photo';%'external';%
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%
whichM = 3; % 1 = mean, 2 = mode/peak, 3 = median
hems = {'rh' 'lh'};
fitSuffix = '';%'_orig';%

switch expt
    case 'fixPRF'
        baseCond = [2]; compConds = [1];
    case 'compPRF'
        baseCond = [3]; compConds = [1 2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins

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
for p = 1:length(roi(1).fits(1).parNames)
    for c = 1:length(compConds)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.2 .1 .8 .8]; end
        niceFig(figSize,fontSize,1);
        numPlots = [2 ceil(length(ROIs)/2)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for r = 1:length(ROIs)
        titleText = [whichModel ' ' roi(1).fits(1).parNames{p} ', Subj: '];
        titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
        
        
        subplot(numPlots(1),numPlots(2),r)
            % get param values for this condition
            pars = vertcat(roi(ROInum(r)).fits(compConds(c)).vox.params);
            cPars = pars(:,p)';
            
            % get param values for base condition
            pars = vertcat(roi(ROInum(r)).fits(baseCond).vox.params);
            bPars = pars(:,p)';
            if convertDVA && ~containsTxt(roi(1).fits(1).parNames{p},'exp') && ~containsTxt(roi(1).fits(1).parNames{p},'gain')
        % rescale some parameters so that they are in DVA units and
        % centered around zero (center of screen)
        if ~containsTxt(roi(1).fits(1).parNames{p},'sd') % don't re-center the SD
            cPars = roi(1).fits(1).res-cPars-roi(1).fits(1).res/2; bPars = roi(1).fits(1).res-bPars-roi(1).fits(1).res/2;
        end
        cPars = cPars./roi(1).fits(1).ppd; bPars = bPars./roi(1).fits(1).ppd;
        end
            
            scatterCent(bPars,cPars,condColors(compConds(c)),...
                roi(ROInum(r)).fits(baseCond).cond,roi(ROInum(r)).fits(compConds(c)).cond,[ROIs{r} ' (' num2str(length(roi(r).fits(1).vox)) ' vox)'],fontSize);
            
            if containsTxt(roi(1).fits(1).parNames{p},'gain') % cut off gain plots
            xlim([0 10]); ylim([0 10]); end
        
            %title(ROIs{r},'fontSize',titleSize,'interpreter','none','FontWeight','bold');
            %xlabel(roi(1).fits(1).parNames{p},'fontSize',titleSize);
            %xlabel(ROIs{r},'fontSize',titleSize,'interpreter','none','FontWeight','bold');
            
        end
        superTitle(titleText,titleSize,.025);
        
        if saveFig == 1
            if length(subjs) == 1 
                txt = ['subj' subjs{1}]; 
            else txt = ['groupN' num2str(length(subjs))]; end
            if exist('task','var')
                txt = [roi(1).fits(1).parNames{p} '_' task 'Task_' txt ];
            else txt = [roi(1).fits(1).parNames{p} '_' txt ];end
            if length(hems) == 1
                txt = [txt '_' hems{1}]; end
            if ~containsTxt(whichStim,'photo')
            txt = [whichStim '_' txt];  end
            txt = ['scatter_' whichModel '_' txt];
            niceSave([dirOf(pwd) 'figures/' expt '/params/'],txt); % just save pngs, since these can be generated pretty quickly
        end
    end
end
if onLaptop playSound; end