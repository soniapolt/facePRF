% loads & plots distributions of XY changes/anything else per each subject

clear all; close all;

subjs = prfSubjs;%{'MG' 'JG' 'TH' 'EM' 'DF' 'SP'};%
expt = 'fixPRF';

saveFig =1;
convertDVA = 1;

minR2 = 'r2-20';          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')]

whichStim = 'outline';%'photo';%'internal';%'eyes';%
whichModel = 'kayCSS';%'inflipCSSn';%'flipCSSn';%'cssExpN';%'cssShift';%
whichM = 'median';

plotPar = 'X';
parTitle = 'X estim';%'Gain Estim';

hems = {'rh' 'lh'};
fitSuffix = '';%'_orig';%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins

% now we load in the data from both hemispheres, and threshold across
niceFig([.4 .4 .4 .5]);


for h = 1:length(hems)
allPars = nan(2*length(ROIs),length(subjs)); n=1;
colors = [];
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems(h),fitSuffix));
for r = 1:length(ROIs)

ROInum = cellNum(ROIs{r},info.ROIs);
  
        for c = [2 1]
            for s = 1:length(subjs)
                fits = subj(s).roi(ROInum).fits(c);
                try
                allPars(n,s) = nanmedian(getPar(plotPar,fits,1));
                catch allPars(n,s) = NaN; end % for missing values
            end
            allFactors{n} = [ROIs{r} '-' fits.cond(1:3)];
            n=n+1;
            if c == 1 mult = .5; else mult = 1; end
            colors = [roiColors(ROIs{r})*mult; colors];
        end
end

if containsTxt(plotPar,'gain') cutY = [0 3]; 
elseif containsTxt(plotPar,'r2') cutY=[minR2 100]; 
else cutY = []; end

% niceBoxplot(data,labels,plotMed,colors,cutY)
niceBoxplot(allPars',allFactors,0,colors,cutY,1); hold on;
title(plotPar); titleText = [expt ' (' hemText(hems) '), voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];
set(gca,'TickDir','out');ylim([-5 5]);pbaspect([1 2 1]);
end
%superTitle(titleText,titleSize,.025);
    
    if saveFig == 1
        txt = ['acrossROIs_' plotPar '_twoHems'];
        niceSave([dirOf(pwd) 'figures/' expt '/crossSubj/'],txt,[],[],{'svg' 'png'}); % just save pngs, since these can be generated pretty quickly
    end

if onLaptop playSound; end