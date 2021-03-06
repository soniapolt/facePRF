% loads & plots distributions of XY changes/anything else per each subject
%%% this is the manuscript fig 2 code

clear all; %close all;

subjs = prfSubjs;%{'MG' 'JG' 'TH' 'EM' 'DF' 'SP'};%
expt = 'fixPRF';

saveFig =1;
convertDVA = 1;
flipX = 1;

minR2 = 'r2-20';          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')]

whichStim = 'internal';%'outline';%'photo';%'internal';%'eyes';%
whichModel = 'kayCSS';%'inflipCSSn';%'flipCSSn';%'cssExpN';%'cssShift';%
whichM = 'mean';

%plotPar = 'X';
%parTitle = 'X estim';%'Gain Estim';
plotPars = {'X' 'Y' 'gain' 'size' 'r2'};%{'gain' 'r2' 'Y' 'X' };
parTitles = {'X Estim.' 'Y Estim.' 'Gain Estim' 'Size [Sigma/sqrt(N)] (dva)' 'r2'};%{'Gain Estim' 'Estimated R^{2}' 'Y Estim.' 'X Estim.' 'Size [Sigma/sqrt(N)] (dva)'};


hems = {'lh' 'rh'};
fitSuffix = '';%'_orig';%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
nBins = 20; % histogram bins

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
niceFig([.1 .4 .2*length(plotPars) .3]);

for p = 1:length(plotPars)
    plotPar = plotPars{p}; parTitle = parTitles{p};
    subplot(1,length(plotPars),p);
    
allPars = nan(2*length(ROIs),length(subjs)); n=1;
colors = [];
for r = 1:length(ROIs)

ROInum = cellNum(ROIs{r},info.ROIs);
  
        for c = [2 1]
            for s = 1:length(subjs)
                fits = subj(s).roi(ROInum).fits(c);
                try
                eval(['allPars(n,s) = nan' whichM '(getPar(plotPar,fits,1,flipX));']);
                catch allPars(n,s) = NaN; end % for missing values
            end
            allFactors{n} = [ROIs{r} '-' fits.cond(1:3)];
            n=n+1;
            if c == 1 mult = .5; else mult = 1; end
            colors = [roiColors(ROIs{r})*mult; colors];
        end
end

if containsTxt(plotPar,'gain') cutY = [0 9]; 
elseif containsTxt(plotPar,'r2') cutY=[0 100]; 
else cutY = []; end
set(gca,'TickDir','out'); 
% niceBoxplot(data,labels,plotMed,colors,cutY)
%niceBoxplotGrouped(allPars',ROIs,fliplr({roi(1).fits.cond}),0,colors,cutY,1);
niceBoxPlusGrouped(allPars',ROIs,fliplr({roi(1).fits.cond}),flipud(colors),cutY,1,0);       

%niceBoxplot2(allPars',allFactors,0,colors,cutY,1);
axis square;
title(plotPar); titleText = [expt ' (' hemText(hems) '), voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];


superTitle(titleText,titleSize,.025);
end

    if saveFig == 1
        txt = ['acrossROIs_parSummary_' whichStim '_' whichM '_' hemText(hems) '_flipX' num2str(flipX)];
        niceSave([dirOf(pwd) 'figures/' expt '/crossSubj/'],txt,[],[],{'svg' 'png'}); % just save pngs, since these can be generated pretty quickly
    end

if onLaptop playSound; end