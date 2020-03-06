% loads & plots changes along the visual heirarchy for cssFit properties
% currently - across all subject in pRFset file

clear all; close all;

task = '';
expt = 'fixPRF';%'invPRF3';%

minR2 = 20;          % cutoff for vox selection
ROIs = standardROIs('face+');%{'hV4' 'IOG_faces'    'pFus_faces'    'mFus_faces'};
subjs = prfSubjs; % assuming either one subject, or all subjects
plotPar = 'gain';

% manual set of baseCond + compConds
% [baseCond, compCond], more flexibly defined
compConds = [2 1];%[1 2; 1 3; 2 3];
fitSuffix = '';
convertDVA = 1;

saveFig = 1;

whichStim = 'photo';
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%

hems = {'lh' 'rh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix,task));
subjNum = cellNum(subjs,info.subjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create supertitle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titleText = [expt ' (' hemText(hems) ') '];
titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1:
% 1) XY shift magnitude along visual heirarchy
% 2) size change along visual heirarchy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niceFig([.1 .1 .8 .55],fontSize);
allPars = nan(2*length(ROIs),length(roi(1).fits(1).vox)); %spNum = [1 length(ROIs)]; 
n = 1;
for r = 1:length(ROIs)
%    subplot(spNum(1),spNum(2),r)
    ROI = ROIs{r};
    ROInum = cellNum(ROI,info.ROIs);
    if length(subjNum) == 1 roi = subj(subjNum).roi; end
    % assuming either one subject, or all subjects
    fits = roi(ROInum).fits;
    
    for cc = 1:length(compConds)
        c = compConds(cc);
        parNum = cellNum(plotPar,fits(1).parNames);
        if ~isempty(parNum)
            pars = vertcat(fits(c).vox.params);
            plPars{cc} = pars(:,parNum)';
        else
            eval(['plPars{cc} = [fits(c).vox.' plotPar '];']);  end

        if convertDVA && containsTxt(plotPar,'Y') || containsTxt(plotPar,'X') || containsTxt(plotPar,'sd')
            % rescale some parameters so that they are in DVA units and
            % centered around zero (center of screen)
            if ~containsTxt(plotPar,'sd') % don't re-center the SD
                plPars{cc} = fits(1).res-plPars{cc}-roi(1).fits(1).res/2;
            end
            plPars{cc} = plPars{cc}./roi(1).fits(1).ppd;
        end
        allPars(n,1:length(plPars{cc})) = plPars{cc}; allFactors{n} = [ROI '-' fits(c).cond(1:3)];
        n = n+1;
    end
end
if containsTxt(plotPar,'gain') cutY = [0 10]; yl = [0 2];
elseif containsTxt(plotPar,'r2') cutY=[minR2 100]; 
else cutY = []; end
% niceBoxplot(data,labels,plotMed,colors,cutY)
niceBoxplot(allPars',allFactors,0,repmat([condColors(4);condColors(2)],length(ROIs),1),cutY);
title(plotPar); 

superTitle([titleText],titleSize,.97);

if saveFig == 1
    txt = [plotPar];
    if length(hems) == 1
        txt = [txt hems{1} '_']; end
    txt = [whichModel '_' txt];
    niceSave([dirOf(pwd) 'figures/' expt '/acrossROIs/' ],txt,ROIs,info.subjs); % just save pngs, since these can be generated pretty quickly
end

