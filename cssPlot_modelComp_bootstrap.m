% plots a model comparison of the bootstrap results - currently set up to
% just work on the delta bootstrapping
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

expt = 'fixPRF';
% assumes all subjects and all ROIs...
% to add later - bootstrapping within subjects

minR2 = ['r2-20'];          % cutoff for vox selection

whichModels = {'kayCSS' 'kayCSS','kayCSS','kayCSS'};
modelStims = {'outline','photo', 'internal','eyes'};

hems = {'lh' 'rh'};
plotPars = {'gain' 'Y' 'size'};
saveFig = 1;
plotROIs = ['V1' standardROIs('face')];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
titleTxt = [];
for mm = 1:length(whichModels)
    whichModel = whichModels{mm};
    whichStim = modelStims{mm};
    
    prfSet = pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems);
    load(prfSet);
    if ~isfield(roi(1),'boot')
    error(sprintf('Missing bootstrapping for %s %s %s!\n',whichModel,whichStim,minR2)); end

%%% grab relevant ROIs for this model
ROInum = cellNum(plotROIs,info.ROIs);
roi = roi(ROInum);

%%% shade the different model estimates
colorSpace = linspace(.3,1,length(whichModels));


    if mm == 1 fig = niceFig([.1 .1 .4 .4],14); end
    for b = 1:length(plotPars)
        bb = cellNum(plotPars{b},{roi(1).boot.parName});
        % niceFig([.1 .1 .5 .5],14);
        subplot(length(plotPars),1,b);
        title(plotPars{b});hold on;
        
        for r = 1:length(plotROIs)
            % gather relevant data
            plotM(r,b) = roi(r).boot(bb).median;
            plotCI(r,b) = roi(r).boot(bb).CI(2) - roi(r).boot(bb).median;
        end
        scatterErr([1:length(plotROIs)]+.1*(mm-1),plotM(:,b),plotCI(:,b),roiColors(plotROIs)*colorSpace(mm));hold on;
    end
   titleTxt = [titleTxt '   \color[rgb]{' num2str(colorSpace(mm)*condColors(1)) '} ' whichModel '-' whichStim '-' minR2];
end

% figure updates - don't need to do with each model
for b = 1:length(plotPars)
        % niceFig([.1 .1 .5 .5],14);
        subplot(length(plotPars),1,b); xticks([1:length(plotROIs)]);
        xticklabels(plotROIs); ylabel(['delta ' plotPars{b} ' (Upr - Inv)']);
        
        xlim([0 length(plotROIs)+1]);
        l = hline(0,'k:');set(l,'LineWidth',1.5);hold on;
end

t = superTitle(titleTxt,14,.05);set(t,'Interpreter','Tex');

if saveFig
    niceSave([raid 'invPRF/figures/fixPRF/modelComp/'],['bootComp_' minR2],[],info.subjs,{'png' 'svg'}); 
end