% creates and saves a bootstrap estimates for relevant params on a pRF-sets
% struct. saves it back in the PRF-sets .mat file. i've spun this out of
% any plotting functions because it will be relatively time intensive
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

expt = 'fixPRF';
% assumes all subjects and all ROIs...
% to add later - bootstrapping within subjects

minR2 = ['r2-20'];          % cutoff for vox selection

whichStim = 'outline';%'internal';
whichModel = 'kayCSS';%'cssExpN';%'inflipCSSn';%%
hems = {'lh' 'rh'};
bootPars = {'gain' 'r2' 'Y' 'X' 'size' 'eccen'};
plotIt = 1;
saveFig = 1;
doBoot = 1; % re-run bootstrapping even if it has already been done with this prfSet
plotROIs = ['V1' 'hV4' standardROIs('face')];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prfSet = pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems);
load(prfSet);

if ~isfield(roi(1),'boot') || doBoot
    
    info.boot = 'delta bootstrapping';
    info.bootPars = bootPars;
    info.numBoot = 1000;
    info.AminusB = [2 1]; % which conditions are we subtracting for the bootstrap estimate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % bootstrap o'clock                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for r = 1:length(info.ROIs)
        fprintf(info.ROIs{r})
        tic
        for b = 1:length(bootPars)
            pars = getPar(bootPars{b},roi(r).fits(info.AminusB(1)),1)-getPar(bootPars{b},roi(r).fits(info.AminusB(2)),1);
            
            roi(r).boot(b).parName = bootPars{b};
            [roi(r).boot(b).median ...
                roi(r).boot(b).CI ...
                roi(r).boot(b).dist] = bootstrapCI(pars,[],info.numBoot,[]);
            
            if plotIt
                if b == 1 niceFig([.1 .1 .8 .8]); pl = 1;
                    superTitle([info.ROIs{r} ' bootstrap, numIter = ' num2str(info.numBoot)],14,.05);end
                subplot(2, round(length(bootPars)/2),pl);
                niceHist(roi(r).boot(b).dist,condColors(r,1),1);
                t = title(sprintf('%s %s %s Bootstrap %s: Median = %.2f, CI = [%.2f %.2f]',whichModel, whichStim, minR2,bootPars{b},roi(r).boot(b).median,...
                    roi(r).boot(b).CI(1),roi(r).boot(b).CI(2)));
                if roi(r).boot(b).CI(1)<0 && roi(r).boot(b).CI(2)>0 set(t,'Color',[168 0 0]/255); end % mark CIs that include zero
                pl = pl+1;
                if b == length(bootPars) && saveFig
                    niceSave([dirOf(pwd) 'figures/fixPRF/bootPars/' info.ROIs{r} '/' ],['deltaBoot_' whichModel '_' whichStim '_' minR2 '_' info.ROIs{r}],[],info.subjs);
                end
            end
        end
        toc
    end
    save(prfSet,'info','subj','roi'); % save pRFset with bootstrapping results
end

%%% make summary figure across ROIs
ROInum = cellNum(plotROIs,info.ROIs);
roi = roi(ROInum);

niceFig([.1 .1 .9 .9],14);
for b = 1:length(bootPars)
    % niceFig([.1 .1 .5 .5],14);
    subplot(length(bootPars),1,b);
    title(bootPars{b});hold on;
    for r = 1:length(plotROIs)
        % gather relevant data
        plotM(r,b) = roi(r).boot(b).median;
        plotCI(r,b) = roi(r).boot(b).CI(2) - roi(r).boot(b).median;
    end
    
    scatterErr(1:length(plotROIs),plotM(:,b),plotCI(:,b),roiColors(plotROIs));hold on;
    
    xlim([0 length(plotROIs)+1]); xticks(1:length(plotROIs)); xticklabels(plotROIs); 
    y = ylabel({['delta \bf' bootPars{b}]; '\rm(Upr - Inv)'},'FontSize',16); set(y,'Interpreter','Tex');
    l = hline(0,'k:');set(l,'LineWidth',1.5);hold on;
    superTitle(sprintf('Model: %s, Stim: %s, Hems: %s R2 cutoff: %s',whichModel,whichStim,hemText(hems),minR2),16,.05);
end




%Save the coverage info
% saveFile = fullfile(savePath, [hems{h} '_coverage_data_ve', num2str(ve_cutoff*100)]);
