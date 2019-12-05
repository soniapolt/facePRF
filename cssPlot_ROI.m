% loads & plots some basics re: invPRF fits (css model)
% SP 10/25/17
% _indiv code allows us to plot individual ROIs in multiple plots
% 10/8/18 edit: made more flexible across expriment with different numbers
% of fit conditions

clear all; close all;

subjs = prfSubjs;%{'SP' 'DF' 'EM' 'TH' 'MG' 'JG'};%{'george'};%;
expt = 'fixPRF';%'nhp';%'

minR2 = 'r2-20';%'perc-50';          % cutoff for vox selection
ROIs= standardROIs;%['hV4' standardROIs('face')]; %{'ML' 'PL'};%{};%('face')
% manual set of baseCond + compConds
% [baseCond, compCond], more flexibly defined
base = 2; comps = [2 1];

saveFig = 1;

whichStim = 'outline';%'edge';%'binary';%'internal';%
whichModel = 'kayCSS';%'inflipCSSn';%'kayCSS';%'tempCSSn';%'%cssShift';%
fitSuffix = '';

hems = {'lh' 'rh'};%{''};%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotBasics = 1;
%plotBaselines = 0;
plotXY =1;
plotSize = 1;

fontSize = 12; titleSize = 14;

% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

for r = 1:length(ROIs)
if length(subjNum) == 1 bFits = subj(subjNum).roi(ROInum(r)).fits; 
elseif length(subjNum)== 0 error('Missing this subject in prfSet!');
else bFits =roi(ROInum(r)).fits; end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  create supertitle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist('task','var')
        titleText = [task ' task, ']; else titleText = [];end
    titleText = [whichModel ' ' titleText hemText(hems) ' ' ROIs{r} ', Subjs: ' strTogether(subjs) ...
    ' (' num2str(length(bFits(1).vox)) ' voxels R^2 > ' num2str(minR2) '), ' whichStim];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure 1:
    % 1) base coverage
    % 2) size x eccen
    % 3) stim sample
    % 4) r2 across conditions
    % 5) comp 1 coverage
    % 6) comp 1 eccen vs. base
    % 7) comp 1 size vs. base
    % 8) comp 1 r2 vs. base
    % 5) ...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotBasics == 1
       if onLaptop niceFig([.1 .1 .8 .8],fontSize); else 
           niceFig([0 .2 .6 .25*(size(comps,1)+1)],fontSize); end
        numPlots = [1+size(comps,1) 4]; pl = 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1) coverage, baseCond
        subplot(numPlots(1),numPlots(2),pl)
        % plotCoverage(vox,color,leg,ppd,res,plotSize,alphaGain,sampleVox,centerMass,plotCirc)
        plotCoverage(bFits(base).vox,condColors(base),bFits(base).cond,roi(1).fits(1).ppd,roi(1).fits(1).res,0);
        pl = pl+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2) size by eccen in all conditions
        subplot(numPlots(1),numPlots(2),pl)
        
        for c = 1:length(bFits)
            s(c) = scatter([bFits(c).vox.eccen],[bFits(c).vox.size],30,condColors(c)); hold on;
        end
        l=lsline; l=fliplr(l); for c = 1:length(bFits) set(l(c),'Color', condColors(c),'LineWidth',2); end
        
        xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (2*SD/sqrt(N)) (dva)'],'FontSize',fontSize); title('Size by Eccentricity','fontSize',titleSize);%ylim([0 2.5]);
        axis square; g = legend([s(1:length(bFits))],{bFits.cond}); set(g,'box','off','FontSize',fontSize,'location','best');
        pl = pl+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3) stim example
        subplot(numPlots(1),numPlots(2),pl)
        load([dirOf(pwd) 'cssFit/' roi(1).fits(1).stim]);
        imshow(condAvg(:,:,1)); title(['Sample pRF Stim Coding (' roi(1).fits(1).cond ')'],'fontSize',titleSize);
        pl = pl+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4) plot of r2 values - all, before the cutoff is implemented
        subplot(numPlots(1),numPlots(2),pl)
        allr2s = [bFits(1).vox.r2];
        for c = 2:length(bFits)
            allr2s = [allr2s; [bFits(c).vox.r2]];
        end
        boxplot(allr2s',{bFits(1:end).cond},'colors',[condColors(1:length(bFits))]); hold on; %hline(minR2,'k:','cutoff');
        title('R^{2} Values for plotted fits','fontSize',titleSize);
        axis square;
        pl = pl+1;
        
        for cc = 1:size(comps,1)
            baseCond = comps(cc,1);
            c = comps(cc,2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % col 1) coverage
            subplot(numPlots(1),numPlots(2),pl)
            plotCoverage(bFits(c).vox,condColors(c),bFits(c).cond,roi(1).fits(1).ppd,roi(1).fits(1).res,0);
            pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % col 4) eccen vs basecond
            subplot(numPlots(1),numPlots(2),pl)
            scatterCent([bFits(baseCond).vox.eccen],[bFits(c).vox.eccen],condColors(c),...
                bFits(baseCond).cond,bFits(c).cond,'Eccen (dva)',fontSize,1,1);
            pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % col 3) size vs basecond
            subplot(numPlots(1),numPlots(2),pl)
            
            scatterCent([bFits(baseCond).vox.size],[bFits(c).vox.size],condColors(c),...
                bFits(baseCond).cond,bFits(c).cond,'Size (2*SD/sqrt(N)) (dva)',fontSize,1,1);
            pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % col 5) r2 vs basecond
            subplot(numPlots(1),numPlots(2),pl)
            scatterCent([bFits(baseCond).vox.r2],[bFits(c).vox.r2],condColors(c),...
                bFits(baseCond).cond,bFits(c).cond,'R^{2}',fontSize,1,1);
            pl = pl+1;
        end
        
        
        superTitle(titleText,titleSize,.05);
        if saveFig 
        txt = [hemText(hems) '_' whichModel  '_' whichStim '_basics'];
        niceSave([dirOf(pwd) 'figures/' expt '/ROIplots/' ROIs{r}  '/'],txt,[],subjs); % just save pngs, since these can be generated pretty quickly
        end
        
        %     for c = 1:2
        %         [~,p,~,stats]=ttest([bFits(c).vox.size],[bFits(end).vox.size]);
        %         fprintf('%s vs. %s Size: p = %1.3f, t=%1.3f\n',bFits(c).cond,bFits(end).cond,p,stats.tstat);
        %         [~,p,~,stats]=ttest([bFits(c).vox.eccen],[bFits(end).vox.eccen]);
        %         fprintf('%s vs. %s Eccen: p = %1.3f, t=%1.3f\n',bFits(c).cond,bFits(end).cond,p,stats.tstat);
        %     end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(strfind(whichModel,'Shift')) && plotBaselines
        if onLaptop niceFig([.1 .1 .8 .8],fontSize); else niceFig([.3 .3 .6 .5],fontSize);end
        numPlots = [1 size(comps,1)+1];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure 2:
        % 1) shift boxplot - all conds
        % 2) shift scatter - cond1 vs cond3
        % 3) shift scatter - cond2 vs cond3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1) plot of r2 values - all, before the cutoff is implemented
        subplot(numPlots(1),numPlots(2),1)
        
        baselines = [];
        for c = 1:length(bFits)
            baselines = [baselines;bFits(c).vox.baseline]; end
        
        boxplot(baselines',{bFits(1:end).cond},'colors',[condColors(1:length(bFits))]);hline(0,'k-');
        title('Baseline Estim. Across Conditions','fontSize',fontSize,1,1);
        axis square;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2&3) shifts vs baseCond
        for cc = 1:size(comps,1)
            baseCond = comps(cc,1);
            c = comps(cc,2);
            subplot(numPlots(1),numPlots(2),1+cc)
            scatterCent([bFits(baseCond).vox.baseline],[bFits(c).vox.baseline],condColors(c),...
                bFits(baseCond).cond,bFits(c).cond,'Baselines (no outlier trim)',fontSize,1,1);
        end
        
        superTitle(titleText,titleSize,.97)
        if saveFig 
        txt = [hemText(hems) '_' whichModel '_' whichStim '_baselines'];
        niceSave([dirOf(pwd) 'figures/' expt '/ROIplots/' ROIs{r}  '/'],txt,[],subjs); % just save pngs, since these can be generated pretty quickly
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotXY == 1
        if onLaptop niceFig([.1 .1 .8 .8],fontSize); else niceFig([.1 .2 .8 .5*size(comps,1)],fontSize); end
        numPlots = [size(comps,1) 5]; pl = 1;
        set(gcf,'renderer','Painters');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure 3:
        % 1) cond 1 shifts by angle
        % 2) cond 1 shifts by in/out
        % 3) cond 1 shifts histogram
        % 4) cond 1 shifts by eccen
        % 5) cond 1 shifts by r2
        % 1) ...
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for cc = 1:size(comps,1)
            baseCond = comps(cc,1);
            c = comps(cc,2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1) shifts colored by angle
            subplot(numPlots(1),numPlots(2),pl)
            [angles{c}, lengths{c}] = plotXYshift(bFits(baseCond).vox,bFits(c).vox,roi(1).fits(1).ppd,roi(1).fits(1).res);
            
            ylabel([bFits(baseCond).cond ' to ' bFits(c).cond],'FontSize',titleSize+8,'FontWeight','bold');
            set(get(gca,'YLabel'),'Visible','on');
            pl = pl+1;
            
            % [left bottom width height]
            pos = get(gca,'Position');
            legSz = .08/size(comps,1);
            axes('Position',[pos(1) pos(2)+pos(4)-2*legSz legSz legSz]);
            drawcolorbarcircular(cmapang,1);
            
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % 2) shifts colored by in/out
%             subplot(numPlots(1),numPlots(2),pl)
%             [angles{c}, lengths{c}, sign{c}] = plotXYshift2(bFits(baseCond).vox,bFits(c).vox,roi(1).fits(1).ppd,roi(1).fits(1).res);
%             title(['{\color{red}Toward\color{black}-\color{blue}Away From\color{black} Center}'],'FontSize',titleSize,'interpreter','tex');
%             set(get(gca,'title'),'Visible','on');
%             pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % 3) magnitude of X change histogram 
            subplot(numPlots(1),numPlots(2),pl)
            shift = [];
            for v = 1:length(bFits(1).vox)
                shift(end+1) = bFits(baseCond).vox(v).XYdeg(1)- bFits(c).vox(v).XYdeg(1);
            end
            niceHist(shift,condColors(c),1);
            xlabel('Magnitude of X Shift (dva)'); ylabel('Count');
            title(['X Shift Magnitude'],'FontSize',titleSize);
            pl = pl+1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 3) magnitude of Y change histogram 
            subplot(numPlots(1),numPlots(2),pl)
            shift = [];
            for v = 1:length(bFits(1).vox)
                shift(end+1) = bFits(baseCond).vox(v).XYdeg(2)- bFits(c).vox(v).XYdeg(2);
            end
            niceHist(shift,condColors(c),1);
            xlabel('Magnitude of Y Shift (dva)'); ylabel('Count');
            title(['Y Shift Magnitude'],'FontSize',titleSize);
            pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 4) shifts by eccen
            subplot(numPlots(1),numPlots(2),pl)
            scatterCent([bFits(baseCond).vox.eccen],shift,condColors(c),...
                'Eccen (dva)','Magnitude of Shift (dva)','',fontSize,0,2);
            title(['Y Shift by Eccen'],'FontSize',titleSize);
            
            pl = pl+1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 5) shifts by r2
            subplot(numPlots(1),numPlots(2),pl);
            deltaR2 = [bFits(2).vox.r2] - [bFits(1).vox.r2];
            scatterCent(deltaR2,shift,[.5 .5 .5],...
                'Delta R^2','Shift','',fontSize,0,2);
            title(['Y Shift by R^2'],'FontSize',titleSize);
            pl = pl+1;
        end
        if saveFig 
        txt = [hemText(hems) '_' whichModel '_' whichStim  '_XYshift'];
        niceSave([dirOf(pwd) 'figures/' expt '/ROIplots/' ROIs{r}  '/'],txt,[],subjs); % just save pngs, since these can be generated pretty quickly
        end
        superTitle(['XY Shift: ' titleText],titleSize,.97);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotSize == 1
        if onLaptop niceFig([.1 .1 .8 .8],fontSize); else niceFig([.15 .2 .85 .5*size(comps,1)],fontSize); end
        numPlots = [size(comps,1) 4]; pl = 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure 4:
        % 1) cond 1 circle plots of prf size changes - cond1 vs cond3
        % 2) cond 1 size change magnitude
        % 3) cond 1 size change by eccen
        % 4) cond 1 size change by r2
        % 5) ...
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for cc = 1:size(comps,1)
            baseCond = comps(cc,1);
            c = comps(cc,2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1) dot plot of size changes
            subplot(numPlots(1),numPlots(2),pl)
            [sizeCh{c}] = plotSizeChange(bFits(baseCond).vox,bFits(c).vox,roi(1).fits(1).ppd,roi(1).fits(1).res);
            title(['{\color{red}Bigger\color{black}-\color{blue}Smaller \color{black}}'],'FontSize',titleSize,'interpreter','tex');
            set(get(gca,'title'),'Visible','on');
            
            ylabel([bFits(baseCond).cond ' to ' bFits(c).cond],'FontSize',titleSize+8,'FontWeight','bold');
            set(get(gca,'YLabel'),'Visible','on');
            pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 2) histogram of size changes
            subplot(numPlots(1),numPlots(2),pl)
            niceHist(sizeCh{c}/roi(1).fits(1).ppd,condColors(c),1)
            xlabel('Size Change (dva)'); ylabel('Count');
            title(['Size Change Magnitude'],'FontSize',titleSize);
            pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 3) size ch vs eccen
            subplot(numPlots(1),numPlots(2),pl)
            scatterCent([bFits(baseCond).vox.eccen],(sizeCh{c}/roi(1).fits(1).ppd),condColors(c),...
                'Eccen (dva)','Size Change (dva)','',fontSize,0,2);
            title(['Size Changeby Eccen'],'FontSize',titleSize);
            pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 4) size change by r2
            subplot(numPlots(1),numPlots(2),pl);
            deltaR2 = [bFits(2).vox.r2] - [bFits(1).vox.r2];
            scatterCent(deltaR2,[sizeCh{c}/roi(1).fits(1).ppd],[.5 .5 .5],...
                'Delta R^2','Size Change (dva)','',fontSize,0,2);
            title(['Size Change by R^2'],'FontSize',titleSize);
            pl = pl+1;
            
        end
        superTitle(['Size Changes: ' titleText],titleSize,.97);
        if saveFig 
        txt = [hemText(hems) '_' whichModel '_' whichStim  '_sizeChange'];
        niceSave([dirOf(pwd) 'figures/' expt '/ROIplots/' ROIs{r}  '/'],txt,[],subjs); % just save pngs, since these can be generated pretty quickly
        end
        
    end
    %if saveFig close all; end
end % ROIs

if onLaptop playSound; end