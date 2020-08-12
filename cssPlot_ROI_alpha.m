% loads & plots some basics re: invPRF fits (css model)
% SP 10/25/17
% _indiv code allows us to plot individual ROIs in multiple plots
% 10/8/18 edit: made more flexible across expriment with different numbers
% of fit conditions

clear all; close all;

subjs = prfSubjs;%{'SP' 'DF' 'EM' 'TH' 'MG' 'JG'};%{'george'};%;
expt = 'fixPRF';%'nhp';%'

minR2 = 'r2-20';%'perc-50';          % cutoff for vox selection
ROIs= {'mFus_faces'};%{'V1'};%standardROIs;%['hV4' standardROIs('face')]; %{'ML' 'PL'};%{};%('face')
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
plotXY =0;
plotSize = 0;

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
    if onLaptop niceFig([.1 .1 .8 .8],fontSize); else 
           niceFig([0 .2 .6 .25*(size(comps,1)+1)],fontSize); end
        numPlots = [1+size(comps,1) 4]; pl = 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1) coverage, baseCond
        subplot(numPlots(1),numPlots(2),pl)
        % plotCoverage(vox,color,leg,ppd,res,plotSize,alphaGain,sampleVox,centerMass,plotCirc)
        plotCoverage(bFits(base).vox,roiColors(ROIs{r}),bFits(base).cond,roi(1).fits(1).ppd,roi(1).fits(1).res,0);
        pl = pl+1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2) size by eccen in all conditions
        subplot(numPlots(1),numPlots(2),pl)
        
        for c = 1:length(bFits)
           s(c) = scatter([bFits(c).vox.eccen],[bFits(c).vox.size],5,[1 .5 1],'filled'); hold on;
        end
        l=lsline; l=fliplr(l); for c = 1:length(bFits) set(l(c),'Color', roiColors(ROIs{r})*(c*.5),'LineWidth',2);  end
        
        for c = 1:length(bFits)
         delete(s(c));
         alpha = ([bFits(c).vox.r2]-min([bFits(c).vox.r2]))./(max([bFits(c).vox.r2])-min([bFits(c).vox.r2]));
         s2{c} = scatterAlpha([bFits(c).vox.eccen],[bFits(c).vox.size],alpha,roiColors(ROIs{r})*(c*.5),5);  hold on;
        end    
        
        xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (2*SD/sqrt(N)) (dva)'],'FontSize',fontSize); title('Size by Eccentricity','fontSize',titleSize);%ylim([0 2.5]);
        axis square; set(gca,'TickDir','out');% g = legend([s2(1:length(bFits))],{bFits.cond}); set(g,'box','off','FontSize',fontSize,'location','best');
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
        axis square;set(gca,'TickDir','out');
        pl = pl+1;
        
        for cc = 1:size(comps,1)
            baseCond = comps(cc,1);
            c = comps(cc,2)
            alpha = ([bFits(c).vox.r2]-min([bFits(c).vox.r2]))./(max([bFits(c).vox.r2])-min([bFits(c).vox.r2]));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % col 1) coverage
            subplot(numPlots(1),numPlots(2),pl)
            plotCoverage(bFits(c).vox,roiColors(ROIs{r})*.5,bFits(c).cond,roi(1).fits(1).ppd,roi(1).fits(1).res,0);
            pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % col 4) gain vs basecond
            subplot(numPlots(1),numPlots(2),pl)
             hold on;
            % scatterCent(x,y,color,xlab,ylab,titleTx,fontSize,equalLims,line,alpha)
            scatterCent([bFits(baseCond).vox.gain],[bFits(c).vox.gain],black,...
                bFits(baseCond).cond,bFits(c).cond,'Gain',fontSize,1,1,alpha);
            xlim([0 5]); ylim([0 5])
            set(gca,'TickDir','out');
            pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % col 3) size vs basecond
            subplot(numPlots(1),numPlots(2),pl)
            
            scatterCent([bFits(baseCond).vox.size],[bFits(c).vox.size],black,...
                bFits(baseCond).cond,bFits(c).cond,'Size (2*SD/sqrt(N)) (dva)',fontSize,1,1,alpha);
            set(gca,'TickDir','out');
            pl = pl+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % col 5) r2 vs basecond
            subplot(numPlots(1),numPlots(2),pl)
            scatterCent([bFits(baseCond).vox.r2],[bFits(c).vox.r2],black,...
                bFits(baseCond).cond,bFits(c).cond,'R^{2}',fontSize,1,1,alpha);
            set(gca,'TickDir','out');
            pl = pl+1;
        end
        
        
        superTitle(titleText,titleSize,.05);
        if saveFig 
        txt = [hemText(hems) '_' whichModel  '_' whichStim '_' minR2 '_basics_alpha'];
        niceSave([dirOf(pwd) 'figures/' expt '/ROIplots/' ROIs{r}  '/'],txt,[],[],{'png' 'svg'}); % just save pngs, since these can be generated pretty quickly
        end
 
end % ROIs

if onLaptop playSound; end