% creates and saves a bootstrap-ed coverage map for every individual
% subject, and then averages them together for the total map.
% hybrid of mrvista coverage generation (bootstrapping) and SP coverage
% generation (no spatial smoothing, other choices...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this is the main coverage function as of mar 10 2020! adding support of
%%% max/mean, binary plotting

clear all; close all;

subjs = prfSubjs;
expt = 'fixPRF';

minR2 = 'r2-20';%['perc-50'];%20;          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')];%+');
sampleVox = 0; % to place centers of voxels. if zero, don't do this.
contourLabels = 0;
plotContours = 1;

whichStim = 'internal';%'outline';%;
whichModel = 'kayCSS';
plotCent = 1; % weighted centroid of each image
hems = {'rh' 'lh'};

%%% how were the bootstraps generated?
boot.iters = 1000;
boot.vox = 0.8; % now implements this as a proportion of total voxels, not an absolute number
boot.method = 'binary';%'mean';%'max';%'binary';% % 'mean' or 'max'
boot.scaleSubjs = 0; % rescale each individual's coverage to [0 1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prfSet = pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems);
load(prfSet);

ppd = roi(1).fits(1).ppd; res = roi(1).fits(1).res;

ROInum = cellNum(ROIs,info.ROIs)';
subjNum = cellNum(subjs,info.subjs);

exptDir = fullfile(raid,'invPRF',expt);
covPath = [exptDir '/coverage/']; checkDir(covPath);

bootOpts = [boot.method '_iters' num2str(boot.iters) '_vox' num2str(boot.vox) '_scale' num2str(boot.scaleSubjs) '.mat'];


tic
grpCov = [covPath fileName(prfSet) '_' bootOpts];
load(grpCov);
niceFig([.02 .02 .4 1]);
for r = ROInum
    
    % plot mean coverage
    [H,P,CI,STATS] = ttest(boot.roi(r).cond(1).area,boot.roi(r).cond(2).area);
    
    for c = 1:length(boot.roi(1).cond)
        rr = find(ROInum==r);
        subplot(length(ROIs),2,2*(rr-1)+abs(c-3));
        meanIm = squeeze(nanmean(boot.roi(r).cond(c).cov));
        if strcmp(boot.method,'max') % an additional thresholding for our max images
            meanIm(find(meanIm<.5)) = 0;
        end
        
        [~,cX(r), cY(r)] = plotCovIm(meanIm,boot.res,boot.ppd,1,0,[]); hold on;
        if plotContours
            contourLines = [.5 .5];  hold on;
            [cc,h] = contour(flipud(meanIm),repmat(max(meanIm(:)),1,length(contourLines)).*contourLines,...
                'Color','w','LineWidth',1);
        end
        if contourLabels
            h.LevelList= round(h.LevelList,3);  %rounds levels to 2nd decimal place
            clabel(cc,h,'LabelSpacing',2000,'Color','w');
        end
        
        %%% sample a subset of voxels
        if sampleVox > 0
                        [~,sv] = sort([roi(r).fits(c).vox.r2],'descend');
                        %sv = randperm(length(roi(r).fits(c).vox));
                        if sampleVox>length(sv) sampleVox = length(sv); end
                        vox = roi(r).fits(c).vox(sv(1:sampleVox));
                        hold on;
%                         col = colorinterpolate([1 1 1; roiColors(ROIs{find(ROInum==r)})],5,1); col = col(4,:); m =1;
%                         col = 'w';
                        scatter([vox.Xdeg].*boot.ppd+boot.res/2,-[vox.Ydeg].*boot.ppd+boot.res/2,1,white,'filled');
        end
        
        x = xlabel({'Mean FWHM:';sprintf('%.2f dva',mean(boot.roi(r).cond(c).area));...
            sprintf('SE (N=%d) = %.2f',length(subjNum),se(boot.roi(r).cond(c).area))});
        set(x,'visible','on');
        %y= ylabel(roi(1).fits(c).cond,'FontSize',20); set(y,'visible','on');
        t = ylabel({roi(1).fits(c).cond;[ROIs{find(ROInum==r)}];sprintf('t(%d) = %.2f',STATS.df,STATS.tstat);sprintf('p = %.3f',P)});
        set(t,'visible','on');if H set(t,'Color',[0 .5 0]); end
    end
    superTitle(fileName(grpCov),14,.05);
end

% centroid calc
for rr = 1:length(ROInum)
    r = ROInum(rr);
    for c = 1:2
        for s = 1:length(subjNum)
            [boot.roi(r).cond(c).centX(s),cY] = FWHMcentroid(squeeze(boot.roi(r).cond(c).cov(s,:,:)),.5);
            boot.roi(r).cond(c).centY(s) = boot.res+cY;
        end
        if sampleVox == 0
        subplot(length(ROIs),2,2*(rr-1)+abs(c-3));
        hold on; plot([boot.roi(r).cond(c).centX],[boot.roi(r).cond(c).centY],'w.','MarkerSize',3);
        hold on; plot(mean([boot.roi(r).cond(c).centX]),mean([boot.roi(r).cond(c).centY]),'r*','MarkerSize',5);
        end
    end
    
    %%% centroid computation
    [h,p,~,stats] = ttest(boot.roi(r).cond(1).centX,boot.roi(r).cond(2).centX); if h==1 st = '***'; else st = ''; end
    fprintf('%s%s X ttest: t(%d) = %.3f, p = %.3f\n',st,ROIs{rr},stats.df, stats.tstat,p);
    [h,p,~,stats] = ttest(boot.roi(r).cond(1).centY,boot.roi(r).cond(2).centY);if h==1 st = '***'; else st = ''; end
    fprintf('%s%s Y ttest: t(%d) = %.3f, p = %.3f\n',st,ROIs{rr},stats.df, stats.tstat,p);
end

niceSave([raid 'invPRF/figures/fixPRF/coverage/centroid/'],[hemText(hems) '_summary_'  whichStim '_' boot.method '_scale' num2str(boot.scaleSubjs) '_clabels' num2str(contourLabels)],[],[],{'png' 'svg'});

if onLaptop playSound; end


% param plot
plotPars = {'centX' 'centY' 'area'}; clear toPlot;
plotType = 'boxplus'; 'box';%'dots';
niceFig([.1 .1 .2 .15*length(plotPars)]);


for n = 1:length(plotPars)
    subplot(length(plotPars),1,n);
    if containsTxt(plotType,'dots')
                
                for c = 1:2
        for rr = 1:length(ROInum)
            r = ROInum(rr);
            
            for s = 1:length(subjNum)
                eval(['toPlot(rr,s) =  boot.roi(r).cond(c).' plotPars{n} '(s);']);
            end
        end
        
        switch plotPars{n}
            case 'centX'
                toPlot = (toPlot-boot.res/2)/boot.ppd;
            case 'centY'
                toPlot = (boot.res/2-toPlot)/boot.ppd;
            case 'area'
                toPlot = toPlot;
        end
        
        if c == 1 mult = .25; else mult = 1; end
        
        scatterErr([1:length(ROInum)]+.1*(c-1),nanmean(toPlot,2),se(toPlot')',roiColors(ROIs)*mult);hold on;
                end
         xticks([1:length(ROInum)]); xticklabels(ROIs);        
    else
            col = 1;
            for rr = 1:length(ROInum)
            r = ROInum(rr);
            for c = [2 1]
                for s = 1:length(subjNum)
            eval(['toPlot(col,s) =  boot.roi(r).cond(c).' plotPars{n} '(s);']);
                end
                 col = col+1;
            end
            end
            switch plotPars{n}
            case 'centX'
                toPlot = (toPlot-boot.res/2)/boot.ppd;
            case 'centY'
                toPlot = (boot.res/2-toPlot)/boot.ppd;
            case 'area'
                toPlot = toPlot;
            end
                % niceBoxplot2(data,labels,plotMed,dotColors,cutY,plotData)
              % niceBoxplot2(toPlot',Expand(ROIs,2,1),0,flipud(Expand(roiColors(ROIs),1,2)).*repmat([.25 1],3,length(ROIs))',[],1)  ;
            % xticks([1:2*length(ROInum)]); xticklabels(Expand(ROIs,2,1));  
            if containsTxt(plotType,'plus')
           niceBoxPlusGrouped(toPlot',ROIs,fliplr({roi(1).fits.cond}),Expand(roiColors(ROIs),1,2).*repmat([1 .25],3,length(ROIs))',[],1,0);       

            else
            niceBoxplotGrouped(toPlot',ROIs,fliplr({roi(1).fits.cond}),0,flipud(Expand(roiColors(ROIs),1,2)).*repmat([.25 1],3,length(ROIs))',[],1);       
            end
            end
         title(plotPars{n});
         ylabel(['Centroid of FWHM (dva)']);

         %set(gca,'tickdir','out');
       % l = hline(0,'k:');set(l,'LineWidth',1.5);hold on;
   
end

niceSave([raid 'invPRF/figures/fixPRF/coverage/centroid/'],[hemText(hems) '_' whichStim '_' boot.method '_barplot_scale' num2str(boot.scaleSubjs)],[],[],{'svg' 'png'});
