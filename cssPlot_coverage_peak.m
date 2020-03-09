% creates and saves a bootstrap-ed coverage map for every individual
% subject, and then averages them together for the total map.
% hybrid of mrvista coverage generation (bootstrapping) and SP coverage
% generation (no spatial smoothing, other choices...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

subjs = prfSubjs;
expt = 'fixPRF';

minR2 = 'r2-20';%['perc-50'];%20;          % cutoff for vox selection
ROIs= ['V1' standardROIs('face')];%+');
sampleVox = 0; % to place centers of voxels. if zero, don't do this.

whichStim = 'outline';
whichModel = 'kayCSS';
plotPeak = 1; %  peak of each image
hems = {'lh' 'rh'};

boot.iters = 1000;
boot.vox = 0.8; % now implements this as a proportion of total voxels, not an absolute number
boot.method = 'mean';% 'max';%% 'mean' or 'max'
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

for r = ROInum
    niceFig([.02*r .02*r .8 .5]);
    % plot mean coverage
    [H,P,CI,STATS] = ttest(boot.roi(r).cond(1).area,boot.roi(r).cond(2).area);
    
    for c = 1:length(boot.roi(1).cond)
        subplot(1,2,abs(c-3));
        switch boot.method
            case 'mean'
                meanIm = squeeze(nanmean(boot.roi(r).cond(c).cov));
            case 'max' % additionally, threshold across participants
                meanIm = squeeze(nanmean(boot.roi(r).cond(c).cov));
                meanIm(find(meanIm<.5)) = 0;
                %meanIm = imgaussfilt(meanIm,boot.ppd/10);
        end
        
        plotCovIm(meanIm,boot.res,boot.ppd,1,1,[]); hold on;
        contourLines = [.5 .9];  hold on;
        [cc,h] = contour(flipud(meanIm),repmat(max(meanIm(:)),1,length(contourLines)).*contourLines,...
            'Color','w','LineWidth',1);
        h.LevelList= round(h.LevelList,3);  %rounds levels to 2nd decimal place
        clabel(cc,h,'LabelSpacing',2000,'Color','w');
        
        if plotPeak
            [peak,idx]=max(meanIm(:));
            [row,col]=ind2sub(size(meanIm), idx);
            hold on; scatter(col,size(meanIm,2)-row,100,'k*');
        end
        
        % plot some pRF centers
        %%%% sample a subset of voxels randomly
        %                 sv = randperm(length(roi(r).fits(c).vox));
        %                 if sampleVox>length(sv) sampleVox = length(sv); end
        %                 vox = roi(r).fits(c).vox(sv(1:sampleVox));
        %                 hold on;
        %                 col = colorinterpolate([1 1 1; roiColors(ROIs{find(ROInum==r)})],5,1); col = col(4,:); m =1;
        %                 col = 'w';
        %                 scatter([vox.Xdeg].*boot.ppd+boot.res/2,m*[vox.Ydeg].*boot.ppd+boot.res/2,centSize,col,'filled');
        %
        
        x = xlabel({'Mean FWHM:';sprintf('%.2f dva',mean(boot.roi(r).cond(c).area));...
            sprintf('SE (N=%d) = %.2f',length(subjNum),se(boot.roi(r).cond(c).area))});
        set(x,'visible','on');
        y= ylabel(roi(1).fits(c).cond,'FontSize',20); set(y,'visible','on');
        t = title({[ROIs{find(ROInum==r)}];sprintf('t(%d) = %.2f',STATS.df,STATS.tstat);sprintf('p = %.3f',P)});
        set(t,'visible','on');if H set(t,'Color',[0 .5 0]); end
    end
    superTitle(fileName(grpCov),14,.05);
    niceSave([raid 'invPRF/figures/fixPRF/coverage/centroid/'],[hemText(hems) '_' ROIs{find(ROInum==r)} '_scale' num2str(boot.scaleSubjs)],[],[],{'svg'});
end
if onLaptop playSound; end

% centroid calc
for rr = 1:length(ROInum)
    r = ROInum(rr);
    for c = 1:2
        for s = 1:length(subjNum)
            sIm = squeeze(boot.roi(r).cond(c).cov(s,:,:));
             [peak,idx]=max(sIm(:));
            [row,col]=ind2sub(size(sIm), idx);
            boot.roi(r).cond(c).centX(s) = col;
            boot.roi(r).cond(c).centY(s) = size(meanIm,2)-row;
        end
    end
    
    %%% centroid computation
    [h,p,~,stats] = ttest(boot.roi(r).cond(1).centX,boot.roi(r).cond(2).centX); if h==1 st = '***'; else st = ''; end
    fprintf('%s%s X ttest: t(%d) = %.3f, p = %.3f\n',st,ROIs{rr},stats.df, stats.tstat,p);
    [h,p,~,stats] = ttest(boot.roi(r).cond(1).centY,boot.roi(r).cond(2).centY);if h==1 st = '***'; else st = ''; end
    fprintf('%s%s Y ttest: t(%d) = %.3f, p = %.3f\n',st,ROIs{rr},stats.df, stats.tstat,p);
end

% param plot
plotPars = {'centX' 'centY' 'area'}; clear toPlot;
niceFig([.1 .1 .4 .4]);
for n = 1:length(plotPars)
   for c = 1:2
       for rr = 1:length(ROInum)
    r = ROInum(rr);

        for s = 1:length(subjNum)
          eval(['toPlot(rr,s) =  boot.roi(r).cond(c).' plotPars{n} '(s);']);
        end
       end
  subplot(length(plotPars),1,n);
  
      switch plotPars{n}
        case 'centX'
            toPlot = (toPlot-boot.res/2)/boot.ppd;
        case 'centY'
            toPlot = -(toPlot-boot.res/2)/boot.ppd;
        case 'area'
            toPlot = toPlot;
    end
  
    if c == 1 mult = .25; else mult = 1; end
  scatterErr([1:length(ROInum)]+.1*(c-1),nanmean(toPlot,2),se(toPlot')',roiColors(ROIs)*mult);hold on;
  title(plotPars{n});
   xticks([1:length(ROInum)]); xticklabels(ROIs); ylabel(['Peak of FWHM (dva)']);

        xlim([0 length(ROIs)+1]);
        l = hline(0,'k:');set(l,'LineWidth',1.5);hold on;
   end
end

niceSave([raid 'invPRF/figures/fixPRF/coverage/peak/'],[hemText(hems) '_heirarchy_scale' num2str(boot.scaleSubjs)],[],[],{'svg'});

