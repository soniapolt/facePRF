% creates and saves a bootstrap-ed coverage map for every individual
% subject, and then averages them together for the total map.
% hybrid of mrvista coverage generation (bootstrapping) and SP coverage
% generation (no spatial smoothing, other choices...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% subjs = prfSubjs; default now is all subjects
expt = 'fixPRF';

minR2 = 'r2-20';%['perc-50'];%20;          % cutoff for vox selection
ROIs= {'mFus_faces'};%standardROIs('face');%('face+');
sampleVox = 500; % to place centers of voxels. if zero, don't do this.

whichStim = 'outline';
whichModel = 'kayCSS';
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

exptDir = fullfile(raid,'invPRF',expt);
covPath = [exptDir '/coverage/']; checkDir(covPath);

bootOpts = [boot.method '_iters' num2str(boot.iters) '_vox' num2str(boot.vox) '_scale' num2str(boot.scaleSubjs) '.mat'];

grpCov = [covPath fileName(prfSet) '_' bootOpts];
load(grpCov);


tic
for r = 1:length(ROIs)
    ROInum = cellNum(ROIs{r},info.ROIs)';

    niceFig([.02*r .02*r .8 .5]);
 
    for c = 1:length(boot.roi(1).cond)
        subplot(1,2,abs(c-3));
        switch boot.method
            case 'mean'
                meanIm = squeeze(nanmean(boot.roi(ROInum).cond(c).cov));
            case 'max' % additionally, threshold across participants
                meanIm = squeeze(nanmean(boot.roi(ROInum).cond(c).cov));
                meanIm(find(meanIm<.5)) = 0;
                %meanIm = imgaussfilt(meanIm,boot.ppd/10);
        end
            [~,h,areaDeg] = plotContour(meanIm,boot.res,boot.ppd,0,0); 
            drawnow;
        
        % plot some pRF centers
                %%%% sample a subset of voxels randomly
                sv = randperm(length(roi(ROInum).fits(c).vox));
                if sampleVox>length(sv) sampleVox = length(sv); end
                
                vox = roi(ROInum).fits(c).vox(sv(1:sampleVox)); 
                hold on;
                % col = colorinterpolate([1 1 1; roiColors(ROIs{r})],5,1); col = col(4,:); m =1;
                col = 'k'; m = 1;
                scatter([vox.Xdeg].*boot.ppd+boot.res/2,m*[vox.Ydeg].*boot.ppd+boot.res/2,2,col,'filled');
            
        
        x = xlabel({'Mean FWHM:';sprintf('%.2f dva',mean(boot.roi(ROInum).cond(c).area));...
            sprintf('SE (N=%d) = %.2f',length(prfSubjs),se(boot.roi(ROInum).cond(c).area));...
            sprintf('Levels: %.2f %.2f %.2f %.2f %.2f %.2f',h.LevelList)});
        
        set(x,'visible','on');
        y= ylabel(roi(1).fits(c).cond,'FontSize',20); set(y,'visible','on');
        t = title(ROIs{r});
        set(t,'visible','on');
    end
    superTitle(fileName(grpCov),14,.05);
    niceSave([raid 'invPRF/figures/fixPRF/coverage/contour/'],['newer_' hemText(hems) '_' ROIs{r} '_vox' num2str(sampleVox) '_scale' num2str(boot.scaleSubjs)],[],[],{'png' 'svg'});
end
if onLaptop playSound; end