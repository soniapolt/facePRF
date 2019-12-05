% loads & plots some basics re: invPRF fits of the cssModel (positive pRF
% constraint) & the cssShift model (added a baseline term)
%
% SP 4/25/18
% this version plots ROIs as rows, params as columns, subjects as separate graphs
%
% version 2 implements standard pRF trimming, which we can hopefully use
% across different visualizations of the fits

clear all; close all;
addUtils;

subjs = {'JG'};
sessNums = [1 1 1];
task = 'fix';
expt = 'invPRF3';

%session = sessions{1};

whichStim = 'photo';
whichCond = 3; % 1 = inverted, 2 = mis-aligned, 3 = upright faces
allModels = {'negCSS','cssShift','cssNegShift','kayCSS'};
comp = [4,2];

minR2 = 50;
ROIs= {'V1','V2' 'V3' 'hV4'};%{'IOG_faces' 'pFus_faces' 'mFus_faces'};%
plotPars = {'r2' 'x' 'y' 'sd' 'n' 'size' };
parTitles = {'Estimated R^{2}' 'X Estim. (pix)' 'Y Estim. (pix)' 'SD Estim. (pix)' 'Exponent Estim' 'Size [2xSD/sqrt(N)] (dva)'};
compModels = {allModels{comp(1)},allModels{comp(2)}};
hems = {'rh' 'lh'};

sp = [length(ROIs) length(plotPars)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fontSize = 14; titleSize = 18;

for ss = 1:length(subjs)
    switch expt
        case 'invPRF3'
            [session, ~] = invPRF3_sessions(subjs{ss},sessNums(ss),task);
        case 'fixPRF'
            [session, ~] = fixPRF_sessions(subjs{ss},sessNums(ss));
    end
    for r = 1:length(ROIs)
        for m = 1:length(compModels)
            for h = 1:length(hems)
                
                thisROI = vpnlROI([hems{h} '_' ROIs{r}],session(1:2));
                
                clear fits;
                [~,fitFile] = fitsDirs(dirOf(pwd),expt,session,whichStim,thisROI,compModels{m});
                
                fprintf('Fits come from: %s\n',fitFile);
                
                load(fitFile);
                
                if h ==1
                    model(m).roi(r) = fits(whichCond); % do this only initially
                else
                    model(m).roi(r).vox = [model(m).roi(r).vox fits(whichCond).vox]; % do this on every iter
                    model(m).roi(r).ROIname = {model(m).roi(r).ROIname;fits(whichCond).ROIname}; % do this on every iter
                end
                
            end % hems
        end % modelFiles
    end % ROIs
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % translate fit params into more meaningful terms for invPRF expt              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for r = 1:length(ROIs)
        for m = 1:length(compModels)
            
            model(m).roi(r).vox = readPRFs(model(m).roi(r).vox,fits(1).ppd,fits(1).res);
            % for this particular plotting
            for v = 1:length(model(m).roi(r).vox)
                model(m).roi(r).vox(v).x = model(m).roi(r).vox(v).params(2);
                model(m).roi(r).vox(v).y = model(m).roi(r).vox(v).params(1);
                model(m).roi(r).vox(v).sd = model(m).roi(r).vox(v).params(3);
                model(m).roi(r).vox(v).n = model(m).roi(r).vox(v).params(5);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% which voxels are we plotting? trim edge values, R2 cutoff                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:length(ROIs)
    plotVox = 1:length(model(1).roi(r).vox);
    
    for m = 1:length(compModels)
        for v = 1:length(model(m).roi(r).vox)
            if model(m).roi(r).vox(v).r2 < minR2 ...
                    || trimPRFs(model(m).roi(r).vox(v).params,model(m).roi(r).vox(v).betas,fits(1).ppd,fits(1).res)
                plotVox(find(plotVox==v)) = []; end
        end
    end
    
    for m = 1:length(compModels)
        model(m).roi(r).vox = model(m).roi(r).vox(plotVox);
    end
end

niceFig([.1 .1 .9 .9],fontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot # 1: estimated params for css & shift models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(ROIs)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) plot of r2 values - post cutoff
    for p = 1:length(plotPars)
        subplot(sp(1),sp(2),(r-1)*sp(2)+p)
        
        for m = 1:length(compModels)
            eval(['in{m} = [model(m).roi(r).vox.' plotPars{p} '];']);
        end
        
        scatterCent(in{1},in{2},condColors(r),...
            compModels{1},{ROIs{r};[num2str(length(model(1).roi(r).vox)) ' Voxels'];' '; compModels{2}},...
            parTitles{p},fontSize);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% titles

titleText = [];
if length(hems)>1
    titleText = ['bilateral ROIs ' titleText];
else titleText = [hems{1} ' ' titleText]; end
titleText = [titleText ', Session = ' session];

titleText = [titleText ' (R^{2} thresh: ' num2str(minR2) ')'];
superTitle(titleText,titleSize,.95);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot #2: distribution of shifts for the shiftCSS model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:length(compModels)
    if ~isempty(strfind(compModels{m},'Shift'))
        niceFig([.1 .1 1 .7],fontSize);
        
        for r = 1:length(ROIs)
            %subplot(2,length(ROIs),r);
            
            shifts = [model(m).roi(r).vox.baseline];
            
            subplot(1,length(ROIs),r)
            
            thresh = shifts(find(shifts<nanmean(shifts)+3*nanstd(shifts)));
            thresh = thresh(find(thresh>nanmean(shifts)-3*nanstd(shifts)));
            
            hold on; hist(thresh,100); axis square;
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor',condColors(r),'FaceAlpha',.5,'EdgeAlpha',0);
            set(gca,'box','off','color','none');
            title([num2str(length(thresh)) ' Voxels']);
            
            vv = vline(nanmean(thresh),'k:',num2str(nanmean(thresh))); set(vv,'Color',condColors(r),'LineWidth',2);
            v = vline(0,'k:');
            xlabel('Outliers Removed');
            
            titleText = [model(m).roi(1).whichModel ' (' compModels{m} '), Subj. ' session(1:2)];
            superTitle(titleText,titleSize);
        end
        
    end % final end of ROIs
end
%end % subjects