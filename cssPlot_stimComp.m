% loads & plots some basics re: invPRF fits of the cssModel (positive pRF
% constraint) & the cssShift model (added a baseline term)
%
% SP 4/25/18
% this version plots ROIs as rows, params as columns, subjects as separate graphs
%
% version 2 implements standard pRF trimming, which we can hopefully use
% across different visualizations of the fits

clear all; close all;

subjs = {'DF' 'EM' 'TH' 'JG' 'MG'};%;%'JG';
task = '';
expt = 'fixPRF';
saveFig = 1;

whichModel = 'kayCSS';
fitSuffix = '';%'_orig';%
whichCond = 1;
allStims = {'binary','photo','edge','internal'};
comp = [2,4];
fitsSuffix = ''; '_orig';

minR2 = 20;
ROIs=standardROIs('face');
plotPars = {'r2' 'x' 'y' 'sd' 'n' 'size' };
parTitles = {'Estimated R^{2}' 'X Estim. (pix)' 'Y Estim. (pix)' 'SD Estim. (pix)' 'Exponent Estim' 'Size [2xSD/sqrt(N)] (dva)'};
compStims = {allStims{comp(1)},allStims{comp(2)}};
hems = {'rh' 'lh'};
fitSuffix = '';%'_orig';%

sp = [length(ROIs) length(plotPars)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fontSize = 14; titleSize = 18;

for r = 1:length(ROIs)
    for m = 1:length(compStims)
        init = 1;
        for s = 1:length(subjs)
            [session, numRuns] = vpnlSessions(expt,subjs{s}); % OPTIONAL: SESSNUM, TASK

        for h = 1:length(hems)
            
            thisROI = [hems{h} '_' ROIs{r}];
            
            clear fits;
            [dataName, fitsName] = fitsDirs(dirOf(pwd),expt,session,compStims{m},vpnlROI(thisROI,session(1:2)),whichModel,fitsSuffix);
            fprintf('Fits come from: %s\n',fitsName);
            
            load(fitsName);
            
            if init ==1
                model(m).roi(r) = fits(whichCond); init = 0;% do this only initially
            else
                model(m).roi(r).vox = [model(m).roi(r).vox fits(whichCond).vox]; % do this on every iter
                model(m).roi(r).ROIname = {model(m).roi(r).ROIname;fits(whichCond).ROIname}; % do this on every iter
            end
            
        end % hems
        end
    end % modelFiles
end % ROIs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate fit params into more meaningful terms for invPRF expt              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:length(ROIs)
    for m = 1:length(compStims)
        
        model(m).roi(r).vox = readPRFs(model(m).roi(r).vox,fits(1).ppd,fits(1).res);
        
        for v = 1:length(model(m).roi(r).vox)
            % for this particular plotting
            model(m).roi(r).vox(v).x = model(m).roi(r).vox(v).params(2);
            model(m).roi(r).vox(v).y = model(m).roi(r).vox(v).params(1);
            model(m).roi(r).vox(v).sd = model(m).roi(r).vox(v).params(3);
            model(m).roi(r).vox(v).n = model(m).roi(r).vox(v).params(5);
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% which voxels are we plotting? trim edge values, R2 cutoff                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:length(ROIs)
    plotVox = 1:length(model(1).roi(r).vox);
    
    for m = 1:length(compStims)
        for v = 1:length(model(m).roi(r).vox)
            if model(m).roi(r).vox(v).r2 < minR2 ...
                    || trimPRFs(model(m).roi(r).vox(v).params,model(m).roi(r).vox(v).betas,fits(1).ppd,fits(1).res)
                plotVox(find(plotVox==v)) = []; end
        end
    end
    
    for m = 1:length(compStims)
        model(m).roi(r).vox = model(m).roi(r).vox(plotVox);
    end
end

figure;set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 .9 .9],'DefaultTextFontSize',fontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot # 1: estimated params for different stims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(ROIs)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) plot of r2 values - post cutoff
    for p = 1:length(plotPars)
        subplot(sp(1),sp(2),(r-1)*sp(2)+p)
        
        for m = 1:length(compStims)
            eval(['in{m} = [model(m).roi(r).vox.' plotPars{p} '];']);
        end
        
        if p == 1 text = {ROIs{r};[num2str(length(model(1).roi(r).vox)) ' Voxels'];' '; compStims{2}};
        else text = compStims{2}; end
        scatterCent(in{1},in{2},condColors(r),...
            compStims{1},text,...
            parTitles{p},fontSize);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% titles

titleText = [];
if length(hems)>1
    titleText = ['bilateral ROIs ' titleText];
else titleText = [hems{1} ' ' titleText]; end
if length(subjs) ==1
titleText = [titleText ', Session = ' session]; 
else titleText = [titleText ', group (N = ' num2str(length(subjs)) ')']; end
titleText = [titleText ' (R2 thresh: ' num2str(minR2) ')'];
superTitle(titleText,titleSize,.95);

if saveFig == 1
    txt = [compStims{1} 'v' compStims{2} '_' fits(whichCond).cond '_'];
    if length(hems) == 1
        txt = [txt '_' hems{1}]; end
    niceSave([dirOf(pwd) 'figures/' expt '/stimComp/'],txt,ROIs,subjs); % just save pngs, since these can be generated pretty quickly
end

