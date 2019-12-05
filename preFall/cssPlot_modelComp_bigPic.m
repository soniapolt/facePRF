% loads & plots some basics re: invPRF fits of the cssModel (positive pRF
% constraint) & the cssShift model (added a baseline term)
%
% SP 4/25/18
% this version plots ROIs as rows, params as columns, subjects as separate graphs
%
% version 2 implements standard pRF trimming, which we can hopefully use
% across different visualizations of the fits

clear all; close all;

subjs = {'DF' 'EM' 'TH' 'JG' 'MG' 'SP'};%;%;
task = '';
expt = 'fixPRF';
saveFig = 0;
convertDVA = 1;

whichModels = {'cssExpN' 'cssExpN' 'intempCSSn' 'inflipCSSn'};
fitSuffix = '';%'_orig';%
whichCond = 1;
allStims = {'photo','internal','internal','internal'};
fitsSuffix = ''; %'_orig';
compConds = [2 1];

minR2 = 50;
ROI=standardROIs(7);
whichM = 'median';
plotPars = {'r2' 'Y' 'size' 'gain'};
parTitles = {'Estimated R^{2}' 'Y Estim.' 'Size [2xSD/sqrt(N)] (dva)' 'Gain Estim'};
plotType = {'scatter' 'distr' 'distr' 'box'}; % 1 = boxplot, 2 = distr, 3 = scatter
hems = {'rh' 'lh'};
fitSuffix = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fontSize = 14; titleSize = 18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveFig && onLaptop figSize = [0 0 1 1]; else figSize = [.1 .1 .9 .9]; end
niceFig(figSize,fontSize);
numPlots = [length(allStims) length(plotPars)]; pl = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:length(allStims)
    whichModel = whichModels{t};
    load(pRFfile(dirOf(pwd),expt,minR2,allStims{t},whichModel,hems,fitSuffix,task));
    ROInum = cellNum(ROI,info.ROIs);
    subjNum = cellNum(subjs,info.subjs);
    if length(subjNum) == 1 roi = subj(subjNum).roi; end
    % assuming either one subject, or all subjects
    fits = roi(ROInum).fits;
    
    for p = 1:length(plotPars)
        
        subplot(numPlots(1),numPlots(2),pl)
        % get param values for this condition
        
        for cc = 1:length(compConds)
            c = compConds(cc);
            parNum = cellNum(plotPars{p},fits(1).parNames);
            if ~isempty(parNum)
                pars = vertcat(fits(c).vox.params);
                plPars{cc} = pars(:,parNum)';
            else
                eval(['plPars{cc} = [fits(c).vox.' plotPars{p} '];']);  end
            
%             if containsTxt(plotPars{p},'gain') && trimGains>0 % for the time being, only look at reasonable-ish gains
%                 z = plPars{c};
%                 z(find(z>trimGains))=NaN;
%                 plPars{c}=z;
%             end
            if convertDVA && containsTxt(plotPars{p},'Y') || containsTxt(plotPars{p},'X') || containsTxt(plotPars{p},'sd')
        % rescale some parameters so that they are in DVA units and
        % centered around zero (center of screen)
        if ~containsTxt(plotPars{p},'sd') % don't re-center the SD
            plPars{cc} = fits(1).res-plPars{cc}-roi(1).fits(1).res/2;
        end
        plPars{cc} = plPars{cc}./roi(1).fits(1).ppd;
        end
        end
        
        switch plotType{p}
            case 'box'
                if containsTxt(plotPars{p},'gain') 
                    niceBoxplot([plPars{1};plPars{2}]',{fits(compConds(1)).cond fits(compConds(2)).cond},1,[condColors(4);condColors(2)],[0 10]);
                    ylim([0 5]);
                elseif containsTxt(plotPars{p},'r2')                     
                    niceBoxplot([plPars{1};plPars{2}]',{fits(compConds(1)).cond fits(compConds(2)).cond},1,[condColors(4);condColors(2)],[minR2 100]);
                else
                niceBoxplot([plPars{1};plPars{2}]',{fits(compConds(1)).cond fits(compConds(2)).cond},1,[condColors(4);condColors(2)]);
                end
                if containsTxt(plotPars{p},'size') 
                   ylim([0 5]); 
                end
            case 'distr'
                nBins = 20;
                plotDistr(plPars,1,{fits(compConds(1)).cond fits(compConds(2)).cond},nBins,3,1);
                if containsTxt(plotPars{p},'size') 
                   xlim([0 10]); 
                else xlim([-5 5]); end
                
                [pv, od, effectsize] = permutationTest(plPars{1}, plPars{2}, 1000, ...
                     'plotresult', 0, 'showprogress', 0);
                 fprintf('%s Permutation Test, %s: p = %.6f, observed difference = %.2f\n ',...
                     whichModel, plotPars{p},pv,od);
            case 'scatter'
                hold on;
                scatterCent(plPars{1},plPars{2},condColors(c),...
                    fits(compConds(1)).cond,fits(compConds(2)).cond,[],fontSize,0,1);
                %hold on; l =lsline; l=fliplr(l);
                if strcmp(plotPars{p},'r2') xlim([minR2 100]); ylim([minR2 100]); end
                hold on; xl = get(gca,'xlim');
                hold on; plot(xl,xl,'k:');
        end
        %title([whichModels{t}  '-' allStims{t} ' stim' ' ' plotPars{p}],'fontSize',titleSize,'interpreter','none','FontWeight','bold'); 
        pl = pl+1;
        axis square;
    end
    
end

superTitle(sprintf('%s Model Comparisons, %s Params, R2 cutoff = %d',ROI{1},whichM,minR2),titleSize);

%subplotresize;
if saveFig == 1
    if length(subjs) == 1
        txt = ['subj' subjs{1}];
    else txt = ['groupN' num2str(length(subjs))]; end
    if length(hems) == 1
        txt = [txt '_' hems{1}]; end
    
    txt = [whichModel '_stimComp_bigPic_' txt];
    
    niceSave([dirOf(pwd) 'figures/' expt '/modelComp/'],txt); % just save pngs, since these can be generated pretty quickly
end

if onLaptop playSound; end