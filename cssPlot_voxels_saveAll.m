clear all; close all;

subj = 'SP';
expt = 'fixPRF';
task = '';

[session, numRuns] = vpnlSessions(expt,subj); % OPTIONAL: SESSNUM, TASK

whichStim = 'internal';
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%
fitSuffix = '';
noCenters = 0;
ROIs = standardROIs;
hems = {'rh'};
sortR2 = 1; % sort voxels by descending R2 or (0) grab random ones
minR2 = 20;

for r = 1:length(ROIs)
    for h = 1:length(hems)
        ROI = [hems{h} '_' ROIs{r}];
        
        [~, fitsName] = fitsDirs(dirOf(pwd),expt,session,whichStim,vpnlROI(ROI,subj),whichModel,fitSuffix);
        if exist(fitsName)>0
            load(fitsName);
            load([dirOf(pwd) 'cssFit/' fits(1).stim]); fontSize = 12;
            
            %%%%% are we choosing voxels by their R2 value?
            if sortR2 == 1
                % trim nan fits, which otherwise get sorted as >100
                for c = 1:length(fits) for n = 1:length(fits(c).vox)
                        if isnan(fits(c).vox(n).r2) fits(c).vox(n).r2=0; end
                    end
                end
                
                sortBy = [fits(1).vox.r2];
                for n = 2:length(fits)
                    sortBy = [sortBy; [fits(n).vox.r2]]; end
                sortBy = mean(sortBy);
                
                [~,plotVox] = sort(sortBy,2,'descend');
            else
                plotVox = 1:length(fits(1).vox);
            end
            
            %%%%% trim fits to plot
            % reasonable parameter estimates, at least one positive response per condition
            for v = 1:length(plotVox)
                for c = 1:length(fits)
                    if isfield(fits,'expN');
                        fits(c).vox = readPRFs(fits(c).vox,fits(1).ppd,fits(1).res,fits(1).expN);
                    else fits(c).vox = readPRFs(fits(c).vox,fits(1).ppd,fits(1).res,[]);end
                
                if  fits(c).vox(v).r2 < minR2 ...
                        || trimPRFs(fits(c).vox(v).params,fits(c).vox(v).betas,fits(1).ppd,fits(1).res);
                    plotVox(find(plotVox==v)) = []; end
            end
        end
        
        for p = 1:length(plotVox)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if onLaptop figSize = [.2 .2 .5 .5];
            else figSize = [.2 .2 .2 .2];
            end
            f = niceFig(figSize,fontSize);
            numPlots = [1 4];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%% plot betas and modelfit
            subplot(numPlots(1),numPlots(2),1);
            
            v = plotVox(p);
            for c = 1:length(fits)
                hbar(c) = niceBars([fits(c).condNums],[fits(c).vox(v).betas],[fits(c).vox(v).sems],condColors(c));
                plot(fits(c).condNums,fits(c).vox(v).modelfit,'b-','LineWidth',2);
            end
            
            %legend([hbar(1) hbar(2) hbar(3)],{fits.cond},'location','Best','FontSize',10,'box','off');
            
            xlabel(strTogether({fits.cond},10)); xlim([0 fits(end).condNums(end)+1]); set(gca,'xticklabel',[]);
            ylabel('Beta Estim.'); axis square;
            
            for tt = 1:length(fits)
                parText{tt} = ['R^{2}=' num2str(round(fits(tt).vox(v).r2)) ' ' fits(tt).cond ];
            end
            t = title(parText);
            
            set(gca,'box','off','FontName','Arial','FontWeight','normal','TickDir','out');
            
            %%%%%%%%% plot betas in space
            subplot(numPlots(1),numPlots(2),2:3);
            betaSpace = [];
            for c = 1:length(fits)
                betas = fits(c).vox(v).betas;
                if noCenters == 1;
                    betas = insertVal(NaN,betas,13);
                end
                betaSpace = [betaSpace reshape(betas,5,5)'];
            end
            
            plotInSpace(betaSpace,'BetaEstimates');
            
            %%%%%%%%% plot comparison of pRF locations
            subplot(numPlots(1),numPlots(2),4);
            
            for c = 1:length(fits)
                % PRF size is defined as S/sqrt(N).
                h(c) = plotCoverage(fits(c).vox(v),condColors(c),[],fits(1).ppd,fits(1).res);
                if c == length(fits) l=legend([h(1:length(fits))],{fits.cond});
                    set(l,'box','off','location','Best');
                    for tt = 1:length(fits)
                        parText{tt} = ['Params' num2str(tt) ': [' num2str(round(fits(tt).vox(v).params)) ']'];
                    end
                    t = title(parText);
                    set(t,'visible','on');
                end
            end
            
            titleText = [ROI ' Voxel #' num2str(v) ', ' whichStim ' stim, ' whichModel ' model'];
            
            titleText = [titleText ', Session = ' session];
            superTitle(titleText,12, .05)
            
            %%% save each voxel as fig
            if noCenters == 0;
                saveDir = ['voxPlots/' subj task '/' ROI '/'];
            else saveDir = ['voxPlots/noCenters/' subj task '/' ROI '/']; end
            checkDir([dirOf(pwd) 'figures/' expt '/' saveDir]);
            txt = [whichModel '_' num2str(p) '_vox' num2str(v) '_' whichStim];
            niceSave([dirOf(pwd) 'figures/' expt '/' saveDir],txt); % just save pngs, since these can be generated pretty quickly
            close(f);
        end
        else
            fprintf('Missing %s!\n',fitsName);
    end
end
end