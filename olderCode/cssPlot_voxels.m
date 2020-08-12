clear all; close all;


subj = 'SP';
task = '';
expt = 'fixPRF';
noCenters = 0;
[session, numRuns] = vpnlSessions(expt,subj); % OPTIONAL: SESSNUM, TASK


whichStim = 'photo';
whichModel = 'kayCSS';%'cssExpN';%'cssShift';%

ROI = 'rh_mFus_faces';%'rh_mFus_faces';%
sortR2 = 1; % sort voxels by descending R2 or (0) grab random ones
minR2 = 10;
numPlot = 1;%[1:10];

[~, fitsName] = fitsDirs(dirOf(pwd),expt,session,whichStim,vpnlROI(ROI,subj),whichModel);
load(fitsName);
load([dirOf(pwd) 'cssFit/' fits(1).stim]); fontSize = 12;

%%%%% are we choosing voxels by their R2 value?
if sortR2 == 1
    sortBy = [];
    % trim nan fits, which otherwise get sorted as >100
    for c = 1:length(fits) for n = 1:length(fits(c).vox)
            if isnan(fits(c).vox(n).r2) fits(c).vox(n).r2=0; end
        end
        sortBy = [sortBy; [fits(c).vox.r2]];
    end
    sortBy = mean(sortBy);
    [~,plotVox] = sort(sortBy,2,'descend');
else
    plotVox = Shuffle(1:length(fits(1).vox));
end
plotVox = 62;

%%%%% trim fits to plot
% reasonable parameter estimates, at least one positive response per condition
for v = 1:length(plotVox)
    for c = 1:length(fits)
        if  fits(c).vox(v).r2 < minR2 ...
                || trimPRFs(fits(c).vox(v).params,fits(c).vox(v).betas,fits(1).ppd,fits(1).res);
            plotVox(find(plotVox==v)) = []; end
    end
end

% load someVox.mat; % if we've previously chosen some voxels to plot


if length(numPlot)>length(plotVox)
    numPlot = plotVox; end

for p = 1:length(numPlot)
    if mod(p,4)==1
        column = 1;
        niceFig([.1 .2 .8 1],fontSize,1);
    end
    
    %%%%%%%%% plot betas and modelfit
    subplot(3,4,column)
    v = plotVox(p);
    
    for c = 1:length(fits)
        hbar(c) = niceBars([fits(c).condNums],[fits(c).vox(v).betas],[fits(c).vox(v).sems],condColors(c));
        plot(fits(c).condNums,fits(c).vox(v).modelfit,'b-','LineWidth',2);
        %errorbar_fix(hbar(c),condColors(c));
    end
    
    legend([hbar(1:length(fits))],{fits.cond},'location','SouthEast','FontSize',10,'box','off');
    
    xlabel('Stimulus Number'); xlim([0 fits(end).condNums(end)+1]);
    ylabel('Beta Estim.');
    
    titleTx{1} = ['Vox ' num2str(v) ', R^{2}=' num2str(round(fits(1).vox(v).r2)) ' ' fits(1).cond ];
    for t = 2:length(fits)
    titleTx{t} = ['R^{2}=' num2str(round(fits(t).vox(v).r2)) ' ' fits(t).cond];
    end
    title(titleTx);
    
    set(gca,'box','off','FontName','Arial','FontWeight','normal','TickDir','out');
    
    %%%%%%%%% plot betas in space
    subplot(3,4,column+4)
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
    subplot(3,4,column+8)
    
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
    
    column = column+1;
    if column ==5
        titleText = [whichModel '/' whichStim ': ' ROI ', Session = ' session];
        superTitle(titleText,18,.97)
    end
    end


