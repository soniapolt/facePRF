% loads & plots some parameter after binning pRF by their centers

clear all; close all;

whichParam = 'size'; % cases: eccen, size, baseline

expt = 'compPRF';
subjs = {'SP'};

saveFig = 1;

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs;
conds = [1 2 3];  %
binSize = 1;      % in deg

whichStim = 'photo';
whichModel = 'kayCSS';%'cssShift';%
fitSuffix = '_new';

hems = {'lh' 'rh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data, sort plotVox           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;
% now we load in the data from both hemispheres, and threshold across
load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create supertitle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titleText = [expt ' ' whichParam ' '  strTogether(subjs) ' (' hemText(hems) ' voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if onLaptop figSize = [0 0 1 .33*length(conds)];
else figSize = [.2 .2 .5 .15*length(conds)];
end
niceFig(figSize,fontSize);
numPlots = [length(conds) length(ROIs)];  pl = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1:
% rows are conditions
% colums are ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for cc = 1:length(conds)
    c  = conds(cc);
    for r = 1:length(ROIs)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overall coverage
        subplot(numPlots(1),numPlots(2),pl); pl = pl+1;
        if length(subjNum) == 1 prfs = subj(subjNum).roi(ROInum(r)); else prfs = roi(ROInum(r)); end

        edges = [0:binSize*prfs.fits(1).ppd:prfs.fits(1).res]; % 2-deg bins
        [bin,counts] = binPRFcenters(prfs.fits(c).vox,edges);
        
        for b = 1:length(bin)
            switch whichParam
                case 'baseline'
                    plotVals{b} = [prfs.fits(c).vox(bin(b).vox).baseline];
                    clim = [-1.5 1.5];
                case 'eccen'
                    plotVals{b} = [prfs.fits(c).vox(bin(b).vox).eccen];
                    clim = [];
                case 'size'
                    plotVals{b} = [prfs.fits(c).vox(bin(b).vox).size];
                    clim = [];
            end
            plotBin(b) = nanmean(plotVals{b});
        end
        valInSpace= reshape(plotBin,length(edges)-1,length(edges)-1);
        plotInSpace(valInSpace,[],[],0,clim);
        
        set(get(gca,'Title'),'Visible','on');
        title({ROIs{r}},'fontSize',titleSize,'Interpreter','none','Color',condColors(r+3,1));
        
        if r == 1 % left edge of plots
            set(get(gca,'YLabel'),'Visible','on');
            ylabel(roi(r).fits(c).cond,'FontSize',titleSize+8,'FontWeight','bold');
        end
    end
end


superTitle([titleText],titleSize,.97);

if saveFig == 1
    txt = [whichParam 'Binned'];
    niceSave([dirOf(pwd) 'figures/' expt '/binned/'],txt,[],subjs); % just save pngs, since these can be generated pretty quickly
end

playSound;
