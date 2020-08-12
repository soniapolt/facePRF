% loads & plots model properties (like size x eccen) along the visual heirarchy for cssFit properties
% 7/27/18: changes from lsline to fitline2derror

clear all; close all;

subjs = {'TH' 'DF' 'EM' 'JG' 'MG' 'SP'};%'SP' 
task = '';
expt = 'fixPRF';

minR2 = 50;          % cutoff for vox selection
ROIs= standardROIs;% 'V3' 'hV4' 'IOG_faces' 'pFus_faces','mFus_faces'
whichCond =2; 
saveFig = 0;

whichStim = 'photo';
whichModel = 'kayCSS';%'cssExpN';%
fitSuffix = '';

hems = {'rh' 'lh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

load(pRFfile(dirOf(pwd),expt,minR2,whichStim,whichModel,hems,fitSuffix,task));
ROInum = cellNum(ROIs,info.ROIs);
subjNum = cellNum(subjs,info.subjs);

if length(subjNum) == 1 roi = subj(subjNum).roi; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create supertitle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titleText = [expt ' (' hemText(hems) ') '];
titleText = [titleText strTogether(subjs) ' (voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
niceFig([.1 .1 .8 .8],fontSize,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1:
% 1) sixe by eccen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:length(ROIs)
    s(r) = scatter([roi(ROInum(r)).fits(whichCond).vox.eccen],[roi(ROInum(r)).fits(whichCond).vox.size],30,condColors(r,1)); hold on;
    %if r~=length(ROIs)
    [h1 h2] = scatterline([roi(ROInum(r)).fits(whichCond).vox.eccen],[roi(ROInum(r)).fits(whichCond).vox.size],[0:.5:7],NaN,100,condColors(r,1),1,1);
    alpha(.2);
    %end
    %[linePar{r} lineR2{r}] = fitl1line([roi(r).fits(whichCond).vox.eccen],[roi(r).fits(whichCond).vox.size]);
    %%% fitline2derror option
    %     [linePar{r} lineR2{r}]=fitline2derror([roi(r).fits(whichCond).vox.eccen],[roi(r).fits(whichCond).vox.size]);
    %     ln = polyval(linePar{r},[roi(r).fits(whichCond).vox.eccen]);
    %     l = plot([roi(r).fits(whichCond).vox.eccen],ln,'Color', condColors(r,1));
end

% l=lsline; l=fliplr(l); for r = 1:length(ROIs) set(l(r),'Color', condColors(r,1),'LineWidth',2); end

xlabel('Eccen (dva)','FontSize',fontSize); ylabel(['Size (2*SD/sqrt(N)) (dva)'],'FontSize',fontSize); title('Size by Eccentricity','fontSize',titleSize);%ylim([0 2.5]);
axis square; g = legend([s(fliplr(1:r))],fliplr(ROIs)); set(g,'box','off','FontSize',fontSize+4,'location','NorthEastOutside','Interpreter','none');
xlim([0 12]); ylim([0 12]);


superTitle(titleText,titleSize);

if saveFig == 1
    txt = ['sizeEccen_' whichModel ];
    if ~containsTxt(whichStim,'photo')
            txt = [whichStim '_' txt];  end
    niceSave([dirOf(pwd) 'figures/' expt '/sizeEccen/'],txt,[],subjs); % just save pngs, since these can be generated pretty quickly
end
 playSound;