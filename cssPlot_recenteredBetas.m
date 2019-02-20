% takes the data for each voxel and plots it in space, but aligns each
% voxel's approximate center

addUtils;

clear all; close all;


subjs = {'MG'};
task = '';
expt = 'compPRF';

saveFig = 1;
plotWhat = 'betas';

minR2 = 20;          % cutoff for vox selection
ROIs= standardROIs('EVC');%{'IOG_faces' 'pFus_faces' 'mFus_faces'};%
conds = [1 2 3];    % 1 = inverted, 2 = misaligned, 3 = normal
binSize = 11/5;      % in deg

whichStim = 'photo';
whichModel = 'cssShift';%'kayCSS';%

hems = {'rh' 'lh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data, sort plotVox           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;

for s = 1:length(subjs)
    for r = 1:length(ROIs)
        [session, numRuns] = vpnlSessions(expt,subj{s},[],task); % OPTIONAL: SESSNUM, TASK

        for h = 1:length(hems)
            
            thisROI = [hems{h} '_' ROIs{r}];
            clear fits;
            [~, fitsName] = fitsDirs(dirOf(pwd),expt,session,whichStim,vpnlROI(thisROI,session(1:2)),whichModel);
            if exist(fitsName)>0
            load(fitsName);
            
            if h ==1 && s == 1
                
                roi(r).fits = fits; % do this only initially
            else for c = 1:length(fits)
                    roi(r).fits(c).vox = [roi(r).fits(c).vox fits(c).vox]; % do this on every iter
                end
            end
            else
                fprintf('Missing %s!\n',fitsName);
            end
        end % hems
    end %ROIs
end % subjs



for r = 1:length(ROIs)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % translate fit params into more meaningful terms for invPRF expt              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for c = 1:length(roi(r).fits)
        roi(r).fits(c).vox = readPRFs(roi(r).fits(c).vox,fits(1).ppd,fits(1).res);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % which voxels are we plotting? trim edge values, R2 cutoff                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    plotVox = 1:length(roi(r).fits(1).vox);
    
    for v = 1:length(roi(r).fits(1).vox)
        for c = 1:length(roi(r).fits)
            if  roi(r).fits(c).vox(v).r2 < minR2 ...
                    || trimPRFs(roi(r).fits(c).vox(v).params,roi(r).fits(c).vox(v).betas,fits(1).ppd,fits(1).res);
                plotVox(find(plotVox==v)) = []; end
        end
    end
    
    %load someVox.mat; % if we've previously chosen some voxels to plot
    
    for c = 1:length(roi(r).fits)
        roi(r).fits(c).vox = roi(r).fits(c).vox(plotVox);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  create supertitle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

titleText = [ strTogether(subjs) ' (' hemText(hems) ' voxels R^2 > ' num2str(minR2) '), ' whichStim ' stim, ' whichModel ' model'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if onLaptop figSize = [0 0 .33*length(ROIs) .33*length(conds)];
else figSize = [.2 .2 .15*length(ROIs) .15*length(conds)];
end
f(1) = niceFig(figSize,fontSize);
numPlots = [length(conds) length(ROIs)];  pl = 0;
f(2) = niceFig(figSize,fontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1: in-space plots
% rows are conditions
% colums are ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for cc = 1:length(conds)
    c  = conds(cc);
    for r = 1:length(ROIs)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overall coverage
        pl = pl+1; figure(f(1));
        subplot(numPlots(1),numPlots(2),pl);
        
        edges = [0:binSize*fits(1).ppd:fits(1).res]; % 5 bins
        [bin,counts,indX,indY] = binPRFcenters(roi(r).fits(c).vox,edges);
        
        matSize = 9; matCent = round(matSize/2);
        clear centBetas
        for v = 1:length(roi(r).fits(c).vox)
            centBetas(v,:,:) = nan(matSize,matSize);
            b = reshape(roi(r).fits(c).vox(v).betas,5,5)';
            % row, column
            % cent = [indY(v) indX(v)];
            % bb(cent(1),cent(2)) is the center of the beta values
            startR = matCent-indY(v)+1;
            startC = matCent-indX(v)+1;
            centBetas(v,startR:startR+size(b,1)-1,startC:startC+size(b,2)-1) = b;
            %figure;subplot(1,2,1);imagesc(bb); subplot(1,2,2);imagesc(bigMat)
        end
        
        switch plotWhat
            case 'betas';
                plotInSpace(squeeze(nanmean(centBetas,1)),'Center-Aligned Betas');
            case 'std'
                plotInSpace(squeeze(nanstd(centBetas,1)),'Center-Aligned Betas');
            case 'count'
                plotInSpace(squeeze(sum(~isnan(centBetas))),'Counts');end
            
            set(get(gca,'Title'),'Visible','on');
            title({ROIs{r}},'fontSize',titleSize,'Interpreter','none');
            
            if r == 1 % left edge of plots
                set(get(gca,'YLabel'),'Visible','on');
                ylabel(roi(r).fits(c).cond,'FontSize',titleSize+8,'FontWeight','bold');
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % figure 2: collapsed by eccen
            % rows are conditions
            % colums are ROI
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure(f(2));
            subplot(numPlots(1),numPlots(2),pl);
            res = fits(1).res; ppd = fits(1).ppd;
            
            %binCents = [binSize/2:binSize:(binSize*(matSize-.5))];
            %binCents = binCents-binCents(ceil(matSize/2));
            
            % hardcoded because the rounding was messing things up before in the flex
            % version ughh
            %b = ([1:9]-5)*binSize;
            binCents = [-8.8 -6.6 -4.4 -2.2 0 2.2 4.4 6.6 8.8];
            
            for x = 1:length(binCents)
                for y = 1:length(binCents)
                    % the
                    ecc(x,y) = sqrt(binCents(x)^2+binCents(y)^2);
                end
            end
            
            uniqE = unique(ecc);
            % pos x voxels
            %ee = reshape(ecc,size(ecc,1)*size(ecc,2),1);
            %bb = reshape(centBetas',size(centBetas,2)*size(centBetas,3),size(centBetas,1));
            
            eccSort = struct('betas',[]);
            % collapse centBetas by these eccens
            for b = 1:length(uniqE)
                eccSort(b).ecc = uniqE(b);
                [i,j] = find(ecc==uniqE(b));
                for n = 1:length(i)
                    eccSort(b).betas = [eccSort(b).betas centBetas(:,i(n),j(n))];
                end
                eccSort(b).mean = nanmean(eccSort(b).betas(:));
                eccSort(b).std = nanstd(eccSort(b).betas(:));
                eccSort(b).count = sum(~isnan(eccSort(b).betas(:)));
            end
            
            errorbar([eccSort.ecc],[eccSort.mean],[eccSort.std]./sqrt([eccSort.count]),'Color',condColors(c,1));
            set(gca,'box','off','color','none');
            set(get(gca,'Title'),'Visible','on');
            title([ROIs{r} ' ' roi(r).fits(c).cond],'fontSize',titleSize,'Interpreter','none');
            xlabel('Distance from Binned Center (dva)','fontSize',fontSize);
            ylabel('Mean Beta Value','fontSize',fontSize);
    end
end

for n = 1:2
    figure(f(n));
    superTitle([titleText],titleSize,.97);
    
    if saveFig == 1
        if n == 1
            txt = ['centBetas_' task '_'  plotWhat '_' ];
            niceSave([dirOf(pwd) 'figures/' expt '/betas/'],txt,ROIs,subjs); % just save pngs, since these can be generated pretty quickly
        else
            txt = ['centBetas_' task '_byEccen_' ];
            niceSave([dirOf(pwd) 'figures/' expt '/betas/'],txt,ROIs,subjs); % just save pngs, since these can be generated pretty quickly
        end
    end
end
playSound;
