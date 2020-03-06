% simulates overlap of pRF coverage and mean face features at different
% locations
% annoyingly, this does require PTB...
% _parallel gives us fewer checkpoints/output but should be faster

clear all; close all;


expt = 'fixPRF'; % assumes all subject in pRFset
minR2 = ['r2-20'];          % cutoff for vox selection
ROIs = {'IOG_faces' 'pFus_faces' 'mFus_faces' 'pSTS_faces'}; % if more than one ROI, combine voxel positions for them

whichStim = 'outline';%'photo';%'internal';%
whichModel = 'kayCSS';%'cssShift';%
hems = {'rh' 'lh'};
recomp = 1;


sim.numSims = 1000;%50;
sim.drawVox = .8; % now a proportion, vs. absolute number
sim.runIndivs = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in data, sort plotVox           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontSize = 11; titleSize = 14;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load coverage maps                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(ROIs)
  sim.name = [raid 'invPRF/fixPRF/behavSim/' ROIs{r} '_' minR2];  
for c = 1:2
    if ~exist([sim.name '_cond' num2str(c) '.mat']) || recomp
        fprintf(['Starting ' [sim.name '_cond' num2str(c)] '.mat: %s...\n'],datestr(now)); tic;
    sim.roi = []; 
% in this version, we load in pre-made coverage ims of all pRFs rather
% than calculating them using the PRF() functions

%%% aggregate voxels across our ROIs

    try
        imFile = [raid 'invPRF/fixPRF/prfIms/' ROIs{r} '_' fileName(pRFfile('',expt,minR2,whichStim,whichModel,hems)) '_cond' num2str(c) '.mat'];
        load(imFile)
    catch
        error(sprintf('Missing %s! Check filepath or run analysis_makePRFIms.m\n',imFile)); end
    vox = cPRFs;
end

sim.ppd = vis.ppd;% now taking directly from the pre-made coverages
sim.res = vis.res; % give a little buffer around our image/make indexing easier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation parameters                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XX,YY]=meshgrid(-sim.res/2:sim.res/2);
sim.steps = -3:.1:3; [X,Y]=(meshgrid(sim.steps));
sim.degPos = [Y(:) X(:)];
sim.centers = sim.degPos*sim.ppd+sim.res/2+1;
sim.faceSize = 3.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% internal features images             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('internalAvg.mat');
faceSize = sim.faceSize*sim.ppd;
face = imresize(avgFace,[faceSize faceSize]);
if c == 1 face = flipud(face); end

%%% 12/18/19 hack to get around nefesh's lack of PTB and slow VPN - remove!
try
    % make simCenters array - just do this once per simulation
    for n = 1:length(sim.centers)
        % make current face im
        fIm = zeros(sim.res+1,sim.res+1); co = CenterRectOnPoint([1 1 faceSize faceSize],sim.centers(n,1),sim.centers(n,2));
        fIm(co(1):co(3),co(2):co(4)) = face;
        faceIm{n} = fIm;
    end
catch
    load('ptbHack_faceIm.mat');
    eval(['faceIm = faceIm' num2str(c) ';']); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run group simulation - CURRENTLY NOT USED                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear draw; draw(sim.numSims) = struct;
% tic
% for s = 1:sim.numSims
%     if size(vox,1)<sim.drawVox drawVox = size(vox,1); else drawVox = sim.drawVox; end
%     covIm = vox(datasample(1:size(vox,1),drawVox,'Replace',true),:,:);
%     
%     % to deal with the issue of plot/imagesc using different coordinate systems
%     covIm = flipud(squeeze(mean(covIm)));
%     
%     %subplot(1,3,1); imshow(covIm); colormap(mrvColorMaps('hot'));title([ROIs{r} ' coverage']);
%     
%     %subplot(1,3,2); imshow(faceIm{1});
%     for n = 1:length(sim.centers)
%         featCov = faceIm{n}.*covIm; %subplot(1,3,3); imshow(featCov);
%         draw(s).result(n) = sum(featCov(:));
%     end
%     
%     draw(s).best = find(draw(s).result==max(draw(s).result));
%     draw(s).bestDegX = [sim.degPos(draw(s).best,2)];
%     draw(s).bestDegY = [sim.degPos(draw(s).best,1)];
%     % toc
%     if mod(s,100)==0
%         %save([sim.name '_cond' num2str(c)],'draw','sim','face');
%         fprintf('Done with %i draws...\n',s);
%         toc
%     end
%  end
% 
%     %%%% aggregate group results
%     sim.roi.name = sim.name;
%     sim.roi.bestPos = [mean([draw.bestDegX]) mean([draw.bestDegY])];
%     sim.roi.sd = [std([draw.bestDegX]) std([draw.bestDegY])];
%     fprintf(['Saving group sim: ' [sim.name '_cond' num2str(c)] '\n']);toc;
%     save([sim.name '_cond' num2str(c)],'draw','sim','face');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run indiv simulation                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sim.runIndivs
    
    allVox = vox;
    tic
    for n = 1:length(unique(vis.subjInd))
        clear draw; draw(sim.numSims) = struct;
        vox = allVox(vis.subjInd==n,:,:);
        % tic
        drawVox = round(size(vox,1)*sim.drawVox);
         if size(vox,1)<sim.drawVox drawVox = size(vox,1); end
        for s = 1:sim.numSims
            covIm = vox(datasample(1:size(vox,1),drawVox,'Replace',true),:,:);

            % to deal with the issue of plot/imagesc using different coordinate systems
            covIm = flipud(squeeze(mean(covIm)));

            for m = 1:length(sim.centers)
                featCov = faceIm{m}.*covIm; 
                draw(s).result(m) = sum(featCov(:));
            end
            
            
            draw(s).best = find(draw(s).result==max(draw(s).result));
            draw(s).bestDegX = [sim.degPos(draw(s).best,2)];
            draw(s).bestDegY = [sim.degPos(draw(s).best,1)];
            
            if mod(s,100)==0
                fprintf('Subj %d, done with %i draws...\n',n,s);
                toc
            end
        end    
            %%%% aggregate subj results
            sim.roi(r).subj(n).bestPos = [mean([draw.bestDegX]) mean([draw.bestDegY])];
            sim.roi(r).subj(n).sd = [std([draw.bestDegX]) std([draw.bestDegY])];
            save([sim.name '_cond' num2str(c)],'draw','sim','face');
            fprintf(['Saving Subj ' num2str(n) ' sim: ' [sim.name '_cond' num2str(c)] '\n']); toc;

    end
end
    else fprintf(['Wont recompute ' [sim.name '_cond' num2str(c)] '...\n']); end
end
% end

%if onLaptop playSound; end

%[bestPos,sim] = prfRec_drawBestFace(draw,sim,face); title([strTogether(ROIs) ' - Cond ' num2str(c)]);
% niceSave([raid 'invPRF/figures/fixPRF/simulation/' num2str(sim.numSims) '_' strTogether(ROIs,0,'-') '_cond' num2str(c)]);
%endf