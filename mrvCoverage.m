function [RFcov, figHandle, all_models, weight, data] = mrvCoverage(vox,titleText)
% translates rmPlotCoverage() to work with prfSet framework
              
%rmPlotCoverage - calulate the visual field coverage within an ROI
% 
% [RFcov figHandle all_models weight data]  = rmPlotCoverage(vw, varargin)
%
%
% Before you run this script, you have to load 'variance explained', 'eccentricity',
% 'polar-angle' and 'prf size' into 'co', 'map', 'ph' and 'amp' fields, respectively
% 
% OUTPUT
%  RFcov
% INPUT
%  prf_size:        0 = plot pRF center; 1 = use pRF size
%  fieldRange:      maximum eccentricity to plot (deg)
%  method:          'sum','max', 'clipped average', 'signed profile'
%  newfig:          make a new figure (1) or not (0). (-1 indicates don't plot
%                       anything, just return the coverage map.)
%  nboot:           the number of bootstrapping (0 means no bootstrapping)
%  normalizeRange:  if true, scale z axis to [0 1]
%  smoothSigma:     median smoothing default: 2 nearest values
%  cothresh:        threshold by variance explained in model
%  eccthresh:       2-vector ecc limits (default = [0 1.5*fieldRange])
%  nsamples:        num samples in square grid (default = 128)
%  weight:          any of {'fixed', 'parameter map', 'variance explained'} (default = fixed) 
%  weightBeta:      use beta values from rmModel to weight pRFs (default = false)
%  addcenters:      1 = superimpose dots for pRF centers; 0 = do not show centers
%

% using DF's default parameters
vfc.title = titleText;
vfc.prf_size = true; 
vfc.fieldRange = 4.5;
vfc.method = 'maximum profile';         
vfc.newfig = 1;                      
vfc.nboot = 100;                          
vfc.normalizeRange = false;              
vfc.smoothSigma = false;                
vfc.cothresh = 0;         % we've already done this at prfSet stage
vfc.eccthresh = [0 1.5*vfc.fieldRange]; % we've already done this at prfSet stage
vfc.nSamples = 500;            
vfc.meanThresh = 0;
vfc.weight = 'variance explained';  
vfc.weightBeta = 0;
vfc.cmap = 'jet';						
vfc.clipn = 'fixed';                    
vfc.addCenters = true;                 
vfc.verbose = prefsVerboseCheck;

compVolume = false;


% Get co and ph (vectors) for the current scan, within the
% current ROI.
co      = [vox.r2];
sigma1  = [vox.sdDeg];
sigma2  = [vox.sdDeg];
theta   = zeros(1,length(vox)); % angle of sigmaMajor (radians, 0=vertical)
beta    = [vox.gain];
xy = vertcat(vox.XYdeg);
x0      = xy(:,1);
y0      = xy(:,2);


% grabbing both (x0, y0) and (pol, ecc) are redundant, and allow for the
% two specifications to get separated (because of the y-flip issue). So,
% re-express ecc and pol using x0 and y0.
[ph , ecc] = cart2pol(x0, y0);

% smooth sigma
if vfc.smoothSigma
    
    % Cannot do sigma smoothing if the ROI has less than 3 voxels. 
    if size(x0,1) < 3
        error('Cannot perform sigma smoothing when the ROI has less than 3 voxels. ')
    end
    
    
    if vfc.smoothSigma == 1
        vfc.smoothSigma = 3; %default
    end
    n = vfc.smoothSigma;
	
    %check sigma1==sigma2
    if sigma1 == sigma2
        for ii = 1:length(sigma1)
            %compute nearest coords
            dev = sqrt(abs(x0(ii) - x0).^2 + abs(y0(ii) - y0).^2);
            [dev, ix] = sort(dev); %#ok<*ASGLU>
            sigma1(ii) = median(sigma1(ix(1:n)));           
        end
        sigma2 = sigma1;
    else
        for ii = 1:length(sigma1)
            %compute nearest coords
            dev = sqrt(abs(x0(ii) - x0).^2 + abs(y0(ii) - y0).^2);
            [dev, ix] = sort(dev);
            sigma1(ii) = median(sigma1(ix(1:n)));
            sigma2(ii) = median(sigma2(ix(1:n)));
        end
    end
end
             
% polar plot
subX = single(ecc .* cos(ph));
subY = single(ecc .* sin(ph));


% visual field
x = single( linspace(-vfc.fieldRange, vfc.fieldRange, vfc.nSamples) );
[X,Y] = meshgrid(x,x);

% gather this data to make accessible in the plot
if vfc.newfig  > -1   % -1 is a flag that we shouldn't plot the results
	data.figHandle = gcf;
	data.co        = co;
	data.ph        = ph;
	data.subCo     = co;
	data.subPh     = ph;
	data.subEcc    = ecc;
	data.subx0     = x0;
	data.suby0     = y0;
    data.subSize1  = sigma1;
    data.subSize2  = sigma2;
    data.X         = X;
	data.Y         = Y;
end

% For the pRF center plot, use a small constant pRF size
if vfc.prf_size==0
   sigma1 = ones(size(sigma1)) * 0.1;
   sigma2 = ones(size(sigma2)) * 0.1;   
   theta = zeros(size(theta));   
end

switch lower(vfc.weight)
    case 'fixed'
        weight = ones(size(co));
        
    case {'variance explained', 'varexp', 've'}
        weight = co;
        
    otherwise 
        error('Unknown weight parameter: %s',vfc.weight);
end

if vfc.weightBeta==1
    weight = weight .* beta(coIndices);
end
weight = single(weight);

%%% special case: for the 'density' coverage option, we don't need to
%%% do a lot of memory-hungry steps like making all pRFs. So, I've set those
%%% computations aside in their own subroutine. (ras)
if isequal( lower(vfc.method), 'density' )
	RFcov = prfCoverageDensityMap(x0, y0, sigma1, X, Y);
	
	all_models = []; % not created for this option	
	if vfc.newfig==-1
		figHandle = [];
	else
		figHandle = createCoveragePlot(RFcov, vfc, roi, data);
	end

	return
end


%% make all pRFs:
% make in small steps so we don't go into swap space for large ROIs
n = numel(subX);
s = [(1:ceil(n./1000):n-2) n+1]; 

% For the line above (which assumes that we have at least 3 voxels,
% probably for median smoothing),  s is an empty vector when n < 3 
% -- and an empty rfcov is  returned when we try to plot the 
% coverage. So modify s accordingly for these edge cases: 
% The definition of s is not very intuitive
if n < 3
    s = [1:n+1]; 
end

all_models = zeros( numel(X), n, 'single' );
fprintf(1,'[%s]:Making %d pRFs:...', mfilename, n);
drawnow;
for n=1:numel(s)-1
    % make rfs
    rf   = rfGaussian2d(X(:), Y(:),...
						sigma1(s(n):s(n+1)-1), ...
						sigma2(s(n):s(n+1)-1), ...
						theta(s(n):s(n+1)-1), ...
						subX(s(n):s(n+1)-1), ...
						subY(s(n):s(n+1)-1));
    all_models(:,s(n):s(n+1)-1) = rf;
end
clear n s rf pred;
fprintf(1, 'Done.\n');
drawnow;

% Correct volume
if compVolume
    tmp = ones(size(all_models, 1), 1, 'single');
    
    vol = sigma1(coIndices).^2;
    vol = vol * (2 * pi);
    
    all_models = all_models ./ (tmp * vol);
end

% For the pRF center plot, put a constant value (1) within each Gaussian
if vfc.prf_size==0
    all_models(all_models>0.1)=1;
end

% weight all models
if isequal( lower(vfc.weight), 'fixed' )
	% if the weights are even, we avoid the redundant, memory-hungry
	% multiplication step that would otherwise be done. 
	all_models_weighted = all_models;
else
	tmp = ones(size(all_models, 1), 1, 'single');
	all_models_weighted = all_models .* (tmp * weight);
	clear tmp
end

%% Different ways of combining them: 
% 1) bootstrap (yes, no) 2) which statistic (sum, max, etc), 
% bootstrap

% If we are only working with 1 voxel, bootstrapping does not do anything. 
% We turn it off because the bootstp function does not handle this case well. 

if vfc.nboot>0
    if isempty(which('bootstrp'))
        warndlg('Bootstrap requires statistics toolbox');
        RFcov = [];
        return;
    end
    all_models(isnan(all_models))=0;

    switch lower(vfc.method)
        case {'sum','add','avg','average everything', 'average'}
            m = bootstrp(vfc.nboot, @mean, all_models');
        
        case {'max','profile','maximum profile' 'maximum'}
            m = bootstrp(vfc.nboot, @max, all_models');
        
        otherwise
            error('Unknown method %s',vfc.method)
    end
    RFcov=median(m,1)';
    
% no bootstrap
else
    switch lower(vfc.method)
                    
        % coverage = sum(pRF(i)*w(i)) / (sum(pRF(i))
        case {'beta-sum','betasum','weight average'}
            RFcov = sum(all_models_weighted, 2) ./ sum(all_models,2);
            
        % coverage = sum(pRF(i)*w(i)) / (sum(pRF(i)) + clipping
        case {'clipped beta-sum','clippedbeta','clipped weight average'}
            % set all pRF beyond 2 sigmas to zero
            clipval = exp( -.5 *((2./1).^2));
            all_models(all_models<clipval) = 0;
            n = all_models > 0;
            
            % recompute all_models_weighted
			tmp = ones( size(all_models,1), 1, 'single' );
            all_models_weighted = all_models .* (tmp*weight);
            
            % compute weighted clipped sum/average
            sumn = sum(n,2);
            mask = sumn==0;
            sumn(mask) = 1; % prevent dividing by 0
            RFcov = sum(all_models_weighted,2) ./ sum(all_models,2);
            RFcov(mask) = 0;
            
            %clip to zero if n<clipn
            if isnumeric(vfc.clipn)
                RFcov(sumn<=vfc.clipn) = 0;
            end            
           
        % coverage = sum(pRF(i)*w(i)) / (sum(w(i))
        case {'sum','add','avg','average','prf average'}
            RFcov = sum(all_models_weighted, 2) ./ sum(weight);
        
        % coverage = sum(pRF(i)*w(i)) / (sum(w(i)) + clipping
        case {'clipped average','clipped','clipped prf average'}
            % set all pRF beyond 2 sigmas to zero
            clipval = exp( -.5 *((2./1).^2));
            all_models(all_models<clipval) = 0;
            n = all_models > 0;
            
            % recompute all_models_weighted
			tmp = ones( size(all_models,1), 1, 'single' );
            all_models_weighted = all_models .* (tmp*weight);
            
            % compute weighted clipped mean
            sumn = sum(weight.*n);
            mask = sumn==0;
            sumn(mask) = 1; % prevent dividing by 0
            RFcov = sum(all_models_weighted,2) ./ sumn;
            RFcov(mask) = 0;
            
            %clip to zero if n<clipn
            if isnumeric(vfc.clipn)
                RFcov(sumn<=vfc.clipn) = 0;
            end
            
        % coverage = max(pRF(i))
        case {'maximum profile', 'max', 'maximum'}
            RFcov = max(all_models_weighted,[],2);
            
        case {'signed profile'}
            RFcov  = max(all_models_weighted,[],2);
            covmin = min(all_models_weighted,[],2);
            ii = RFcov<abs(covmin);
            RFcov(ii)=covmin(ii);
            
        case {'p','probability','weighted statistic corrected for upsampling'}
            RFcov = zeros(vfc.nSamples);
			
			% I guess this upsample factor assumes your functional data are
			% 2.5 x 2.5 x 3 mm?
            upsamplefactor = 2.5*2.5*3; % sigh.....
            for ii = 1:size(all_models,1)
                s = wstat(all_models(ii,:),weight,upsamplefactor);
                if isfinite(s.tval)
                    RFcov(ii) = 1 - t2p(s.tval,1,s.df);
                end
            end

        otherwise
            error('Unknown method %s',vfc.method)
    end
end

% convert 1D to 2D
RFcov = reshape( RFcov, [1 1] .* sqrt(numel(RFcov)) );

% When no voxels exceed threshold, return nan matrix rather than empty
% matrix
if sum(size(RFcov))==0
    RFcov=nan(nSamples,nSamples);
end

% if the newfig flag is set to -1, just return the image
if vfc.newfig==-1
    figHandle = [];
else
	figHandle = createCoveragePlot(RFcov, vfc, data);
end


return
% /--------------------------------------------------------------------/ %
end


% /--------------------------------------------------------------------/ %
function figHandle = createCoveragePlot(RFcov, vfc,  data)
% plotting subroutine for rmPlotCoverage. Broken off by ras 10/2009.
% if vfc.newfig
%     figHandle = figure('Color', 'w');
% else
	figHandle = selectGraphWin;
%end

%set(gcf, 'Name', vfc.title);


% normalize the color plots to 1
if vfc.normalizeRange
	rfMax = max(RFcov(:)); 
else
	rfMax = 1; 
end

img = RFcov ./ rfMax;
mask = makecircle(length(img));
img = img .* mask;
imagesc(data.X(1,:), data.Y(:,1), img);
set(gca, 'YDir', 'normal');
grid on

colormap(vfc.cmap);
colorbar;

% start plotting
hold on;

% add polar grid on top
p.ringTicks = 1:ceil(vfc.fieldRange);%(1:3)/3*vfc.fieldRange;
p.color = 'w';
polarPlot([], p);

% add pRF centers if requested
if vfc.addCenters
    inds = data.subEcc < vfc.fieldRange;
    plot(data.subx0(inds), data.suby0(inds), '.', ...
		'Color', [.5 .5 .5], 'MarkerSize', 4); 
end


% scale z-axis
if vfc.normalizeRange
	if isequal( lower(vfc.method), 'maximum profile' )
		caxis([.5 1]);
	else
	    caxis([0 1]);
	end
else
    if min(RFcov(:))>=0
        caxis([0 ceil(max(RFcov(:)))]);
    else
        caxis([-1 1] * ceil(max(abs(RFcov(:)))));
    end
end
axis image;   % axis square;
xlim([-vfc.fieldRange vfc.fieldRange])
ylim([-vfc.fieldRange vfc.fieldRange])

title(vfc.title, 'FontSize', 24, 'Interpreter', 'none');

% Save the data in gca('UserData')
set(gca, 'UserData', data);

return;
% /------------------------------------------------------------------/ %

end


% /------------------------------------------------------------------/ %
function RFcov = prfCoverageDensityMap(x0, y0, sigma, X, Y) 
% for each point (x, y) in visual space, this returns
% the proportion of voxels in the ROI for which (x, y) is
% within one standard deviation of the pRF center.
mask = NaN( size(X, 1), size(X, 2), length(x0) );

for v = 1:length(x0)
	% make a binary mask within one sigma of the center
	R = sqrt( (X - x0(v)) .^ 2 + (Y - y0(v)) .^ 2 );
	mask(:,:,v) = ( R < 2*sigma(v) );
end

% average (sum?) across all masks
RFcov = nansum(mask, 3);

return
% /------------------------------------------------------------------/ %



end