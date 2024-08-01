function [beim,cp]=fingerenhance2(fpath) 

    %close all;clc;    
    %dname = 'Dataset\database1\49\fingerprint';
    %fname = '107_1.tif';
    %fpath = sprintf('%s\\%s',dname,fname);
    
    [im,map]=imread(fpath);
    %figure(1);imshow(im,map); title('Input fingerprint Image');
    im = im2double(im);
    [M,N,C]=size(im);
    if(C>1)
    im= rgb2gray(im);
    end

    % Identify ridge-like regions and normalise image
    blksze = 16; thresh = 0.2;
    [normim, mask] = ridgesegment(im, blksze, thresh);
    se = strel('disk',32);  
    mask=imopen(mask,se);
 
    %figure;subplot(1,2,1);imshow(normim);title('Input Fingerprint');
    %subplot(1,2,2);imshow(mask); title('Fingerprint mask (ROI)');
    
    % Determine ridge orientations
    [orientim, reliability] = ridgeorient(normim, 1, 5, 5);
    %plotridgeorient(orientim, 20, im, 2)
    %show(reliability,6)
    
    % Determine ridge frequency values across the image
    blksze = 36; 
    [freq, medfreq] = ridgefreq(normim, mask, orientim, blksze, 5, 5, 15);
    %figure;imshow(freq);title('Ridge Frequency');
    %show(freq,3) 
    
    % Actually I find the median frequency value used across the whole
    % fingerprint gives a more satisfactory result...
    freq = medfreq.*mask;   
    % Now apply filters to enhance the ridge pattern
    newim = ridgefilter(normim, orientim, freq, 0.5, 0.5);
    %show(newim,4);
    %figure;imshow(newim);title('Ridge filtered Image');
        
   
    CC = bwconncomp(mask);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    mask = zeros(size(newim));
    mask(CC.PixelIdxList{idx}) = 1;  
    
    beim = (newim > 0).*mask;    
    stats=regionprops(beim,'BoundingBox'); 
    bbox=stats.BoundingBox;  
    beim=imcrop(beim,bbox); 

    orientim=imcrop(orientim,bbox); 
    mask=imcrop(mask,bbox); 
    
    a=cos(orientim);
    b=sin(orientim);
    [afx,afy] = gradient(a);
    [bfx,bfy] = gradient(b);
    
    ma = sqrt( afx.^2 + afy.^2);
    mb = sqrt( bfx.^2 + bfy.^2);    
    mab = min(ma,mb);
    
    Vx= max(mab(:));    
    Cp=mab>0.3*Vx;
    Cp=Cp.*mask;
    
    [I,J]=find(Cp);
    ds=sqrt( (I-M/2).^2 + (J-N/2).^2);
    [~,mi]=min(ds);
    cp= [I(mi),J(mi)];
    
    

    

function [orientim, reliability, coherence] = ridgeorient(im, gradientsigma, blocksigma, orientsmoothsigma)
        
    if ~exist('orientsmoothsigma', 'var'), orientsmoothsigma = 0; end
    [rows,cols] = size(im);
    
    % Calculate image gradients.
    [Gx, Gy] = derivative5(gaussfilt(im, gradientsigma), 'x', 'y');
    
    % Estimate the local ridge orientation at each point by finding the
    % principal axis of variation in the image gradients.
    Gxx = Gx.^2;       % Covariance data for the image gradients
    Gxy = Gx.*Gy;
    Gyy = Gy.^2;
    
    % Now smooth the covariance data to perform a weighted summation of the
    % data.
    sze = fix(6*blocksigma);   if ~mod(sze,2); sze = sze+1; end    
    f = fspecial('gaussian', sze, blocksigma);
    Gxx = filter2(f, Gxx); 
    Gxy = 2*filter2(f, Gxy);
    Gyy = filter2(f, Gyy);
    
    % Analytic solution of principal direction
    denom = sqrt(Gxy.^2 + (Gxx - Gyy).^2) + eps;
    sin2theta = Gxy./denom;            % Sine and cosine of doubled angles
    cos2theta = (Gxx-Gyy)./denom;

    if orientsmoothsigma
        sze = fix(6*orientsmoothsigma);   if ~mod(sze,2); sze = sze+1; end    
        f = fspecial('gaussian', sze, orientsmoothsigma);    
        cos2theta = filter2(f, cos2theta); % Smoothed sine and cosine of
        sin2theta = filter2(f, sin2theta); % doubled angles
    end
    
    orientim = pi/2 + atan2(sin2theta,cos2theta)/2;      
    
    Imin = (Gyy+Gxx)/2 - (Gxx-Gyy).*cos2theta/2 - Gxy.*sin2theta/2;
    Imax = Gyy+Gxx - Imin;
    
    reliability = 1 - Imin./(Imax+.001);
    coherence = ((Imax-Imin)./(Imax+Imin)).^2;    
    reliability = reliability.*(denom>.001);
    
    

function [normim, mask, maskind] = ridgesegment(im, blksze, thresh)
    
    im = normalise(im,0,1);  % normalise to have zero mean, unit std dev
    
    fun = inline('std(x(:))*ones(size(x))');
    
    stddevim = blkproc(im, [blksze blksze], fun);
    
    mask = stddevim > thresh;
    maskind = find(mask);
    
    % Renormalise image so that the *ridge regions* have zero mean, unit
    % standard deviation.
    im = im - mean(im(maskind));
    normim = im/std(im(maskind));    

    

function [freq, medianfreq] = ridgefreq(im, mask, orient, blksze, windsze,minWaveLength, maxWaveLength)                                                
    [rows, cols] = size(im);
    freq = zeros(size(im));
    
    for r = 1:blksze:rows-blksze
        for c = 1:blksze:cols-blksze
          blkim = im(r:r+blksze-1, c:c+blksze-1);    
          blkor = orient(r:r+blksze-1, c:c+blksze-1);       
          
          freq(r:r+blksze-1,c:c+blksze-1) =  ...
              freqest(blkim, blkor, windsze, minWaveLength, maxWaveLength);
        end
    end

    % Mask out frequencies calculated for non ridge regions
    freq = freq.*mask;
    
    % Find median freqency over all the valid regions of the image.
    medianfreq = median(freq(find(freq>0)));  
    
    
    
function freqim =  freqest(im, orientim, windsze, minWaveLength, maxWaveLength)
    
    
    [rows,cols] = size(im);
    
    % Find mean orientation within the block. This is done by averaging the
    % sines and cosines of the doubled angles before reconstructing the
    % angle again.  This avoids wraparound problems at the origin.
    orientim = 2*orientim(:);    
    cosorient = mean(cos(orientim));
    sinorient = mean(sin(orientim));    
    orient = atan2(sinorient,cosorient)/2;

    % Rotate the image block so that the ridges are vertical
    rotim = imrotate(im,orient/pi*180+90,'nearest', 'crop');
    
    % Now crop the image so that the rotated image does not contain any
    % invalid regions.  This prevents the projection down the columns
    % from being mucked up.
    cropsze = fix(rows/sqrt(2)); offset = fix((rows-cropsze)/2);
    rotim = rotim(offset:offset+cropsze, offset:offset+cropsze);

    % Sum down the columns to get a projection of the grey values down
    % the ridges.
    proj = sum(rotim);
    
    % Find peaks in projected grey values by performing a greyscale
    % dilation and then finding where the dilation equals the original
    % values. 
    dilation = ordfilt2(proj, windsze, ones(1,windsze));
    maxpts = (dilation == proj) & (proj > mean(proj));
    maxind = find(maxpts);

    % Determine the spatial frequency of the ridges by divinding the
    % distance between the 1st and last peaks by the (No of peaks-1). If no
    % peaks are detected, or the wavelength is outside the allowed bounds,
    % the frequency image is set to 0
    if length(maxind) < 2
	freqim = zeros(size(im));
    else
	NoOfPeaks = length(maxind);
	waveLength = (maxind(end)-maxind(1))/(NoOfPeaks-1);
	if waveLength > minWaveLength & waveLength < maxWaveLength
	    freqim = 1/waveLength * ones(size(im));
	else
	    freqim = zeros(size(im));
	end
    end

    
function newim = ridgefilter(im, orient, freq, kx, ky)

    %if nargin == 5
    %    showfilter = 0;
    %end
    
    angleInc = 3;  % Fixed angle increment between filter orientations in
                   % degrees. This should divide evenly into 180
    
    im = double(im);
    [rows, cols] = size(im);
    newim = zeros(rows,cols);
    
    [validr,validc] = find(freq > 0);  % find where there is valid frequency data.
    ind = sub2ind([rows,cols], validr, validc);

    % Round the array of frequencies to the nearest 0.01 to reduce the
    % number of distinct frequencies we have to deal with.
    freq(ind) = round(freq(ind)*100)/100;
    
    % Generate an array of the distinct frequencies present in the array
    % freq 
    unfreq = unique(freq(ind)); 
    
    % Generate a table, given the frequency value multiplied by 100 to obtain
    % an integer index, returns the index within the unfreq array that it
    % corresponds to
    freqindex = ones(100,1);
    for k = 1:length(unfreq)
        freqindex(round(unfreq(k)*100)) = k;
    end
    
    % Generate filters corresponding to these distinct frequencies and
    % orientations in 'angleInc' increments.
    filter = cell(length(unfreq),180/angleInc);
    sze = zeros(length(unfreq),1);
    
    for k = 1:length(unfreq)
        sigmax = 1/unfreq(k)*kx;
        sigmay = 1/unfreq(k)*ky;
        
        sze(k) = round(3*max(sigmax,sigmay));
        [x,y] = meshgrid(-sze(k):sze(k));
        reffilter = exp(-(x.^2/sigmax^2 + y.^2/sigmay^2)/2)...
                .*cos(2*pi*unfreq(k)*x);

        % Generate rotated versions of the filter.  Note orientation
        % image provides orientation *along* the ridges, hence +90
        % degrees, and imrotate requires angles +ve anticlockwise, hence
        % the minus sign.
        for o = 1:180/angleInc
            filter{k,o} = imrotate(reffilter,-(o*angleInc+90),'bilinear','crop'); 
        end
    end

    %if showfilter % Display largest scale filter for inspection
    %    figure(7), imshow(filter{1,end},[]); title('filter'); 
    %end
    
    % Find indices of matrix points greater than maxsze from the image
    % boundary
    maxsze = sze(1);    
    finalind = find(validr>maxsze & validr<rows-maxsze & ...
                    validc>maxsze & validc<cols-maxsze);
    
    % Convert orientation matrix values from radians to an index value
    % that corresponds to round(degrees/angleInc)
    maxorientindex = round(180/angleInc);
    orientindex = round(orient/pi*180/angleInc);
    i = find(orientindex < 1);   orientindex(i) = orientindex(i)+maxorientindex;
    i = find(orientindex > maxorientindex); 
    orientindex(i) = orientindex(i)-maxorientindex; 

    % Finally do the filtering
    for k = 1:length(finalind)
        r = validr(finalind(k));
        c = validc(finalind(k));

        % find filter corresponding to freq(r,c)
        filterindex = freqindex(round(freq(r,c)*100));
        
        s = sze(filterindex);   
        newim(r,c) = sum(sum(im(r-s:r+s, c-s:c+s).*filter{filterindex,orientindex(r,c)}));
    end
    
    
function plotridgeorient(orient, spacing, im, figno)
    
    [rows, cols] = size(orient);
    
    lw = 2;             % linewidth
    len = 0.8*spacing;  % length of orientation lines

    % Subsample the orientation data according to the specified spacing
    s_orient = orient(spacing:spacing:rows-spacing, ...
		      spacing:spacing:cols-spacing);

    xoff = len/2*cos(s_orient);
    yoff = len/2*sin(s_orient);    
    
    
    % Determine placement of orientation vectors
    [x,y] = meshgrid(spacing:spacing:cols-spacing, ...
		     spacing:spacing:rows-spacing);
    
    x = x-xoff;
    y = y-yoff;
    
    % Orientation vectors
    u = xoff*2;
    v = yoff*2;
    
    
    figure;imshow(im); hold on;
    quiver(x,y,u,v,0,'.','linewidth',1, 'color','r');    
    axis equal, axis ij,  hold off
    title('Ridge Orientation');
        
    
function varargout = derivative5(im, varargin)

    varargin = varargin(:);
    varargout = cell(size(varargin));
    
    % Check if we are just computing 1st derivatives.  If so use the
    % interpolant and derivative filters optimized for 1st derivatives, else
    % use 2nd derivative filters and interpolant coefficients.
    % Detection is done by seeing if any of the derivative specifier
    % arguments is longer than 1 char, this implies 2nd derivative needed.
    secondDeriv = false;    
    for n = 1:length(varargin)
        if length(varargin{n}) > 1
            secondDeriv = true;
            break
        end
    end
    
    if ~secondDeriv
        % 5 tap 1st derivative cofficients.  These are optimal if you are just
        % seeking the 1st deriavtives
        p = [0.037659  0.249153  0.426375  0.249153  0.037659];
        d1 =[0.109604  0.276691  0.000000 -0.276691 -0.109604];
    else         
        % 5-tap 2nd derivative coefficients. The associated 1st derivative
        % coefficients are not quite as optimal as the ones above but are
        % consistent with the 2nd derivative interpolator p and thus are
        % appropriate to use if you are after both 1st and 2nd derivatives.
        p  = [0.030320  0.249724  0.439911  0.249724  0.030320];
        d1 = [0.104550  0.292315  0.000000 -0.292315 -0.104550];
        d2 = [0.232905  0.002668 -0.471147  0.002668  0.232905];
    end

    % Compute derivatives.  Note that in the 1st call below MATLAB's conv2
    % function performs a 1D convolution down the columns using p then a 1D
    % convolution along the rows using d1. etc etc.
    
    for n = 1:length(varargin)
      if strcmpi('x', varargin{n})
          varargout{n} = conv2(p, d1, im, 'same');    
      elseif strcmpi('y', varargin{n})
          varargout{n} = conv2(d1, p, im, 'same');
      elseif strcmpi('xx', varargin{n})
          varargout{n} = conv2(p, d2, im, 'same');    
      elseif strcmpi('yy', varargin{n})
          varargout{n} = conv2(d2, p, im, 'same');
      elseif strcmpi('xy', varargin{n}) || strcmpi('yx', varargin{n})
          varargout{n} = conv2(d1, d1, im, 'same');
      else
          error('''%s'' is an unrecognized derivative option',varargin{n});
      end
    end
  
    

function n = normalise(im, reqmean, reqvar)
   
	if ~isa(im,'double'), im = double(im); end	
	im = im - mean(im(:));    
	im = im/std(im(:));      % Zero mean, unit std dev

	n = reqmean + im*sqrt(reqvar);
    
    
function smim = gaussfilt(im, sigma)
         
    sze = max(ceil(6*sigma), 1);
    if ~mod(sze,2)    % Ensure filter size is odd
        sze = sze+1;
    end
    
    h = fspecial('gaussian', [sze sze], sigma);

    % Apply filter to all image channels   
    smim = filter2(h, im);
    