function [positions, time] = tracker(video_path, img_files, pos, target_sz, ...
	padding, kernel, lambda, output_sigma_factor, interp_factor, cell_size, ...
	features, show_visualization)
%TRACKER Kernelized/Dual Correlation Filter (KCF/DCF) tracking.
%   This function implements the pipeline for tracking with the KCF (by
%   choosing a non-linear kernel) and DCF (by choosing a linear kernel).
%
%   It is meant to be called by the interface function RUN_TRACKER, which
%   sets up the parameters and loads the video information.
%
%   Parameters:
%     VIDEO_PATH is the location of the image files (must end with a slash
%      '/' or '\').
%     IMG_FILES is a cell array of image file names.
%     POS and TARGET_SZ are the initial position and size of the target
%      (both in format [rows, columns]).
%     PADDING is the additional tracked region, for context, relative to 
%      the target size.
%     KERNEL is a struct describing the kernel. The field TYPE must be one
%      of 'gaussian', 'polynomial' or 'linear'. The optional fields SIGMA,
%      POLY_A and POLY_B are the parameters for the Gaussian and Polynomial
%      kernels.
%     OUTPUT_SIGMA_FACTOR is the spatial bandwidth of the regression
%      target, relative to the target size.
%     INTERP_FACTOR is the adaptation rate of the tracker.
%     CELL_SIZE is the number of pixels per cell (must be 1 if using raw
%      pixels).
%     FEATURES is a struct describing the used features (see GET_FEATURES).
%     SHOW_VISUALIZATION will show an interactive video if set to true.
%
%   Outputs:
%    POSITIONS is an Nx2 matrix of target positions over time (in the
%     format [rows, columns]).
%    TIME is the tracker execution time, without video loading/rendering.
%
%   Joao F. Henriques, 2014


% ¡åconfigurations¡å %
    failureThreshold = 0.6;
    restartThreshold = 0.6;
    windowExpAmount = 1;
    visualizeSearchGrid = 0;
    
    sizeOfCorrPeaksList = 30;
    numOfCorrPeaksTaken = 0;
    corrPeaks = linspace( 0, 0, sizeOfCorrPeaksList );
    isTracking = 1;
% ¡ãconfigurations¡ã %

% ¡åconfigurations¡å %
    sizeOfAveNeighborCorrList = 30;
    numOfAveNeighborCorrTaken = 0;
    aveNeighborCorrs = linspace( 0, 0, sizeOfAveNeighborCorrList );
% ¡ãconfigurations¡ã %


    peakValues = [];
    maxNumOfMvs = 3;       % ¢Ð configure
    numOfValidMvs = 0;
    mvVert = [];
    mvHor = [];
    aveMvVert = 0;
    aveMvHor = 0;
    for idx = 1:maxNumOfMvs
        mvVert = [mvVert, 0]
        mvHor = [mvHor, 0]        
    end


	%if the target is large, lower the resolution, we don't need that much
	%detail
	resize_image = ( sqrt( prod( target_sz )) >= 100 );  %diagonal size >= threshold
	if resize_image
		pos = floor( pos / 2 );
		target_sz = floor( target_sz / 2 );
	end


	%window size, taking padding into account
	window_sz = floor( target_sz * ( 1 + padding ));
	
% 	%we could choose a size that is a power of two, for better FFT
% 	%performance. in practice it is slower, due to the larger window size.
% 	window_sz = 2 .^ nextpow2(window_sz);

	
	%create regression labels, gaussian shaped, with a bandwidth
	%proportional to target size
	output_sigma = sqrt( prod( target_sz )) * output_sigma_factor / cell_size;
	yf = fft2( gaussian_shaped_labels( output_sigma, floor( window_sz / cell_size )));
    
    
	%store pre-computed cosine window
	cos_window = hann( size( yf, 1 )) * hann( size( yf, 2 ))';	
	
	
	if show_visualization   %create video interface
		update_visualization = show_video( img_files, video_path, resize_image );
	end
	
	
	%note: variables ending with 'f' are in the Fourier domain.

	time = 0;  %to calculate FPS
	positions = zeros( numel( img_files ), 2 );  %to calculate precision

	for frame = 1:numel( img_files )
		%load image
		im = imread([video_path img_files{frame}]);
		if size( im, 3 ) > 1
			im = rgb2gray( im );
		end
		if resize_image
			im = imresize( im, 0.5 );
		end

		tic()

		if frame > 1
            if isTracking == 1
                               
        		%obtain a subwindow for detection at the position from last
    			%frame, and convert to Fourier domain (its size is unchanged)
    			patch = get_subwindow( im, pos, window_sz );
    			zf = fft2( get_features( patch, features, cell_size, cos_window ));
			
    			%calculate response of the classifier at all shifts
    			switch kernel.type
    			case 'gaussian'
                	kzf = gaussian_correlation( zf, model_xf, kernel.sigma );
            	case 'polynomial'
        			kzf = polynomial_correlation( zf, model_xf, kernel.poly_a, kernel.poly_b );
    			case 'linear'
    				kzf = linear_correlation( zf, model_xf );
                end
                response = real( ifft2( model_alphaf .* kzf ));  %equation for fast detection
            
                %target location is at the maximum response. we must take into
                %account the fact that, if the target doesn't move, the peak
                %will appear at the top-left corner, not at the center (this is
                %discussed in the paper). the responses wrap around cyclically.
                [vert_delta, horiz_delta] = find( response == max( response(:)), 1 );                
            else  % if isTracking == 1
                
                % In this method, we use always Gaussian correlation. We
                % omit 'switch'.
                                
                priorVirCorrPeak = GetPeakInPriorSearchWindow( pos, aveMvVert, aveMvHor, im, window_sz, features, cell_size, cos_window, model_xf, model_alphaf, kernel );
                fprintf('Prior Search Window Y: %d, Prior Search Window X: %d\n', priorVirCorrPeak(6), priorVirCorrPeak(7));
                maxVirCorrPeak = [ 0, 0, 0, 0, 0 ];
                if(( priorVirCorrPeak(5) / averCorrPeaks ) >= restartThreshold )
                    maxVirCorrPeak(1) = priorVirCorrPeak(1);
                    maxVirCorrPeak(2) = priorVirCorrPeak(2);
                    maxVirCorrPeak(3) = priorVirCorrPeak(3);
                    maxVirCorrPeak(4) = priorVirCorrPeak(4);
                    maxVirCorrPeak(5) = priorVirCorrPeak(5);
                else
                    virCorrPeaks = [];
                    for winBlkVer = -windowExpAmount:windowExpAmount
                        for winBlkHor = -windowExpAmount:windowExpAmount
                            if(( winBlkVer == priorVirCorrPeak( 6 )) & ( winBlkHor == priorVirCorrPeak( 7 )))
                               virCorrPeaks = [ virCorrPeaks; priorVirCorrPeak(1), priorVirCorrPeak(2), priorVirCorrPeak(3), priorVirCorrPeak(4), priorVirCorrPeak(5)];
                               continue;
                            end
                            virPos = [ pos( 1, 1 ) + ( winBlkVer * window_sz( 1, 1 )), pos( 1, 2 ) + ( winBlkHor * window_sz( 1, 2 ))];
                            patch = get_subwindow( im, virPos, window_sz );
                            zf = fft2( get_features( patch, features, cell_size, cos_window ));
                            kzf = gaussian_correlation( zf, model_xf, kernel.sigma );
                            response = real( ifft2( model_alphaf .* kzf ));
                            [ vert_delta, horiz_delta ] = find( response == max( response(:)), 1 );
                            virCorrPeaks = [ virCorrPeaks ; virPos, vert_delta, horiz_delta, response( vert_delta, horiz_delta )];
                        end
                    end                    
                    
                    for virCorrPeaksIdx = 1:(( 2 * windowExpAmount + 1 ) * ( 2 * windowExpAmount + 1 ))
                        if ( virCorrPeaks( virCorrPeaksIdx, 5 ) > maxVirCorrPeak( 1, 5 ))
                            maxVirCorrPeak = virCorrPeaks( virCorrPeaksIdx, : );
                        end
                    end
                end
                
            end  % if isTracking == 1
            
            
            if isTracking == 1
                aveNeighborCorr = GetAveNeighborCorr( response, horiz_delta, vert_delta, size( zf, 2 ), size( zf, 1 ));
                aveNeighborCorrs( mod( numOfAveNeighborCorrTaken, sizeOfAveNeighborCorrList ) + 1 ) = aveNeighborCorr;
                numOfAveNeighborCorrTaken = numOfAveNeighborCorrTaken + 1;
                if( numOfAveNeighborCorrTaken < sizeOfAveNeighborCorrList )
                    aveAveNeighborCorrs = sum( aveNeighborCorrs ) / numOfAveNeighborCorrTaken;
                else
                    aveAveNeighborCorrs = sum( aveNeighborCorrs ) / sizeOfAveNeighborCorrList;
                end
            end
            
            
            if isTracking == 1
                curCorrPeak = response( vert_delta, horiz_delta );
            else
                curCorrPeak = maxVirCorrPeak( 1, 5 );
            end
            
            
            % Gathering peak values.
            if isTracking == 1
                corrPeaks( mod( numOfCorrPeaksTaken, sizeOfCorrPeaksList ) + 1 ) = curCorrPeak;
                numOfCorrPeaksTaken = numOfCorrPeaksTaken + 1;
                if( numOfCorrPeaksTaken < sizeOfCorrPeaksList ) 
                    averCorrPeaks = sum( corrPeaks ) / numOfCorrPeaksTaken;
                else
                    averCorrPeaks = sum( corrPeaks ) / sizeOfCorrPeaksList;
                end
            end
            
            
            % Calculating standard deviation.
            sumOfPeakDevSq = 0;
            if( numOfCorrPeaksTaken < sizeOfCorrPeaksList )
                for peakIdx = 1:numOfCorrPeaksTaken
                    sumOfPeakDevSq = sumOfPeakDevSq + (( corrPeaks( peakIdx ) - averCorrPeaks ) * ( corrPeaks( peakIdx ) - averCorrPeaks ));
                end
                peakVariation = sumOfPeakDevSq / numOfCorrPeaksTaken;
            else
                for peakIdx = 1:sizeOfCorrPeaksList
                    sumOfPeakDevSq = sumOfPeakDevSq + (( corrPeaks( peakIdx ) - averCorrPeaks ) * ( corrPeaks( peakIdx ) - averCorrPeaks ));
                end
                peakVariation = sumOfPeakDevSq / sizeOfCorrPeaksList;
            end
            peakStdDev = sqrt( peakVariation );
            
            
            curPeakToAverPeak = curCorrPeak / averCorrPeaks;
            curPeakToAveAveNeighborCorr = curCorrPeak / aveAveNeighborCorrs;
            
            
            fprintf('curPeakToAveAveNeighborCorr: %f\n', curPeakToAveAveNeighborCorr );
            
            
            if 1
                if(( isTracking == 1 ) & (( curPeakToAverPeak < failureThreshold ) | ( curPeakToAveAveNeighborCorr < 0.6 )))
                    fprintf('============================== Stop Tracking ==============================\n');
                    isTracking = 0;
                    fprintf("Debug point\n");
                end
            end
                        
            
            if 1
                if ( isTracking == 0 ) & ( curPeakToAverPeak >= restartThreshold )
                    fprintf('============================== Restart Tracking ==============================\n');
                    pos( 1, 1 ) = maxVirCorrPeak( 1, 1 );
                    pos( 1, 2 ) = maxVirCorrPeak( 1, 2 );
                    vert_delta = maxVirCorrPeak( 1, 3 );
                    horiz_delta = maxVirCorrPeak( 1, 4 );
                    isTracking = 1;
                    fprintf("Debug point\n");
                end
            end
            
            
            fprintf('%d.  CurCorr : %f || AverCorr : %f || CurToAver : %f || AverAverNeighborCorr: %f || CurToAverAverNeighborCorr : %f\n', frame, curCorrPeak, averCorrPeaks, curCorrPeak / averCorrPeaks, aveAveNeighborCorrs, curCorrPeak / aveAveNeighborCorrs );
            
            
            if frame == 0
                fprintf('DEBUG \n', frame );
            end

            
            
			if vert_delta > size( zf, 1 ) / 2  %wrap around to negative half-space of vertical axis
				vert_delta = vert_delta - size( zf, 1 );
            end
            
			if horiz_delta > size(zf,2) / 2  %same for horizontal axis
				horiz_delta = horiz_delta - size(zf,2);
            end
            
            
            if( isTracking == 1 )
                mvVert( mod( frame, maxNumOfMvs ) + 1 ) = cell_size * ( vert_delta - 1 );
                mvHor( mod( frame, maxNumOfMvs ) + 1 ) = cell_size * ( horiz_delta - 1 );
                                
                if( numOfValidMvs < maxNumOfMvs )
                    numOfValidMvs = numOfValidMvs + 1;
                end
                
                aveMvVert = int32( sum( mvVert ) / numOfValidMvs );
                aveMvHor = int32( sum( mvHor ) / numOfValidMvs );                    
                
                if( frame == numel( img_files ))
                    fprintf('<JUNGSUP> Last Frame.');
                end
            end
            
            
            if isTracking == 1
                pos = pos + cell_size * [vert_delta - 1, horiz_delta - 1];
            end
        end     % frame > 1

        
		%obtain a subwindow for training at newly estimated target position
		patch = get_subwindow(im, pos, window_sz);
		xf = fft2(get_features(patch, features, cell_size, cos_window));
        
        
		%Kernel Ridge Regression, calculate alphas (in Fourier domain)
		switch kernel.type
		case 'gaussian'
			kf = gaussian_correlation( xf, xf, kernel.sigma );
		case 'polynomial'
			kf = polynomial_correlation(xf, xf, kernel.poly_a, kernel.poly_b);
		case 'linear'
			kf = linear_correlation(xf, xf);
		end
		alphaf = yf ./ (kf + lambda);   %equation for fast training

		if frame == 1  %first frame, train with a single image
			model_alphaf = alphaf;
			model_xf = xf;
		else
			%subsequent frames, interpolate model
            if isTracking == 1
                model_alphaf = (1 - interp_factor) * model_alphaf + interp_factor * alphaf;
                model_xf = (1 - interp_factor) * model_xf + interp_factor * xf;
            end
        end
        

		%save position and timing
		positions(frame,:) = pos;
		time = time + toc();
        

		%visualization
		if show_visualization
			box = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
			stop = update_visualization(frame, box);
			if stop, break, end  %user pressed Esc, stop early
            
            
            if ( visualizeSearchGrid == 1 ) & ( isTracking == 0 )
                basePosX = pos( 1, 2 ) - 0.5 * window_sz( 1, 2 );
                basePosY = pos( 1, 1 ) - 0.5 * window_sz( 1, 1 );
                for visSearchGridVer = -windowExpAmount:windowExpAmount
                    for visSearchGridHor = -windowExpAmount:windowExpAmount
                        box = [ basePosX + visSearchGridHor * window_sz( 1, 2 ), basePosY + visSearchGridVer * window_sz( 1, 1 ), window_sz( 1, 2 ), window_sz( 1, 1 )];
                        stop = update_visualization(frame, box);               
                    end
                end
            end
            
            
			drawnow
            
            % Image file dump during tracking process
            if 1
                outFilePath = sprintf("output/%05d.jpg", frame );
                saveas( gcf, outFilePath );
            end            
            
        end
        
        
        % Write estimated target position
        if 0
            estimatedXofTopLeft = pos(2) - ( target_sz(2) / 2 );
            estimatedYofTopLeft = pos(1) - ( target_sz(1) / 2 );
            estimatedWidth = target_sz(2);
            estimatedHeight = target_sz(1);
            
            estimatedXofTopLeft = round( estimatedXofTopLeft );
            estimatedYofTopLeft = round( estimatedYofTopLeft );
            estimatedWidth = round( estimatedWidth );
            estimatedHeight = round( estimatedHeight );            
            
            resultFilePath = sprintf('output/estimatedPosition.txt');
            fpWrite = fopen( resultFilePath, 'a');
            fprintf( fpWrite, '%d,%d,%d,%d\n', estimatedXofTopLeft, estimatedYofTopLeft, estimatedWidth, estimatedHeight );
            fclose( fpWrite );
        end
        
        
        % Drawing search window
        if 0
            if( frame > 1 )
                if( frame > 2 )
                    delete( searchWindowRect );
                end
                
                searchWindowX = pos(2) - ( window_sz(2) / 2 );
                searchWindowY = pos(1) - ( window_sz(1) / 2 );
                searchWindowRect = rectangle('Position', [searchWindowX searchWindowY window_sz(2) window_sz(1)]);
                searchWindowRect.LineWidth = 3;
            end
        end
    end     % frame = 1:numel( img_files )
    
    
	if resize_image
		positions = positions * 2;
	end
end

