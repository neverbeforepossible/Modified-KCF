function result = GetPeakInPriorSearchWindow(...
    pPos, pMvVer, pMvHor,...
    pIm, pWindow_sz, pFeatures, pCell_size, pCos_window, pModel_xf, pModel_alphaf, pKernel )
    if(( pMvVer < 2 ) & ( pMvVer > -2 ) & ( pMvHor < 2 ) & ( pMvHor > -2 ))
        searchWinPosY = 0;
        searchWinPosX = 0;
    else
        pMvVer = typecast( pMvVer, 'single');
        pMvHor = typecast( pMvHor, 'single');
        slope = pMvVer / pMvHor;
        if( pMvHor > 0 )
            if( slope > 2.41 )
                searchWinPosY = -1;
                searchWinPosX = 0;
            elseif( slope > 0.41 )
                searchWinPosY = -1;
                searchWinPosX = 1;
            elseif( slope < -2.41 )
                searchWinPosY = 1;
                searchWinPosX = 0;                
            elseif( slope < -0.41 )
                searchWinPosY = 1;
                searchWinPosX = 1;
            else
                searchWinPosY = 0;
                searchWinPosX = 1;
            end
        else
            if( slope > 2.41 )
                searchWinPosY = 1;
                searchWinPosX = 0;
            elseif( slope > 0.41 )
                searchWinPosY = 1;
                searchWinPosX = -1;
            elseif( slope < -2.41 )
                searchWinPosY = -1;
                searchWinPosX = 0;                
            elseif( slope < -0.41 )
                searchWinPosY = -1;
                searchWinPosX = -1;
            else
                searchWinPosY = 0;
                searchWinPosX = -1;
            end
        end
    end
    
    virPos = [ pPos( 1, 1 ) + ( searchWinPosY * pWindow_sz( 1, 1 )), pPos( 1, 2 ) + ( searchWinPosX * pWindow_sz( 1, 2 ))];
    patch = get_subwindow( pIm, virPos, pWindow_sz );
    zf = fft2( get_features( patch, pFeatures, pCell_size, pCos_window ));
    kzf = gaussian_correlation( zf, pModel_xf, pKernel.sigma );
    response = real( ifft2( pModel_alphaf .* kzf ));
    [ vertDelta, horDelta ] = find( response == max( response(:)), 1 );
    
    virPosY = virPos(1);
    virPosX = virPos(2);
    peak = response( vertDelta, horDelta );
    
    result = [ virPosY, virPosX, vertDelta, horDelta, peak, searchWinPosY, searchWinPosX ];
end