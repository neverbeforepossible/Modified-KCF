function aveNeighborCorr = GetAveNeighborCorr( pRespMap, pPosX, pPosY, pRespMapWidth, pRespMapHeight )
    sumOfCorr = 0;
    for rowNum = ( pPosY - 2 ) : ( pPosY + 2 )
        for colNum = ( pPosX - 2 ) : ( pPosX + 2 )
            curX = colNum;
            curY = rowNum;
            if( curX < 1 )
                curX = curX + pRespMapWidth;                
            end
            if( curX > pRespMapWidth )
                curX = curX - pRespMapWidth;
            end
            if( curY < 1 )
                curY = curY + pRespMapHeight;                
            end
            if( curY > pRespMapHeight )
                curY = curY - pRespMapHeight;
            end
            
            if(( rowNum ~= pPosY ) | ( colNum ~= pPosX ))
                sumOfCorr = sumOfCorr + pRespMap( curY, curX );
            end
        end
    end
    aveNeighborCorr = sumOfCorr / 24;
end