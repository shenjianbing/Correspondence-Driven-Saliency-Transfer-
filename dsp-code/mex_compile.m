if ( ispc )
%     mex -g ./mexopts.bat SPNodeMatchMex.cpp SPNodeMatch.cpp DT.cpp
%     mex -g ./mexopts.bat SPBPMex.cpp SPBP.cpp SPGraph.cpp DT.cpp
%     mex -g ./mexopts.bat PixelMatchMex.cpp PixelMatch.cpp 
    mex   SPNodeMatchMex.cpp SPNodeMatch.cpp DT.cpp
    mex   SPBPMex.cpp SPBP.cpp SPGraph.cpp DT.cpp
    mex   PixelMatchMex.cpp PixelMatch.cpp 
elseif ( isunix )
    mex -f ./gccopts.sh SPNodeMatchMex.cpp SPNodeMatch.cpp DT.cpp
    mex -f ./gccopts.sh SPBPMex.cpp SPBP.cpp SPGraph.cpp DT.cpp
    mex -f ./gccopts.sh PixelMatchMex.cpp PixelMatch.cpp 
else
    disp('NOT SUPPORTED');
end
    