function [mapName] = TNC_CreateRBColormap(numValues,type)
% FUNCTION DETAILS: Simple utility to create a RWB style color map
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

switch lower(type)
    
    case 'rb'
        mapLength = numValues;
        halfMap = floor(numValues./2);
        incr = 1./halfMap;
        increaser = (0:incr:1-incr)';
        tmp  = [ increaser, increaser, ones(halfMap,1)];
        tmp2 = [ ones(halfMap,1) , 1-increaser , 1-increaser ];
        mapName = [tmp;tmp2];

    case 'cb'
        mapLength = numValues;
        halfMap = floor(numValues./2);
        incr = 1./halfMap;
        incr2 = 0.33./halfMap;
        increaser = (0:incr:1-incr)';
        increaser2 = (0:incr2:0.33-incr2)';        
        tmp  = [ increaser, 0.67+increaser2, ones(halfMap,1)];
        tmp2 = [ ones(halfMap,1) , 1-increaser , 1-increaser ];
        mapName = [tmp;tmp2];

    case 'wblack'
        incr = 1./numValues;
        increaser = (0:incr:1-incr)';
        tmp  = [ increaser, increaser, increaser ];
        mapName = tmp;
        
    case 'wred'
        incr = 1./numValues;
        increaser = (0:incr:1-incr)';
        tmp  = [ ones(numValues,1), increaser, increaser ];
        mapName = tmp;

    case 'wblue'
        incr = 0.75./numValues;
        increaser = (0.75-incr:-incr:0)';
        tmp  = [ increaser, increaser, ones(numValues,1)*0.75 ];
        mapName = tmp;

end
