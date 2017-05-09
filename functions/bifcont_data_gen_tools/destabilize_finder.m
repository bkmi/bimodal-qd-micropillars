function [ indDestable, foldNum, hopfNum ] = destabilize_finder( nunst )
% Given a list of nunst values, the destabilize_finder finds the indicies
% of bifurcations which go from stable to unstable. In addition, it will
% output the "number" of that bifurcation for fold or hopf.

% find which ones
indDestable = find(abs(diff(nunst == 0)));
indFold = find(abs(diff(nunst))==1);
indHopf = find(abs(diff(nunst))==2);

% get the number for each bif
[~,foldNum,~] = intersect(indFold,indDestable);
[~,hopfNum,~] = intersect(indHopf,indDestable);

end

