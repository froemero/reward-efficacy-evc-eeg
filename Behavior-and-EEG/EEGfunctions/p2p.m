function [ P2P ] = p2p( MAT, low, up, el, win )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% make subset of data choosing relevant time window
tempMAT = mean(MAT(el,low:up,:),1);% possibility to avearge across electrodes
tempMAT = reshape(tempMAT,[size(tempMAT,2) size(tempMAT,3)]);
tempMAT = tempMAT'; % make structure trials * timepoints

% get minimum
[MIN, INT] = min(tempMAT, [],2);
% 
tempMAT2 = zeros(1, win,size(MAT,3));
for n= 1:length (INT)
    
    tempMAT2(1,:,n) = mean(MAT(el,low+INT(n)-win:low + INT(n)-1,n),1);
    
end
  
tempMAT2 = reshape(tempMAT2,[size(tempMAT2,2) size(tempMAT2,3)]);
tempMAT2 = tempMAT2'; % make structure trials * timepoints

MAX= max(tempMAT2, [],2);

P2P = MIN-MAX;
% 
end

