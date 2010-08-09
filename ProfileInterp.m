function [ OutMat  ] = ProfileInterp(ProfileMat, NumProfiles)
%PROFILEINTERP Summary of this function goes here
%   Detailed explanation goes here
OutMat = zeros(length(ProfileMat(:,1),NumProfiles);
for i = 1:NumProfiles;
    OutMat(:,i*3) = i/NumProfiles;
    for j = 1:length(ProfileMat(:,1))
       OutMat( 
    end
end


end

