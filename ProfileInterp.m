function [ OutMat  ] = ProfileInterp(ProfileMat, NumProfiles)
    %PROFILEINTERP Interpolates and generates fuselage profiles
    %   Usage:
    % OutMat = ProfileInter(ProfileMat, NumProfiles)
    %   Where:
    % ProfileMat is a matrix of size p x N*3 where p is the number of points
    % in one fuselage profile and N is the number of profiles being
    % interpolated between. Colums are in sets of 3, where the first column
    % is the list of x points in the profile, column 2 is the list of y
    % points in the profile, and column 3 is the location along the body
    % (nondimensionalized - i.e. between 0 and 1 inclusive) of the profile.
    % Only the first line of the third (sixth, ninth, etc) column is used.
    OutMat = zeros(length(ProfileMat(:,1)),NumProfiles);
    
    %%Fill Z locations (Non-Dimensionally, ranging from 0 to 1)
    for i = 1:NumProfiles
        OutMat(:,i*3) = (i-1)/(NumProfiles-1);
    end
    for i = 1:NumProfiles
        under = true;
        profile = 1; 
        while under == true
            if ProfileMat(1,profile*3+3) >= (i-1)/(NumProfiles-1);
                under = false;
            else
                profile = profile + 1;
            end
        end
        scale = 1-(OutMat(1,i*3)-ProfileMat(1,profile*3))/(ProfileMat(1,profile*3+3)-ProfileMat(1,profile*3));
        OutMat(:,i*3-2) = (ProfileMat(:,profile).*scale + ProfileMat(:,profile+3).*(1-scale));
        OutMat(:,i*3-1) = (ProfileMat(:,profile+1).*scale + ProfileMat(:,profile+4).*(1-scale));        
    end
    hold off
    plot3(OutMat(:,1), OutMat(:,2), OutMat(:,3));
    hold on 
    for i = 1:NumProfiles-1
        plot3(OutMat(:,i*3+1), OutMat(:,i*3+2), OutMat(:,i*3+3));
    end
    hold off
    
end

