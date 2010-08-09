function [ OutMat  ] = ProfileInterp(ProfileMat, NumProfiles)
    %PROFILEINTERP Summary of this function goes here
    %   Detailed explanation goes here
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

