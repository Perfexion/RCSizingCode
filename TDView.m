% function [] = TDView()
%     

    NumProfiles = 40;
    
    %% Body
    z = circle([0,0], 100, .5);
    profin(:,1:2) = z(:, 1:2);
    profin(:,3) = 0;
    profin(:,4:5) = z(:, 1:2);
    profin(:,6) = 1;
    profiles = profileinterp(profin,NumProfiles);
    
    
    bodyprofile = profiles(1,3:3:end);
    bodyprofile = sin(bodyprofile.*pi);
    for i = 1:length(profiles(1,:))/3
       profiles(:,i*3-2:i*3-1) = profiles(:,i*3-2:i*3-1).*(bodyprofile(i))*AirDevilsOut{2,2};
       profiles(:,i*3) = profiles(:,i*3)*AirDevilsOut{1,2};
    end
    hold off
    plot3(profiles(:,1), profiles(:,2), profiles(:,3));
    hold on 
    for i = 1:length(profiles(1,:))/3-1
        plot3(profiles(:,i*3+1), profiles(:,i*3+2), profiles(:,i*3+3));
    end
    axis equal
    
    %% Main Wing
    foilprofile = AirfoilBuilder('NACA 0010', 100);
    foilprofin = foilprofile(:,1:2);
    foilprofin(:,3) = 0;
    foilprofin(:,4:5) = foilprofile(:,1:2);
    foilprofin(:,6) = 1;
    foilprofiles = profileinterp(foilprofin,NumProfiles);
    
    for i = 1:length(foilprofiles(1,:))/3
       foilprofiles(:,i*3-2:i*3-1) = foilprofiles(:,i*3-2:i*3-1).*(AirDevilsOut{8,2}-(AirDevilsOut{8,2}-AirDevilsOut{9,2}).*foilprofiles(1,i*3)./foilprofiles(1,end)); %Scale Airfoil Profiles by chord
       foilprofiles(:,i*3) = foilprofiles(:,i*3)*AirDevilsOut{7,2}/2; %Scale the wing span
       foilprofiles(:,i*3-2) = foilprofiles(:,i*3-2)+(AirDevilsOut{8,2}-AirDevilsOut{9,2})./AirDevilsOut{7,2}.*foilprofiles(1,i*3); %adjust for leading edge sweep
       foilprofiles(:,i*3-2) = foilprofiles(:,i*3-2)+(foilprofiles(1,i*3)*(AirDevilsOut{13,2}-atan(AirDevilsOut{8,2}-AirDevilsOut{9,2})./AirDevilsOut{7,2})/2); %adjust for quarter chord sweep
       foilprofiles(:,i*3-2) = foilprofiles(:,i*3-2)+AirDevilsOut{5,2}-.25*AirDevilsOut{8,2}; %Adjust back to the located wing AC
    end
    foilprofiles2 = foilprofiles;
    foilprofiles2(:, 3:3:end) = foilprofiles(:,3:3:end)*-1;
    foilprofiles = [foilprofiles foilprofiles2];
   
    wingprofiles = foilprofiles; 
    for i = 1:length(foilprofiles(1,:))/3
       wingprofiles(:,i*3-2) = foilprofiles(:,i*3);
       wingprofiles(:,i*3) = foilprofiles(:, i*3-2);
    end
    

    plot3(wingprofiles(:,1), wingprofiles(:,2), wingprofiles(:,3));
    hold on 
    for i = 1:length(wingprofiles(1,:))/3-1
        plot3(wingprofiles(:,i*3+1), wingprofiles(:,i*3+2), wingprofiles(:,i*3+3));
    end
    axis equal
  
    %% Horizontal Tail
    hfoilprofile = AirfoilBuilder('NACA 0010', 100);
    hfoilprofin = hfoilprofile(:,1:2);
    hfoilprofin(:,3) = 0;
    hfoilprofin(:,4:5) = hfoilprofile(:,1:2);
    hfoilprofin(:,6) = 1;
    hfoilprofiles = profileinterp(hfoilprofin,NumProfiles);
    
    for i = 1:length(hfoilprofiles(1,:))/3
       hfoilprofiles(:,i*3-2:i*3-1) = hfoilprofiles(:,i*3-2:i*3-1).*(AirDevilsOut{23,2}-(AirDevilsOut{23,2}-AirDevilsOut{24,2}).*hfoilprofiles(1,i*3)./hfoilprofiles(1,end)); %Scale Airfoil Profiles by chord
       hfoilprofiles(:,i*3) = hfoilprofiles(:,i*3)*AirDevilsOut{22,2}/2; %Scale the wing span
       hfoilprofiles(:,i*3-2) = hfoilprofiles(:,i*3-2)+(AirDevilsOut{23,2}-AirDevilsOut{24,2})./AirDevilsOut{22,2}.*hfoilprofiles(1,i*3); %adjust for leading edge sweep
       hfoilprofiles(:,i*3-2) = hfoilprofiles(:,i*3-2)+(hfoilprofiles(1,i*3)*(AirDevilsOut{28,2}-atan(AirDevilsOut{23,2}-AirDevilsOut{24,2})./AirDevilsOut{22,2})/2); %adjust for quarter chord sweep
       hfoilprofiles(:,i*3-2) = hfoilprofiles(:,i*3-2)+AirDevilsOut{20,2}-.25*AirDevilsOut{23,2}; %Adjust back to the located wing AC
    end
    hfoilprofiles2 = hfoilprofiles;
    hfoilprofiles2(:, 3:3:end) = hfoilprofiles(:,3:3:end)*-1;
    hfoilprofiles = [hfoilprofiles hfoilprofiles2];
   
    htailprofiles = hfoilprofiles;
    for i = 1:length(hfoilprofiles(1,:))/3
       htailrofiles(:,i*3-2) = hfoilprofiles(:,i*3);
       htailprofiles(:,i*3) = hfoilprofiles(:, i*3-2);
    end
    

    plot3(htailprofiles(:,1), htailprofiles(:,2), htailprofiles(:,3));
    hold on 
    for i = 1:length(htailprofiles(1,:))/3-1
        plot3(htailprofiles(:,i*3+1), htailprofiles(:,i*3+2), htailprofiles(:,i*3+3));
    end
    axis equal
    
% end