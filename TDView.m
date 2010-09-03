function [] = TDView(solidFlag, wingAirfoil, vTailAirfoil, hTailAirfoil)
    format long;

    Aircraft = evalin('base', 'Aircraft');
    NumProfiles = 40;
    
    %% Body
    try 
        profin = evalin('base', 'BodyProfiles');
    catch exception
        z = circle([0,0], 100, .5);
        profin(:,1:2) = z(:, 1:2);
        profin(:,3) = 0;
        profin(:,4:5) = z(:, 1:2);
        profin(:,6) = 1;
    end
        profiles = ProfileInterp(profin,NumProfiles);
    try 
        bodyprofile = evalin('base', 'BodySpline');
    catch exception        
        bodyprofile = profiles(1,3:3:end);
        bodyprofile = sin(sqrt(bodyprofile).*pi);
    end
    if length(bodyprofile) ~= length(profiles(1,3:3:end))
       x = linspace(0,1,length(bodyprofile));
       bodyprofile = interp1(x, bodyprofile, profiles(1,3:3:end));
    end
    for i = 1:length(profiles(1,:))/3
       profiles(:,i*3-2:i*3-1) = profiles(:,i*3-2:i*3-1).*(bodyprofile(i))*Aircraft.FuselageDiameter;
       profiles(:,i*3) = profiles(:,i*3)*Aircraft.FuselageLength;
    end
    hold off
    plot3(profiles(:,1), profiles(:,2), profiles(:,3));
    hold on 
    
    for i = 1:length(profiles(1,:))/3-1
        plot3(profiles(:,i*3+1), profiles(:,i*3+2), profiles(:,i*3+3));
    end
    axis equal
    
    %Import to SolidWorks if solidFlag = 1.
    if solidFlag ==1
        solidcurve(profiles*1200);
    else
    end
    
    rotmat = zeros(2,2);
    %% Main Wing
    foilprofile = wingAirfoil;
    foilprofin = foilprofile(:,1:2);
    foilprofin(:,3) = 0;
    foilprofin(:,4:5) = foilprofile(:,1:2);
    foilprofin(:,6) = 1;
    foilprofiles = ProfileInterp(foilprofin,NumProfiles);
    rotmat = [cos(Aircraft.WingIncidence) -sin(Aircraft.WingIncidence);...
        sin(Aircraft.WingIncidence) cos(Aircraft.WingIncidence)];
    
    
    for i = 1:length(foilprofiles(1,:))/3
       foilprofiles(:,i*3-2:i*3-1) = foilprofiles(:,i*3-2:i*3-1).*(Aircraft.WingCR-(Aircraft.WingCR-Aircraft.WingCT).*foilprofiles(1,i*3)./foilprofiles(1,end)); %Scale Airfoil Profiles by chord
       foilprofiles(:,i*3) = foilprofiles(:,i*3)*Aircraft.WingB/2; %Scale the wing span
       % the following incidence angle rotation rotates about the leading
       % edge of the airfoil. It would make more aerodynamic sense to
       % rotate about the aerodynamic center of the airfoil. Maybe someone
       % can add an adjustment term to do this? 
       for j = 1:length(foilprofiles(:,1))
        foilprofiles(j,i*3-2:i*3-1) = foilprofiles(j,i*3-2:i*3-1)*rotmat; %Adjust for Incidence Angle
       end
       foilprofiles(:,i*3-2) = foilprofiles(:,i*3-2)+((Aircraft.WingCR-Aircraft.WingCT)/4/Aircraft.WingB*2+tan(Aircraft.WingQCSweep))*foilprofiles(1,i*3);%adjust for quarter chord sweep
       foilprofiles(:,i*3-2) = foilprofiles(:,i*3-2)+Aircraft.WingAC-.25*Aircraft.WingCR; %Adjust back to the located wing AC
       foilprofiles(:,i*3-1) = foilprofiles(:,i*3-1)+tan(Aircraft.WingDihedral)*foilprofiles(1,i*3); %Adjust for the dihadrial
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
    
    %Import to SolidWorks if solidFlag = 1.
    if solidFlag ==1
        solidcurve(wingprofiles*1200);
    else
    end
  
    %% Horizontal Tail
    hfoilprofile = hTailAirfoil;
    hfoilprofin = hfoilprofile(:,1:2);
    hfoilprofin(:,3) = 0;
    hfoilprofin(:,4:5) = hfoilprofile(:,1:2);
    hfoilprofin(:,6) = 1;
    hfoilprofiles = profileinterp(hfoilprofin,NumProfiles/4);
    rotmat = [cos(Aircraft.HTailIncidence) -sin(Aircraft.HTailIncidence);...
        sin(Aircraft.HTailIncidence) cos(Aircraft.HTailIncidence)]; 
    
    for i = 1:length(hfoilprofiles(1,:))/3
       hfoilprofiles(:,i*3-2:i*3-1) = hfoilprofiles(:,i*3-2:i*3-1).*(Aircraft.HTailCR-(Aircraft.HTailCR-Aircraft.HTailCT).*hfoilprofiles(1,i*3)./hfoilprofiles(1,end)); %Scale Airfoil Profiles by chord
       hfoilprofiles(:,i*3) = hfoilprofiles(:,i*3)*Aircraft.HTailB/2; %Scale the wing span
       for j = 1:length(foilprofiles(:,1))
        hfoilprofiles(j,i*3-2:i*3-1) = hfoilprofiles(j,i*3-2:i*3-1)*rotmat; %Adjust for Incidence Angle
       end
       hfoilprofiles(:,i*3-2) = hfoilprofiles(:,i*3-2)+((Aircraft.HTailCR-Aircraft.HTailCT)/4/Aircraft.HTailB*2+tan(Aircraft.HTailQCSweep))*hfoilprofiles(1,i*3);%adjust for quarter chord sweep
       hfoilprofiles(:,i*3-2) = hfoilprofiles(:,i*3-2)+Aircraft.HTailAC-.25*Aircraft.HTailCR; %Adjust back to the located wing AC
       hfoilprofiles(:,i*3-1) = hfoilprofiles(:,i*3-1)+tan(Aircraft.HTailDihedral)*hfoilprofiles(1,i*3); %Adjust for the dihadrial
    end
    hfoilprofiles2 = hfoilprofiles;
    hfoilprofiles2(:, 3:3:end) = hfoilprofiles(:,3:3:end)*-1;
    hfoilprofiles = [hfoilprofiles hfoilprofiles2];
   
    htailprofiles = hfoilprofiles;
    for i = 1:length(hfoilprofiles(1,:))/3
       htailprofiles(:,i*3-2) = hfoilprofiles(:,i*3);
       htailprofiles(:,i*3) = hfoilprofiles(:, i*3-2);
    end
    

    plot3(htailprofiles(:,1), htailprofiles(:,2), htailprofiles(:,3));
    hold on 
    for i = 1:length(htailprofiles(1,:))/3-1
        plot3(htailprofiles(:,i*3+1), htailprofiles(:,i*3+2), htailprofiles(:,i*3+3));
    end
    axis equal
    
    %Import to SolidWorks if solidFlag = 1.
    if solidFlag ==1
        solidcurve(htailprofiles*1200);
    else
    end

    %% Vertical Tail
    
    vfoilprofile = vTailAirfoil;
    vfoilprofin = vfoilprofile(:,1:2);
    vfoilprofin(:,3) = 0;
    vfoilprofin(:,4:5) = vfoilprofile(:,1:2);
    vfoilprofin(:,6) = 1;
    vfoilprofiles = profileinterp(vfoilprofin,NumProfiles/4);
    
    for i = 1:length(vfoilprofiles(1,:))/3
       vfoilprofiles(:,i*3-2:i*3-1) = vfoilprofiles(:,i*3-2:i*3-1).*(Aircraft.VTailCR-(Aircraft.VTailCR-Aircraft.VTailCT).*vfoilprofiles(1,i*3)./vfoilprofiles(1,end)); %Scale Airfoil Profiles by chord
       vfoilprofiles(:,i*3) = vfoilprofiles(:,i*3)*Aircraft.VTailB; %Scale the wing span
       vfoilprofiles(:,i*3-2) = vfoilprofiles(:,i*3-2)+((Aircraft.VTailCR-Aircraft.VTailCT)/4/Aircraft.VTailB+tan(Aircraft.VTailQCSweep))*vfoilprofiles(1,i*3); %adjust for quarter chord sweep
       vfoilprofiles(:,i*3-2) = vfoilprofiles(:,i*3-2)+Aircraft.VTailAC-.25*Aircraft.VTailCR; %Adjust back to the located wing AC
    end
   
    vtailprofiles = vfoilprofiles;
    for i = 1:length(vfoilprofiles(1,:))/3
       vtailprofiles(:,i*3-1) = vfoilprofiles(:,i*3);
       vtailprofiles(:,i*3-2) = vfoilprofiles(:, i*3-1);
       vtailprofiles(:,i*3) = vfoilprofiles(:, i*3-2);
    end
    

    plot3(vtailprofiles(:,1), vtailprofiles(:,2), vtailprofiles(:,3));
    hold on 
    for i = 1:length(vtailprofiles(1,:))/3-1
        plot3(vtailprofiles(:,i*3+1), vtailprofiles(:,i*3+2), vtailprofiles(:,i*3+3));
    end
    axis equal
    
    %Import to SolidWorks if solidFlag = 1.
    if solidFlag ==1
        solidcurve(vtailprofiles*1200);
    else
    end
    
end