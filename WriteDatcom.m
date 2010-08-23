function [] = WriteDatcom()
%%This function will write a datcom input file

    %% Get Input
    try
        AirDevilsOut = evalin('base', 'AirDevilsOut');
    catch exception
        error('No cell array was found in the workspace. Did you run the sizing module?');
        rethrow(exception)
    end
    
    %% Generate or Input Body
    NumProfiles = 40;
    try 
        profin = evalin('base', 'BodyProfiles');
    catch exception
        z = circle([0,0], 100, .5);
        profin(:,1:2) = z(:, 1:2);
        profin(:,3) = 0;
        profin(:,4:5) = z(:, 1:2);
        profin(:,6) = 1;
    end
        profiles = profileinterp(profin,NumProfiles);
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
       profiles(:,i*3-2:i*3-1) = profiles(:,i*3-2:i*3-1).*(bodyprofile(i))*AirDevilsOut{2,2};
       profiles(:,i*3) = profiles(:,i*3)*AirDevilsOut{1,2};
    end
    
    %% Determine Datcom body Inputs at each profile
    matrix = zeros(1,NumProfiles);
    zu = matrix;
    zl = matrix; 
    p = matrix;
    x = matrix;
    for i = 1:NumProfiles
        x(i) = profiles(1, i*3);
        zu(i) = max(profiles(:,i*3-1));
        zl(i) = min(profiles(:,i*3-1));
        r(i) = (zu(i)-zl(i))/2;
        s(i) = pi*r(i)^2;
        p(i) = pi*2*r(1);
    end
    
    
    %% Print Output
    fid = fopen(strcat('test','.dat'),'w');
    fprintf(fid, '%s\n', strcat('CASE ID AirDevilsOut ', datestr(clock)));
    %Flight Conditions
    fprintf(fid, '%s\n', strcat('$FLTCON NMACH=1.0,MACH(1)=0.1,NALT=1,ALT(1)=',num2str(AirDevilsOut{57,2},'%6.2f'),',NALPHA=20.0,ALSHD(1)=-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,WT=',num2str(AirDevilsOut{61,2},'%6.2f'),',STMACH=0.95$'));
    %Options
    fprintf(fid, '%s\n', strcat('$OPTINS SREF=', num2str(AirDevilsOut{6,2},'%6.2f'),',CBARR=',num2str(AirDevilsOut{10,2},'%6.2f'),',BLREF=',num2str(AirDevilsOut{7,2},'%6.2f'),'$'));
    %Synths
    fprintf(fid, '%s\n', strcat('$SYNTHS XCG=',num2str(AirDevilsOut{66,2},'%6.2f'),',ZCG=',num2str(AirDevilsOut{67,2},'%6.2f'),',XW=',num2str(AirDevilsOut{5,2}-1/4*AirDevilsOut{8,2},'%6.2f'),',ZW=',num2str(AirDevilsOut{63,2},'%6.2f'),',ALIW=',num2str(AirDevilsOut{62,2},'%6.2f'),',XH=',num2str(AirDevilsOut{21,2}-1/4*AirDevilsOut{24},'%6.2f'),',ZH=',num2str(AirDevilsOut{65,2},'%6.2f'),',ALIH=',num2str(AirDevilsOut{64,2},'%6.2f'),',XV=',num2str(AirDevilsOut{36,2}-1/4*AirDevilsOut{39,2},'%6.2f'),',ZV=0$'));
    
    %Body
    fprintf(fid, '%s', strcat('$BODY NX=',num2str(NumProfiles,'%6.2f'),',ITYPE=1.0,X(1)='));
    for i = 1:NumProfiles
        fprintf(fid, '%s', strcat(num2str(x(i),'%6.2f'),','));
    end
    fprintf(fid, '\n%s','ZU(1)=');
    for i = 1:NumProfiles
       fprintf(fid, '%s', strcat(num2str(zu(i),'%6.2f'),',')); 
    end
    fprintf(fid, '\n%s','ZL(1)=');
    for i = 1:NumProfiles
       fprintf(fid, '%s', strcat(num2str(zl(i),'%6.2f'),',')); 
    end
    fprintf(fid, '\n%s','S(1)=');
    for i = 1:NumProfiles
       fprintf(fid, '%s', strcat(num2str(s(i),'%6.2f'),',')); 
    end
    fprintf(fid, '\n%s','P(1)=');
    for i = 1:NumProfiles
       fprintf(fid, '%s', strcat(num2str(p(i),'%6.2f'),',')); 
    end
    fprintf(fid, '\n%s','R(1)=');
    for i = 1:NumProfiles
       fprintf(fid, '%s', strcat(num2str(r(i),'%6.2f'),',')); 
    end
    fprintf(fid, '%s\n', '$');
    %Wing
    fprintf(fid,'%s\n', AirDevilsOut{16,2});
    fprintf(fid,'%s\n', strcat('$WGPLNF CHRDR=',num2str(AirDevilsOut{8,2},'%6.2f'),',CHRDTP=',num2str(AirDevilsOut{9,2},'%6.2f'),',SSPN=', num2str(AirDevilsOut{7,2}/2,'%6.2f'),',SSPNE=',num2str(AirDevilsOut{7,2}/2-AirDevilsOut{3,2}*.45,'%6.2f'),',SAVSI=',num2str(AirDevilsOut{13,2}*180/pi,'%6.2f'),',CHSTAT=0.25,TWISTA=0.0,DHDADI=',num2str(AirDevilsOut{14,2}*180/pi,'%6.2f'),',TYPE=1.0','$'));
    %HTail
    fprintf(fid,'%s\n', AirDevilsOut{32,2});
    fprintf(fid,'%s\n', strcat('$HTPLNF CHRDR=',num2str(AirDevilsOut{24,2},'%6.2f'),',CHRDTP=',num2str(AirDevilsOut{25,2},'%6.2f'),',SSPN=', num2str(AirDevilsOut{23,2}/2,'%6.2f'),',SSPNE=',num2str(AirDevilsOut{23,2}/2-AirDevilsOut{3,2}*.2,'%6.2f'),',SAVSI=',num2str(AirDevilsOut{29,2}*180/pi,'%6.2f'),',CHSTAT=0.25,TWISTA=0.0,DHDADI=',num2str(AirDevilsOut{30,2}*180/pi,'%6.2f'),',TYPE=1.0','$'));

    %VTail
    fprintf(fid,'%s\n', AirDevilsOut{32,2});
    fprintf(fid,'%s\n', strcat('$VTPLNF CHRDR=',num2str(AirDevilsOut{39,2},'%6.2f'),',CHRDTP=',num2str(AirDevilsOut{40,2},'%6.2f'),',SSPN=', num2str(AirDevilsOut{38,2},'%6.2f'),',SSPNE=',num2str(AirDevilsOut{38,2},'%6.2f'),',SAVSI=',num2str(AirDevilsOut{44,2}*180/pi,'%6.2f'),',CHSTAT=0.25,TYPE=1.0$'));
    
    fprintf(fid, '%s\n', 'NEXT CASE');
    fclose(fid);
    
end

