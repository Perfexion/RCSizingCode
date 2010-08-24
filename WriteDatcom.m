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
    NumProfiles = 20;
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
        p(i) = pi*2*r(i);
    end
    
    %% Print Output
    fid = fopen(strcat('test','.dcm'),'w');
    fprintf(fid, '%s\n', 'DIM FT');
    %Flight Conditions
    fprintf(fid, '%s\n', strcat(' $FLTCON LOOP=2.0,NMACH=1.0,MACH(1)=0.1,TR=1.0,'));
    fprintf(fid, '%s\n', strcat('  NALT=1,ALT(1)=',num2str(AirDevilsOut{57,2},'%6.2f'),',NALPHA=20.0,'));
    fprintf(fid, '%s\n', strcat('  ALSCHD(1)=-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,'));
    fprintf(fid, '%s\n', strcat('   -1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0$'));
    %Options
    fprintf(fid, '%s\n', strcat(' $OPTINS SREF=', num2str(AirDevilsOut{6,2},'%6.2f'),','));
    fprintf(fid, '%s\n', strcat('  CBARR=',num2str(AirDevilsOut{10,2},'%6.2f'),',BLREF=',num2str(AirDevilsOut{7,2},'%6.2f'),'$'));
    %Synths
    fprintf(fid, '%s\n', strcat(' $SYNTHS VERTUP=.TRUE., XCG=',num2str(AirDevilsOut{66,2},'%6.2f'),','));
    fprintf(fid, '%s\n', strcat('  ZCG=',num2str(AirDevilsOut{67,2},'%6.2f'),',XW=',num2str(AirDevilsOut{5,2}-1/4*AirDevilsOut{8,2},'%6.2f'),','));
    fprintf(fid, '%s\n', strcat('  ZW=',num2str(AirDevilsOut{63,2},'%6.2f'),',ALIW=',num2str(AirDevilsOut{62,2},'%6.2f'),',XH=',num2str(AirDevilsOut{21,2}-1/4*AirDevilsOut{24,2},'%6.2f'),','));
    fprintf(fid, '%s\n', strcat('  ZH=',num2str(AirDevilsOut{65,2},'%6.2f'),',ALIH=',num2str(AirDevilsOut{64,2},'%6.2f'),',XV=',num2str(AirDevilsOut{36,2}-1/4*AirDevilsOut{39,2},'%6.2f'),',ZV=0$'));
    
    %Body
    fprintf(fid, '%s', strcat(' $BODY NX=',num2str(NumProfiles,'%6.2f'),',ITYPE=1.0,'));
    fprintf(fid, '\n%s', strcat('  X(1)='));
    count = 0;
    for i = 1:NumProfiles
        count = count+1;
        if count > 5
            fprintf(fid, '\n    %s', strcat(num2str(x(i),'%6.2f'),','));
            count = 0;
        else
            fprintf(fid, '%s', strcat(num2str(x(i),'%6.2f'),','));
        end
    end
    count = 0;
    fprintf(fid, '\n%s',' ZU(1)=');
    for i = 1:NumProfiles
        count = count+1;
       if count > 5
           fprintf(fid, '\n    %s', strcat(num2str(zu(i),'%6.2f'),','));
           count = 0;
       else
           fprintf(fid, '%s', strcat(num2str(zu(i),'%6.2f'),','));
       end
    end
    count = 0;
    fprintf(fid, '\n%s',' ZL(1)=');
    for i = 1:NumProfiles
        count = count+1;
       if count > 5
           fprintf(fid, '\n    %s', strcat(num2str(zl(i),'%6.2f'),','));
           count = 0;
       else
           fprintf(fid, '%s', strcat(num2str(zl(i),'%6.2f'),',')); 
       end
    end
    count=0;
    fprintf(fid, '\n%s',' S(1)=');
    for i = 1:NumProfiles
        count=count+1;
       if count > 5
           fprintf(fid, '\n    %s', strcat(num2str(s(i),'%6.2f'),','));
           count=0;
       else
           fprintf(fid, '%s', strcat(num2str(s(i),'%6.2f'),','));
       end
    end
    count=0;
    fprintf(fid, '\n%s',' P(1)=');
    for i = 1:NumProfiles
        count=count+1; 
       if count > 5
           fprintf(fid, '\n    %s', strcat(num2str(p(i),'%6.2f'),','));
           count = 0;
       else
           fprintf(fid, '%s', strcat(num2str(p(i),'%6.2f'),','));
       end
    end
    count = 0;
    fprintf(fid, '\n%s', ' R(1)=');
    for i = 1:NumProfiles
        count = count + 1;
       if count > 5
           fprintf(fid, '\n    %s', strcat(num2str(r(i),'%6.2f'),','));
           count = 0;
       else
           fprintf(fid, '%s', strcat(num2str(r(i),'%6.2f'),','));
       end
    end
    count = 0;
    fprintf(fid, '%s\n', '$');
    %Wing
    %fprintf(fid,'%s\n', AirDevilsOut{16,2});
    if length(AirDevilsOut{16,2}) == 4
        fprintf(fid, '%s\n', strcat('NACA-W-4-',AirDevilsOut{16,2}));
    elseif length(AirDevilsOut{16,2}) == 5
        fprintf(fid, '%s\n', strcat('NACA-W-5-',AirDevilsOut{16,2}));
    else
        fprintf(fid, '%s\n', 'NACA-W-5-14412');
    end
    
    fprintf(fid,'%s\n', strcat(' $WGPLNF CHRDR=',num2str(AirDevilsOut{8,2},'%6.2f'),',CHRDTP=',num2str(AirDevilsOut{9,2},'%6.2f'),',SSPN=', num2str(AirDevilsOut{7,2}/2,'%6.2f'),','));
    fprintf(fid,'%s\n', strcat('  SSPNE=',num2str(AirDevilsOut{7,2}/2-AirDevilsOut{3,2}*.45,'%6.2f'),',SAVSI=',num2str(AirDevilsOut{13,2}*180/pi,'%6.2f'),','));
    fprintf(fid,'%s\n', strcat('  CHSTAT=0.25,TWISTA=0.0,DHDADI=',num2str(AirDevilsOut{14,2}*180/pi,'%6.2f'),',TYPE=1.0','$'));
    %HTail
    %fprintf(fid,'%s\n', AirDevilsOut{32,2});
    if length(AirDevilsOut{32,2}) == 4
        fprintf(fid, '%s\n', strcat('NACA-H-4-',AirDevilsOut{32,2}));
    elseif length(AirDevilsOut{32,2}) == 5
        fprintf(fid, '%s\n', strcat('NACA-H-5-',AirDevilsOut{32,2}));
    else
        fprintf(fid, '%s\n', 'NACA-H-4-0010');
    end
    fprintf(fid,'%s\n', strcat(' $HTPLNF CHRDR=',num2str(AirDevilsOut{24,2},'%6.2f'),',CHRDTP=',num2str(AirDevilsOut{25,2},'%6.2f'),',SSPN=', num2str(AirDevilsOut{23,2}/2,'%6.2f'),','));
    fprintf(fid,'%s\n', strcat('  SSPNE=',num2str(AirDevilsOut{23,2}/2-AirDevilsOut{3,2}*.2,'%6.2f'),',SAVSI=',num2str(AirDevilsOut{29,2}*180/pi,'%6.2f'),','));
    fprintf(fid,'%s\n', strcat('  CHSTAT=0.25,TWISTA=0.0,DHDADI=',num2str(AirDevilsOut{30,2}*180/pi,'%6.2f'),',TYPE=1.0','$'));

    %VTail
    if length(AirDevilsOut{46,2}) == 4
        fprintf(fid, '%s\n', strcat('NACA-V-4-',AirDevilsOut{46,2}));
    elseif length(AirDevilsOut{46,2}) == 5
        fprintf(fid, '%s\n', strcat('NACA-V-5-',AirDevilsOut{46,2}));
    else
        fprintf(fid, '%s\n', 'NACA-V-4-0010');
    end
    fprintf(fid,'%s\n', strcat(' $VTPLNF CHRDR=',num2str(AirDevilsOut{39,2},'%6.2f'),',CHRDTP=',num2str(AirDevilsOut{40,2},'%6.2f'),',SSPN=', num2str(AirDevilsOut{38,2},'%6.2f'),','));
    fprintf(fid, '%s\n', strcat('  SSPNE=',num2str(AirDevilsOut{38,2},'%6.2f'),',SAVSI=',num2str(AirDevilsOut{44,2}*180/pi,'%6.2f'),','));
    fprintf(fid, '%s\n', strcat('  CHSTAT=0.25,TYPE=1.0$'));

    fclose(fid);
    
    system('test.dcm');
    axis equal
    
end

