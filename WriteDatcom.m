function [] = WriteDatcom()
%%This function will write a datcom input file

    try
        AirDevilsOut = evalin('base', 'AirDevilsOut');
    catch exception
        error('No cell array was found in the workspace. Did you run the sizing module?');
        rethrow(exception)
    end
    [rows cols] = size(AirDevilsOut);
    
    fid = fopen('test.test','w');
    fprintf(fid, '%s\n', strcat('CASE ID AirDevilsOut ', datestr(clock)));
    fprintf(fid, '%s\n', strcat('$FLTCON NMACH=1.0,MACH(1)=0.1,NALT=1,ALT(1)=',num2str(AirDevilsOut{57,2}),'.0,NALPHA=20.0,ALSHD(1)=-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,WT=',num2str(AirDevilsOut{61,2}),',STMACH=0.95$'));
    fprintf(fid, '%s\n', strcat('$OPTINS SREF=', num2str(AirDevilsOut{6,2}),',CBARR=',num2str(AirDevilsOut{10,2}),',BLREF=',num2str(AirDevilsOut{7,2}),'$'));
    
    fprintf(fid, '%s\n', 'Next Case');
    fclose(fid);
    
end

