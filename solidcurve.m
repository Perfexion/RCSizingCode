function solidcurve(N)
%this function accepts a matrix "N" with columns of xyz data of
%the form n=[x1,y1,z1,x2,y2,z2, ... xn,yn,zn];
%writes it to an excel sheet and calls the import curve macro
%which draws a curve through xyz points into SolidWorks.

[m,n]=size(N);
t(:,1)=1:m;

Excel = actxserver ('Excel.Application');
Folder=pwd;
File='\xyzPathMacro.xls';
File = strcat(Folder, File);
 
if ~exist(File,'file') 
    ExcelWorkbook = Excel.workbooks.Add; 
    ExcelWorkbook.SaveAs(File,1); 
    ExcelWorkbook.Close(false); 
end 
wkbk = Excel.Workbooks;
invoke(Excel.Workbooks,'Open',File);

for i = 1:n/3
    j = 1:3:n-1;
    curve = N(:,j(i):j(i)+2);
    xlswrite1(File, t, 'Sheet1', (['A2:A' int2str(eval('m+1'))]));
    xlswrite1(File, curve, 'Sheet1', (['B2:D' int2str(eval('m+1'))])); 
    b=cell(10,n+1);
    xlswrite1(File, b, 'Sheet1', (['A' int2str(eval('m+2')) ':D' int2str(eval('m+102'))]));
    invoke(Excel.ActiveWorkbook,'Save');
    Excel.ExecuteExcel4Macro('!MakeCurve()');
end 
Excel.Quit 
Excel.delete 
clear Excel
