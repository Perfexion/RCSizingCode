%this code will hopefully be used to optimize our plane, right now
%it gets the values of ideal lift coefficient from the datcom output.

%read datcom output
text = fileread('test.out');
cl = regexp(text, '\W(?<=IDEAL LIFT COEFFICIENT =    .)\d\d\d\d\d', 'match');

cl_ideal_wing=str2double(cl(1));
cl_ideal_htail=str2double(cl(2));
cl_ideal_vtail=str2double(cl(3));
