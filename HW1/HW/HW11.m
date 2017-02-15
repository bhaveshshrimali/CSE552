%Nonlinear Finite Element Method : Fall 2016

%HW Assignment 1
%Date: 09/02/2016

clear all; close all; clc;
disp ('%-------------------------------------------------------------------------%')
disp ('%         CEE 576: Nonlinear Finite Element Method                 %')
disp ('%      HW Assignment #1  (Newton Raphson Code )              %')
disp ('%-------------------------------------------------------------------------%')

err_in = 0.9;
err = err_in;
c = 1;
k= 1.2;
count = 1;

while err > 0.036
    error(count,1) = err;
    err = err^1.2;
    count = count +1;
end
disp (['In case (a) The Newton Raphson takes ' num2str(count) ' steps to converge to an error of 0.036' ]);
% error(:,2) = (linspace(1,count-1))';
filename = 'error1.xlsx';
xlswrite(filename,error);

err = err_in;
count =1;

while err > 0.036
    error2(count,1) = err;
    err = 0.7*err;
    count = count+1;
end


disp (['In case (b) The Newton Raphson takes ' num2str(count) ' steps to converge to an error of 0.036' ]);
% error2(:,2) = (linspace(1,count-1))';
filename = 'error2.xlsx';
xlswrite(filename,error2);
