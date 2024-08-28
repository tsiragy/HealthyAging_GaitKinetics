%% clear workspace %% 
clear all 
close all
clc

%% Call the Data %% 
COM_Old = {};
COM_Young = {};
time = 1/120;


for i = 1:length(data_O);
    for k = 1:length(data_O{i}.markers);
        COM_Old{i}{k} = data_O{i}.markers{k}.CentreOfMass(:,1);
    end 
end 


COMv_O = {};
avg_COMv_O = {};

for i = 1:length(COM_Old);
    for k = 1:length(COM_Old{i});
        COMv_O{i}{k} = abs((diff(COM_Old{i}{k}(:,:)/time)));  
    end 
end 


for i = 1:length(COMv_O);
    for k = 1:length(COMv_O{i});
      avg_COMv_O{i}{k} = mean(COMv_O{i}{k});
    end 
end 

for i = 1:length(avg_COMv_O);
    for k = 1:length(avg_COMv_O{i});
    COM_Vel_O(:,i) = mean(avg_COMv_O{i}{k});
    end
end 

COM_Vel_O = (COM_Vel_O/1000)'; %convert mm/s to m/s


%% COM Velocity Young %% 

for i = 1:length(data_Y);
    for k = 1:length(data_Y{i}.markers);
        COM_Young{i}{k} = data_Y{i}.markers{k}.CentreOfMass(:,1);
    end 
end 


COMv_Y = {};
avg_COMv_Y = {};

for i = 1:length(COM_Young);
    for k = 1:length(COM_Young{i});
        COMv_Y{i}{k} = abs((diff(COM_Young{i}{k}(:,:)/time)));  
    end 
end 


for i = 1:length(COMv_Y);
    for k = 1:length(COMv_Y{i});
      avg_COMv_Y{i}{k} = mean(COMv_Y{i}{k});
    end 
end 

for i = 1:length(avg_COMv_Y);
    for k = 1:length(avg_COMv_Y{i});
    COM_Vel_Y(:,i) = mean(avg_COMv_Y{i}{k});
    end
end 

COM_Vel_Y = (COM_Vel_Y/1000)'; %convert mm/s to m/s

%% independent t-test %% 

[h,p] = ttest2(COM_Vel_O,COM_Vel_Y)
