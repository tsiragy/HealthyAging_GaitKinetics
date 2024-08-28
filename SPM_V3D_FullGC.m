%% clear workspace %% 
%clear all 
%close all
%clc


%% Call Data %%
%mat = dir('*.mat'); 

%% Gait Events for Left Trials with 1 Full Stride 
for i = 1:length(LON);
if length(LON{i}) == 2; 
    LON_t{i} = LON{i};
    LOFF_t{i} = LOFF{i};
else 
end
end 


%for i = 1:length(LOFF);
%if length(LOFF{i}) == 2; 
 %   LOFF_t{i} = LOFF{i};
%else 
%end
%end 



if exist('LON_t','var')
LON_id= find(~cellfun(@isempty,LON_t));
%LOFF_id = find(~cellfun(@isempty,LOFF_t));
LON_f = {};
%LOFF_f = {};
Left_Ankle_Moment_f = {};
Left_Knee_Moment_f = {};
Left_Hip_Moment_f = {};
Left_Ankle_Power_f = {};
Left_Knee_Power_f = {};
Left_Hip_Power_f = {};
for i = 1:length(LON_id);
    LON_f{i} = LON{LON_id(i)};
    Left_Ankle_Moment_f{i} = Left_Ankle_Moment{LON_id(i)};
    Left_Knee_Moment_f{i} = Left_Knee_Moment{LON_id(i)};
    Left_Hip_Moment_f{i} = Left_Hip_Moment{LON_id(i)};

    Left_Ankle_Power_f{i} = Left_Ankle_Power{LON_id(i)};
    Left_Knee_Power_f{i} = Left_Knee_Power{LON_id(i)};
    Left_Hip_Power_f{i} = Left_Hip_Power{LON_id(i)};
end 

%for i = 1:length(LOFF_id);
    % LOFF_f{i} = LOFF{LOFF_id(i)};
%end 
LHS_FP = round(cell2mat(cellfun(@(x) x*120,LON_f,'un',0))); %convert time (sec) to frames by multiplying each by 120 Vicon frame rate, for conversion to force plates multiply by 1200
%LTO_FP = round(cell2mat(cellfun(@(x) x*120,LOFF_f,'un',0)));
else 
end 





%% Gait Events for Right Trials with 1 Full Stride




for i = 1:length(RON);
if length(RON{i}) > 1; 
    RON_t{i} = RON{i};
else 
end
end 





if exist('RON_t','var')
RON_id= find(~cellfun(@isempty,RON_t));
%ROFF_id = find(~cellfun(@isempty,ROFF_t));
RON_f = {};
%ROFF_f = {};
Right_Ankle_Moment_f = {};
Right_Knee_Moment_f = {};
Right_Hip_Moment_f = {};
Right_Ankle_Power_f = {};
Right_Knee_Power_f = {};
Right_Hip_Power_f = {};
for i = 1:length(RON_id);
    RON_f{i} = RON{RON_id(i)};
    Right_Ankle_Moment_f{i} = Right_Ankle_Moment{RON_id(i)};
    Right_Knee_Moment_f{i} = Right_Knee_Moment{RON_id(i)};
    Right_Hip_Moment_f{i} = Right_Hip_Moment{RON_id(i)};

    Right_Ankle_Power_f{i} = Right_Ankle_Power{RON_id(i)};
    Right_Knee_Power_f{i} = Right_Knee_Power{RON_id(i)};
    Right_Hip_Power_f{i} = Right_Hip_Power{RON_id(i)};
end 
%for i = 1:length(ROFF_id)
   % ROFF_f{i} = ROFF{ROFF_id(i)};
%end 
RHS_FP = round(cell2mat(cellfun(@(x) x*120,RON_f,'un',0))); %convert time (sec) to frames by multiplying each by 120 Vicon frame rate, for conversion to force plates multiply by 1200
%RTO_FP = round(cell2mat(cellfun(@(x) x*120,ROFF_f,'un',0)));
else 
end 



%% Left Kinetics %% 

L_HipMom = {};
L_KneeMom = {};
L_AnkMom = {};

L_HipMom_t = {};
L_KneeMom_t = {};
L_AnkMom_t = {};




L_HipPwr_t = {};
L_KneePwr_t = {};
L_AnkPwr_t = {};


L_HipPwr = {};
L_KneePwr = {};
L_AnkPwr = {};

if exist('LON_t','var')
for i = 1:size(LHS_FP,2); %Split the Moments for Left Force Plate Hits
        L_HipMom_t{i} = Left_Hip_Moment_f{i}(LHS_FP(1,i):LHS_FP(2,i),:);
        L_KneeMom_t{i} = Left_Knee_Moment_f{i}(LHS_FP(1,i):LHS_FP(2,i),:);
        L_AnkMom_t{i} = Left_Ankle_Moment_f{i}(LHS_FP(1,i):LHS_FP(2,i),:);
end 

% Interpolate to 101 points for Left Variables %
for i = 1:size(LHS_FP,2);
        L_HipMom{i} = interpft(L_HipMom_t{i},101);
        L_KneeMom{i} = interpft(L_KneeMom_t{i},101);
        L_AnkMom{i} = interpft(L_AnkMom_t{i},101);
     
end 

for i = 1:size(LHS_FP,2); %Split the Moments for Left Force Plate Hits
        L_HipPwr_t{i} = Left_Hip_Power_f{i}(LHS_FP(1,i):LHS_FP(2,i),:);
        L_KneePwr_t{i} = Left_Knee_Power_f{i}(LHS_FP(1,i):LHS_FP(2,i),:);
        L_AnkPwr_t{i} = Left_Ankle_Power_f{i}(LHS_FP(1,i):LHS_FP(2,i),:);
     
end 
% Interpolate to 101 points for Left Variables %
for i = 1:size(LHS_FP,2);
        L_HipPwr{i} = interpft(L_HipPwr_t{i},101);
        L_KneePwr{i} = interpft(L_KneePwr_t{i},101);
        L_AnkPwr{i} = interpft(L_AnkPwr_t{i},101);
 
end 
else 
end 


%% Right Kinetics %% 

R_HipMom = {};
R_KneeMom = {};
R_AnkMom = {};

R_HipMom_t = {};
R_KneeMom_t = {};
R_AnkMom_t = {};




R_HipPwr_t = {};
R_KneePwr_t = {};
R_AnkPwr_t = {};


R_HipPwr = {};
R_KneePwr = {};
R_AnkPwr = {};

if exist('RON_t','var')
for i = 1:size(RHS_FP,2); %Split the Moments for Left Force Plate Hits
        R_HipMom_t{i} = Right_Hip_Moment_f{i}(RHS_FP(1,i):RHS_FP(2,i),:);
        R_KneeMom_t{i} = Right_Knee_Moment_f{i}(RHS_FP(1,i):RHS_FP(2,i),:);
        R_AnkMom_t{i} = Right_Ankle_Moment_f{i}(RHS_FP(1,i):RHS_FP(2,i),:);
end 
% Interpolate to 101 points for Left Variables %
for i = 1:size(RHS_FP,2);
        R_HipMom{i} = interpft(R_HipMom_t{i},101);
        R_KneeMom{i} = interpft(R_KneeMom_t{i},101);
        R_AnkMom{i} = interpft(R_AnkMom_t{i},101);
     
end 
for i = 1:size(RHS_FP,2); %Split the Moments for Left Force Plate Hits
        R_HipPwr_t{i} = Right_Hip_Power_f{i}(RHS_FP(1,i):RHS_FP(2,i),:);
        R_KneePwr_t{i} = Right_Knee_Power_f{i}(RHS_FP(1,i):RHS_FP(2,i),:);
        R_AnkPwr_t{i} = Right_Ankle_Power_f{i}(RHS_FP(1,i):RHS_FP(2,i),:);
     
end 
% Interpolate to 101 points for Left Variables %
for i = 1:size(RHS_FP,2);
        R_HipPwr{i} = interpft(R_HipPwr_t{i},101);
        R_KneePwr{i} = interpft(R_KneePwr_t{i},101);
        R_AnkPwr{i} = interpft(R_AnkPwr_t{i},101);
 
end 
else 
end 



%% Create Plots of Grand Mean for Left Foot %% 

avg_LAnkM_AP = {}; 
avg_LKneeM_AP = {}; 
avg_LHipM_AP = {}; 

avg_LAnkM_ML = {}; 
avg_LKneeM_ML = {}; 
avg_LHipM_ML = {}; 

avg_LAnkM_VT = {}; 
avg_LKneeM_VT = {}; 
avg_LHipM_VT = {}; 



if exist('LON_t','var')
for i = 1:length(L_AnkMom);
        avg_LAnkM_AP{i} = mean(L_AnkMom{i}(:,1),2);
        avg_LKneeM_AP{i} = mean(L_KneeMom{i}(:,1),2);
        avg_LHipM_AP{i} = mean(L_HipMom{i}(:,1),2);

        avg_LAnkM_ML{i} = mean(L_AnkMom{i}(:,2),2);
        avg_LKneeM_ML{i} = mean(L_KneeMom{i}(:,2),2);
        avg_LHipM_ML{i} = mean(L_HipMom{i}(:,2),2);

        avg_LAnkM_VT{i} = mean(L_AnkMom{i}(:,3),2);
        avg_LKneeM_VT{i} = mean(L_KneeMom{i}(:,3),2);
        avg_LHipM_VT{i} = mean(L_HipMom{i}(:,3),2);
end 
avg_LAnkM_AP = mean(cell2mat(avg_LAnkM_AP),2);
avg_LKneeM_AP = mean(cell2mat(avg_LKneeM_AP),2);
avg_LHipM_AP = mean(cell2mat(avg_LHipM_AP),2);
avg_LAnkM_ML = mean(cell2mat(avg_LAnkM_ML),2);
avg_LKneeM_ML = mean(cell2mat(avg_LKneeM_ML),2);
avg_LHipM_ML = mean(cell2mat(avg_LHipM_ML),2);
avg_LAnkM_VT = mean(cell2mat(avg_LAnkM_VT),2);
avg_LKneeM_VT = mean(cell2mat(avg_LKneeM_VT),2);
avg_LHipM_VT = mean(cell2mat(avg_LHipM_VT),2);
else 
end

%% Create Plots of Grand Mean for Right Foot %% 

avg_RAnkM_AP = {}; 
avg_RKneeM_AP = {}; 
avg_RHipM_AP = {}; 


avg_RAnkM_ML = {}; 
avg_RKneeM_ML = {}; 
avg_RHipM_ML = {}; 

avg_RAnkM_VT = {}; 
avg_RKneeM_VT = {}; 
avg_RHipM_VT = {}; 


if exist('RON_t','var')
for i = 1:length(R_AnkMom);
        avg_RAnkM_AP{i} = mean(R_AnkMom{i}(:,1),2);
        avg_RKneeM_AP{i} = mean(R_KneeMom{i}(:,1),2);
        avg_RHipM_AP{i} = mean(R_HipMom{i}(:,1),2);


        avg_RAnkM_ML{i} = mean(R_AnkMom{i}(:,2),2);
        avg_RKneeM_ML{i} = mean(R_KneeMom{i}(:,2),2);
        avg_RHipM_ML{i} = mean(R_HipMom{i}(:,2),2);

        avg_RAnkM_VT{i} = mean(R_AnkMom{i}(:,3),2);
        avg_RKneeM_VT{i} = mean(R_KneeMom{i}(:,3),2);
        avg_RHipM_VT{i} = mean(R_HipMom{i}(:,3),2);
end 
avg_RAnkM_AP = mean(cell2mat(avg_RAnkM_AP),2);
avg_RKneeM_AP = mean(cell2mat(avg_RKneeM_AP),2);
avg_RHipM_AP = mean(cell2mat(avg_RHipM_AP),2);
avg_RAnkM_ML = mean(cell2mat(avg_RAnkM_ML),2);
avg_RKneeM_ML = mean(cell2mat(avg_RKneeM_ML),2);
avg_RHipM_ML = mean(cell2mat(avg_RHipM_ML),2);
avg_RAnkM_VT = mean(cell2mat(avg_RAnkM_VT),2);
avg_RKneeM_VT = mean(cell2mat(avg_RKneeM_VT),2);
avg_RHipM_VT = mean(cell2mat(avg_RHipM_VT),2);
else 
end 


%% Trial Grant Mean (average of both sides) %% 


Ankle_Mom_t = [L_AnkMom, R_AnkMom];
Knee_Mom_t = [L_KneeMom, R_KneeMom];
Hip_Mom_t = [L_HipMom, R_HipMom];


Ankle_Pwr_t = [L_AnkPwr, R_AnkPwr];
Knee_Pwr_t = [L_KneePwr, R_KneePwr];
Hip_Pwr_t = [L_HipPwr, R_HipPwr];



%P = sort(randperm(numel(Ankle_Mom_t),3));

P = [1,2,3];

%x1 = Ankle_Mom(P(1));
%x2 = Ankle_Mom(P(2));
%x3 = Ankle_Mom(P(3));


for i = 1:length(P);
    Ankle_Mom{i} = Ankle_Mom_t{P(i)};
    Knee_Mom{i} = Knee_Mom_t{P(i)};
    Hip_Mom{i} = Hip_Mom_t{P(i)};


    Ankle_Pwr{i} = Ankle_Pwr_t{P(i)};
    Knee_Pwr{i} = Knee_Pwr_t{P(i)};
    Hip_Pwr{i} = Hip_Pwr_t{P(i)};

end 



avg_AnkM_AP = (Ankle_Mom{1}(:,1) + Ankle_Mom{2}(:,1) + Ankle_Mom{3}(:,1))/3;
avg_AnkM_ML = (Ankle_Mom{1}(:,2) + Ankle_Mom{2}(:,2) + Ankle_Mom{3}(:,2))/3;
avg_AnkM_VT = (Ankle_Mom{1}(:,3) + Ankle_Mom{2}(:,3) + Ankle_Mom{3}(:,3))/3;



avg_KneeM_AP = (Knee_Mom{1}(:,1) + Knee_Mom{2}(:,1) + Knee_Mom{3}(:,1))/3;
avg_KneeM_ML = (Knee_Mom{1}(:,2) + Knee_Mom{2}(:,2) + Knee_Mom{3}(:,2))/3;
avg_KneeM_VT = (Knee_Mom{1}(:,3) + Knee_Mom{2}(:,3) + Knee_Mom{3}(:,3))/3;


avg_HipM_AP = (Hip_Mom{1}(:,1) + Hip_Mom{2}(:,1) + Hip_Mom{3}(:,1))/3;
avg_HipM_ML = (Hip_Mom{1}(:,2) + Hip_Mom{2}(:,2) + Hip_Mom{3}(:,2))/3;
avg_HipM_VT = (Hip_Mom{1}(:,3) + Hip_Mom{2}(:,3) + Hip_Mom{3}(:,3))/3;




% Graph the Moments
subplot(3,1,1)
plot(avg_HipM_AP)
xlim([0 101])
subplot(3,1,2)
plot(avg_KneeM_AP)
xlim([0 101])
subplot(3,1,3)
plot(avg_AnkM_AP)
xlim([0 101])
%% Average Power Values in each plane %% 
avg_AnkP_AP = (Ankle_Pwr{1}(:,1) + Ankle_Pwr{2}(:,1) + Ankle_Pwr{3}(:,1))/3;
avg_AnkP_ML = (Ankle_Pwr{1}(:,2) + Ankle_Pwr{2}(:,2) + Ankle_Pwr{3}(:,2))/3;
avg_AnkP_VT = (Ankle_Pwr{1}(:,3) + Ankle_Pwr{2}(:,3) + Ankle_Pwr{3}(:,3))/3;



avg_KneeP_AP = (Knee_Pwr{1}(:,1) + Knee_Pwr{2}(:,1) + Knee_Pwr{3}(:,1))/3;
avg_KneeP_ML = (Knee_Pwr{1}(:,2) + Knee_Pwr{2}(:,2) + Knee_Pwr{3}(:,2))/3;
avg_KneeP_VT = (Knee_Pwr{1}(:,3) + Knee_Pwr{2}(:,3) + Knee_Pwr{3}(:,3))/3;


avg_HipP_AP = (Hip_Pwr{1}(:,1) + Hip_Pwr{2}(:,1) + Hip_Pwr{3}(:,1))/3;
avg_HipP_ML = (Hip_Pwr{1}(:,2) + Hip_Pwr{2}(:,2) + Hip_Pwr{3}(:,2))/3;
avg_HipP_VT = (Hip_Pwr{1}(:,3) + Hip_Pwr{2}(:,3) + Hip_Pwr{3}(:,3))/3;


%% Margin of Stability

difference_l = {};
squared_l = {};
total_sum_l = {};
Sqrt_dis_l = {};
l_l = {};
w_l = {};

COG_R_t = {};




difference_r = {};
squared_r = {};
total_sum_r = {};
Sqrt_dis_r = {};
l_r = {};
w_r = {};


if exist('RON_id','var')
for i = 1:length(RON_id);
    COG_R{i} = COG{RON_id(i)};

    for k = 1:length(RHS_FP(i));
        COG_R_t{i} = COG_R{i}(RHS_FP(1,k):RHS_FP(2,k),:);
    end 
    COG_Rint{i} = interpft(COG_R_t{i},101);
    

    


for i = 1:length(RON_id);
    RANK_t{i} = RANK{RON_id(i)};

    for k = 1:length(RHS_FP(i));
        RANK_tp{i} = RANK_t{i}(RHS_FP(1,k):RHS_FP(2,k),:);
    end 
    RANK_int{i} = interpft(RANK_tp{i},101);

end


end 


for i = 1:length(RANK_int);
    difference_r{i} = minus(RANK_int{i}(:,1:3),COG_Rint{i});
    squared_r{i} = difference_r{i}.^2; %square the differences 
    total_sum_r{i} = sum(squared_r{i},2); % sum the differences 
    Sqrt_dis_r{i} = sqrt(total_sum_r{i});
    l_r{i} = mean(Sqrt_dis_r{i}); %length of the inverted pendulum in meters based on pythagorem theorem 
    w_r{i} = sqrt(9.81/l_r{i}); %Calculate eigenfrequency of the inverted pendulum 
end 

else 
end 


if exist('LON_id','var')
for i = 1:length(LON_id);
    COG_L{i} = COG{LON_id(i)};
    LANK_t{i} = LANK{LON_id(i)};

    for k = 1:length(LHS_FP(i));
        COG_L_t{i} = COG_L{i}(LHS_FP(1,k):LHS_FP(2,k),:);
        LANK_tp{i} = LANK_t{i}(LHS_FP(1,k):LHS_FP(2,k),:);
    end 

    COG_Lint{i} = interpft(COG_L_t{i},101);
    LANK_int{i} = interpft(LANK_tp{i},101);


end 




for i = 1:length(LANK_int);
    difference_l{i} = minus(LANK_int{i}(:,1:3),COG_Lint{i});
    squared_l{i} = difference_l{i}.^2; %square the differences 
    total_sum_l{i} = sum(squared_l{i},2); % sum the differences 
    Sqrt_dis_l{i} = sqrt(total_sum_l{i});
    l_l{i} = mean(Sqrt_dis_l{i}); %length of the inverted pendulum in meters based on pythagorem theorem 
    w_l{i} = sqrt(9.81/l_l{i}); %Calculate eigenfrequency of the inverted pendulum 
end 


else 
end 




l = [l_l, l_r];
l = mean(cell2mat(l));
w = [w_l, w_r];
w = mean(cell2mat(w));


%% Parse the position and velocity of the COM

pCOM_ML = {};
pCOM_AP = {};

pCOM_ML_int = {};
pCOM_AP_int = {};
pCOM_V = {};

vCOM_ML_int = {};
vCOM_AP_int = {};
vCOM_V = {};

for i = 1:length(COG);
    pCOM_ML{i} = COG{i}(:,2);
    pCOM_AP{i} = COG{i}(:,1);
    pCOM_V{i} = COG{i}(:,3);

    vCOM_ML{i} = COG_VELOCITY{i}(:,2);
    vCOM_AP{i} = COG_VELOCITY{i}(:,1);
    vCOM_V{i} = COG_VELOCITY{i}(:,3);

end 


%% Calculate the Extrapolated Center of Mass Variable in AP and ML

AP_xCOM = {};
ML_xCOM = {};
MOS_AP_Left_f ={};
MOS_ML_Left_f ={};


MOS_AP_Right_f ={};
MOS_ML_Right_f ={};


Left_MOS_AntP = {};
Right_MOS_AntP = {};
Left_MOS_MedL = {};
Right_MOS_MedL = {};


%for i = 1:length(vCOM_AP);
     %   nanIndices_vCOM = isnan(vCOM_AP{i});
        % Replace NaN values with an empty array
     %   vCOM_AP{i}(nanIndices_vCOM) = [];
       % vCOM_AP_int{i} = interpft(vCOM_AP{i},101);

     %   nanIndices_vCOM_ml = isnan(vCOM_ML{i});
        % Replace NaN values with an empty array
     %   vCOM_ML{i}(nanIndices_vCOM_ml) = [];
       % vCOM_ML_int{i} = interpft(vCOM_ML{i},101);





      %  nanIndices_pCOM = isnan(pCOM_AP{i});
        % Replace NaN values with an empty array
       % pCOM_AP{i}(nanIndices_pCOM) = [];
       % pCOM_AP_int{i} = interpft(pCOM_AP{i},101);

       % nanIndices_pCOM_ml = isnan(pCOM_ML{i});
        % Replace NaN values with an empty array
       % pCOM_ML{i}(nanIndices_pCOM_ml) = [];
       % pCOM_ML_int{i} = interpft(pCOM_ML{i},101);
%end 









for i = 1:length(Left_COP);
Left_CP_AP{i} = Left_COP{i}(:,1);
Left_CP_ML{i} = Left_COP{i}(:,2);
end

           


if exist('LON_id', 'var')

for i = 1:length(LON_id);
    Left_CP_AP_t{i} = Left_CP_AP{LON_id(i)};
    Left_CP_ML_t{i} = Left_CP_ML{LON_id(i)};

    pCOM_ML_Left_t{i} = pCOM_ML{LON_id(i)};
    pCOM_AP_Left_t{i} = pCOM_AP{LON_id(i)};

    vCOM_ML_Left_t{i} = vCOM_ML{LON_id(i)};
    vCOM_AP_Left_t{i} = vCOM_AP{LON_id(i)};

end 
        for i = 1:size(LHS_FP,2);
        pCOM_ML_Left{i} = pCOM_ML_Left_t{i}(LHS_FP(1,i):LHS_FP(2,i),:);
        pCOM_ML_Left_int{i} = interpft(pCOM_ML_Left{i},101);
        pCOM_AP_Left{i} = pCOM_AP_Left_t{i}(LHS_FP(1,i):LHS_FP(2,i),:);
        pCOM_AP_Left_int{i} = interpft(pCOM_AP_Left{i},101);


        vCOM_ML_Left{i} = vCOM_ML_Left_t{i}(LHS_FP(1,i):LHS_FP(2,i),:);
        vCOM_ML_Left_int{i} = interpft(vCOM_ML_Left{i},101);
        vCOM_AP_Left{i} = vCOM_AP_Left_t{i}(LHS_FP(1,i):LHS_FP(2,i),:);
        vCOM_AP_Left_int{i} = interpft(vCOM_AP_Left{i},101);

        AP_xCOM_Left{i} = pCOM_AP_Left_int{i} + (vCOM_AP_Left_int{i}/w); %Calculate xCOM in the AP direction
        ML_xCOM_Left{i} = pCOM_ML_Left_int{i} + (vCOM_ML_Left_int{i}/w); %Calculate xCOM in the AP direction





        Left_CP_AP_temp{i} = Left_CP_AP_t{i}(LHS_FP(1,i):LHS_FP(2,i),:);
        Left_CP_ML_temp{i} = Left_CP_ML_t{i}(LHS_FP(1,i):LHS_FP(2,i),:);



        nanIndices_LCP = isnan(Left_CP_AP_temp{i});
        % Replace NaN values with an empty array
        Left_CP_AP_temp{i}(nanIndices_LCP) = [];
        Left_CP_ML_temp{i}(nanIndices_LCP) = [];
        end

else 
end 



if exist('LON_id', 'var')
    for i = 1:length(LON_id);

        Left_CP_AP_t{i} = Left_CP_AP{LON_id(i)};
        nanIndices_LCP = isnan(Left_CP_AP_t{i});
        % Replace NaN values with an empty array
        Left_CP_AP_t{i}(nanIndices_LCP) = [];


        Left_CP_ML_t{i} = Left_CP_ML{LON_id(i)};
        nanIndices_LCP = isnan(Left_CP_ML_t{i});
        % Replace NaN values with an empty array
        Left_CP_ML_t{i}(nanIndices_LCP) = [];
        
        

    end 
else 
end 



if exist('LON_id', 'var')
for i = 1:length(Left_CP_AP_t);

    Left_CP_AP_int{i} = interpft(Left_CP_AP_temp{i}(:,:),101);
    Left_CP_ML_int{i} = interpft(Left_CP_ML_temp{i}(:,:),101);

end 

for i = 1:length(Left_CP_AP_int);
    Left_MOS_AntP{i} = minus(Left_CP_AP_int{i},AP_xCOM_Left{i});
    Left_MOS_MedL{i} = minus(Left_CP_ML_int{i},ML_xCOM_Left{i});
end 
else 
end 









for i = 1:length(Right_COP);
Right_CP_AP{i} = Right_COP{i}(:,1);
Right_CP_ML{i} = Right_COP{i}(:,2);
end

   

if exist('RON_id', 'var')

for i = 1:length(RON_id);
    Right_CP_AP_t{i} = Right_CP_AP{RON_id(i)};
    Right_CP_ML_t{i} = Right_CP_ML{RON_id(i)};

    pCOM_ML_Right_t{i} = pCOM_ML{RON_id(i)};
    pCOM_AP_Right_t{i} = pCOM_AP{RON_id(i)};

    vCOM_ML_Right_t{i} = vCOM_ML{RON_id(i)};
    vCOM_AP_Right_t{i} = vCOM_AP{RON_id(i)};

end 

        for i = 1:size(RHS_FP,2);
        pCOM_ML_Right{i} = pCOM_ML_Right_t{i}(RHS_FP(1,i):RHS_FP(2,i),:);
        pCOM_ML_Right_int{i} = interpft(pCOM_ML_Right{i},101);
        pCOM_AP_Right{i} = pCOM_AP_Right_t{i}(RHS_FP(1,i):RHS_FP(2,i),:);
        pCOM_AP_Right_int{i} = interpft(pCOM_AP_Right{i},101);


        vCOM_ML_Right{i} = vCOM_ML_Right_t{i}(RHS_FP(1,i):RHS_FP(2,i),:);
        vCOM_ML_Right_int{i} = interpft(vCOM_ML_Right{i},101);
        vCOM_AP_Right{i} = vCOM_AP_Right_t{i}(RHS_FP(1,i):RHS_FP(2,i),:);
        vCOM_AP_Right_int{i} = interpft(vCOM_AP_Right{i},101);

        AP_xCOM_Right{i} = pCOM_AP_Right_int{i} + (vCOM_AP_Right_int{i}/w); %Calculate xCOM in the AP direction
        ML_xCOM_Right{i} = pCOM_ML_Right_int{i} + (vCOM_ML_Right_int{i}/w); %Calculate xCOM in the AP direction


        Right_CP_AP_temp{i} = Right_CP_AP_t{i}(RHS_FP(1,i):RHS_FP(2,i),:);
        Right_CP_ML_temp{i} = Right_CP_ML_t{i}(RHS_FP(1,i):RHS_FP(2,i),:);



        nanIndices_RCP = isnan(Right_CP_AP_temp{i});
        % Replace NaN values with an empty array
        Right_CP_AP_temp{i}(nanIndices_RCP) = [];
        Right_CP_ML_temp{i}(nanIndices_RCP) = [];

        end

else 
end 


     


if exist('RON_id', 'var')
for i = 1:length(Right_CP_AP_t);
    Right_CP_AP_int{i} = interpft(Right_CP_AP_temp{i}(:,:),101);
    Right_CP_ML_int{i} = interpft(Right_CP_ML_temp{i}(:,:),101);

end 

for i = 1:length(Right_CP_AP_int);
    Right_MOS_AntP{i} = minus(Right_CP_AP_int{i},AP_xCOM_Right{i});
    Right_MOS_MedL{i} = minus(Right_CP_ML_int{i},ML_xCOM_Right{i});
end 
else 
end 

for i = 1:length(Left_MOS_MedL);
    Left_MOS_MedL{i} = Left_MOS_MedL{i}.*-1;
end 

MOS_AntP_t = [Left_MOS_AntP, Right_MOS_AntP];
MOS_MedL_t = [Left_MOS_MedL, Right_MOS_MedL];



for i = 1:length(P);
    MOS_AntP{i} = MOS_AntP_t{P(i)};
    MOS_MedL{i} = MOS_MedL_t{P(i)};
end



MOS_AntP = mean(cell2mat(MOS_AntP),2);
MOS_MedL = mean(cell2mat(MOS_MedL),2);




%% Plot MOS and xCOM against COP %% 
% Loop through creating figures
for i = 1:length(Left_MOS_MedL);
    % Create a new figure
    figure;
    
    % Plot something (optional)
    subplot(3,1,1)
    plot(Left_MOS_MedL{i});
    subplot(3,1,2)
    plot(ML_xCOM_Left{i});
    subplot(3,1,3)
    plot(Left_CP_ML_int{i});


    % Customize the figure (optional)
    title(['Left ' num2str(i)]);
    %xlabel('X-axis');
    %ylabel('Y-axis');
    
    % You can add more plots, annotations, etc. here
    
    % Save the figure if needed (optional)
     saveas(gcf, ['left_' num2str(i) '.png']);
end






for i = 1:length(Right_MOS_MedL);
    % Create a new figure
    figure;
    
    % Plot something (optional)
    subplot(3,1,1)
    plot(Right_MOS_MedL{i});
    subplot(3,1,2)
    plot(ML_xCOM_Right{i});
    subplot(3,1,3)
    plot(Right_CP_ML_int{i});


    % Customize the figure (optional)
    title(['Right ' num2str(i)]);
    %xlabel('X-axis');
    %ylabel('Y-axis');
    
    % You can add more plots, annotations, etc. here
    
    % Save the figure if needed (optional)
     saveas(gcf, ['right_' num2str(i) '.png']);
end

%% Save the Data for SPM %% 

filename = "Old_S014.mat";
folderToSaveFileTo = uigetdir ("C:\Users\tsira\Documents\New Baseline Data\SPM Data\Visual 3D");
save (fullfile(folderToSaveFileTo, filename),"avg_HipM_AP", "avg_KneeM_AP", "avg_AnkM_AP", "avg_HipM_ML", "avg_KneeM_ML", "avg_AnkM_ML", "avg_HipM_VT", "avg_KneeM_VT", "avg_AnkM_VT", "avg_HipP_AP", "avg_KneeP_AP", "avg_AnkP_AP", "avg_HipP_ML", "avg_KneeP_ML", "avg_AnkP_ML", "avg_HipP_VT", "avg_KneeP_VT", "avg_AnkP_VT", "MOS_MedL", "MOS_AntP");




%% clear workspace %% 
clear all
close all
clc