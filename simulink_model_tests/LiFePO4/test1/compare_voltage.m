%compare simulink to model 2 and measured voltage.

num_tests = 7;

% index into current file of start point in each test
starts = [423 1516 2278 2869 3403 3907 4395];

% index into current file for end point in each test
stops = [1088 1855 2450 2985 3490 3977 4431];

% names of destination folders for each test

folders = {'C05', 'C1', 'C2', 'C3', 'C4', 'C5', 'C10'};

% currents corresponding to each test
curs = [-0.55, -1.1, -2.2, -3.3, -4.4, -5.5, -11];

load('voltage_m2.mat'); % voltage_predict
voltage_predict = voltage_predict(3280:end);
load('v_IPmodel_simulink_R05.mat'); % v_modelled
voltage_IPmodel_simulink = v_modelled(3280:end);
load('voltage_measured.mat'); % voltage
voltage = voltage(3280:end);
load('test_current_formatted.mat'); % newarray
current = newarray(2,:);

MAVE_m2 = zeros(1,num_tests);
MAVE_simulink = zeros(1,num_tests);
MAVE_simulink_bms = zeros(1,num_tests);
MAVE_IP_simulink = zeros(1,num_tests);

for i=1:num_tests

    load([cell2mat(folders(i)) '/voltage_simulink.mat']);
    v_simulink = v(2,:);
    load([cell2mat(folders(i)) '/voltage_simulink_bms.mat']);
    v_simulink_bms = v(2,:);
    
    v_m2 = voltage_predict(starts(i):stops(i));
    
    v_IP_sim = voltage_IPmodel_simulink(starts(i):stops(i));
    
    v_measured = voltage(starts(i):stops(i));

    c_measured = current(starts(i):stops(i));
    
%     figure;
%     plot(v_measured);
%     hold on;
%     plot(v_m2);
%     hold on;
%     plot(v_simulink);
%     hold on;
%     plot(v_simulink_bms)
%     hold on;
%     plot(v_IP_sim);
%     
%     xlabel('Time')
%     ylabel('Voltage')
%     legend('Measured', 'IP model', 'TD model', 'TDB model', 'IP simulink')
%     set(gca, 'FontSize', 15)
    
    % check only those currents which correspond to the current being
    % tested
    cur_tested = curs(i);
    
    residual_m2 = zeros(1,length(v_measured));
    residual_simulink = zeros(1,length(v_measured));
    residual_simulink_bms = zeros(1,length(v_measured));
    residual_IPmodel_simulink = zeros(1,length(v_measured));
    num_points = 0;
    
    for j=1:length(v_measured)
        
        if (abs(c_measured(j)-cur_tested) < 0.05)
            
            residual_m2(j) = abs(v_measured(j)-v_m2(j));
            residual_simulink(j) = abs(v_measured(j)-v_simulink(j));
            residual_simulink_bms(j) = abs(v_measured(j)-v_simulink_bms(j));
            residual_IPmodel_simulink(j) = abs(v_measured(j)-v_IP_sim(j));
            num_points = num_points+1;
            
        end
            
    end
    
    MAVE_m2(i) = sum(residual_m2)/num_points;
    MAVE_simulink(i) = sum(residual_simulink)/num_points;
    MAVE_simulink_bms(i) = sum(residual_simulink_bms)/num_points;
    MAVE_IP_simulink(i) = sum(residual_IPmodel_simulink)/num_points;
%     figure;
%     plot(residual_m2, 'linewidth', 2);
%     hold on;
%     plot(residual_simulink, 'linewidth', 2);
% 
%     legend('Model 2', 'Simulink')
%     xlabel('Time');
%     ylabel('Voltage Residual (V)');
%     %xlim([0,40100])
%     set(gca, 'FontSize', 15)
    
end    

figure;
%plot(curs, MAVE_m2, '*b');
%hold on;
plot(curs, MAVE_simulink, 'xr');
hold on;
plot(curs, MAVE_simulink_bms, 'ok');
hold on;
plot(curs, MAVE_IP_simulink, '*b');

legend('TD model', 'TDB model', 'IP model')
xlabel('Current');
ylabel('MAVE');
set(gca, 'FontSize', 15)

save('MAVE_lifepo4_discharging', 'MAVE_m2', 'MAVE_simulink', 'curs');

% load(['voltage_simulink_bms.mat']);
% v_simulink_bms = v(2,:);

% figure;
% plot(v(1,600:end), v_simulink_bms(600:end), 'LineWidth', 2)
% hold on;
% plot([v(1,600), v(1,end)], [2,2], '--')
% xlabel('Time (sec)')
% ylabel('Voltage')
% set(gca, 'FontSize', 15)
% xlim([v(1,600), v(1,end)]);
