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

i = 2;

load([cell2mat(folders(i)) '/voltage_simulink_2s.mat']);
v_simulink = v(2,:);

load([cell2mat(folders(i)) '/voltage_simulink_bms_2s.mat']);
v_simulink_bms = v(2,:);

v_m2 = voltage_predict(starts(i):stops(i));

v_IP_sim = voltage_IPmodel_simulink(starts(i):stops(i));

v_measured = voltage(starts(i):stops(i));

time = [1:numel(v_measured)]/6;
time_simbms = [1:numel(v_simulink_bms)]/30;

fig1 = figure(1);
plot(time(1:end-1),v_measured(1:end-1), 'LineWidth', 3);
hold on;
plot(time(1:end-1),v_IP_sim(1:end-1), '-.', 'LineWidth', 3);
hold on;
plot(time_simbms,v_simulink_bms, '-', 'LineWidth', 0.5, 'color', [23, 86, 44]/255);
hold on;
plot(time_simbms,v_simulink, ':', 'LineWidth', 1, 'color', [84, 54, 107]/255);


hold on;
plot([min(time), max(time)+1], [3.7, 3.7], '-', 'color', [0.4 0.4 0.4])
hold on;
plot([min(time), max(time)+1], [2, 2], '-', 'color', [0.4 0.4 0.4])

xlabel('Time (min.)')
ylabel('Voltage (V)')
legend('Measured', 'IP model', 'TDB model', 'TD model');
set(gca, 'FontSize', 15)
xlim([0, max(time)+1]);
ylim([1.4,3.9])

%plot the inset
fig2 = figure(2);

plot(time(1:end-1),v_measured(1:end-1), 'LineWidth', 3);
hold on;
plot(time(1:end-1),v_IP_sim(1:end-1), '-.', 'LineWidth', 3);
hold on;
plot(time_simbms,v_simulink_bms, '-', 'LineWidth', 0.5, 'color', [23, 86, 44]/255);
hold on;
plot(time_simbms,v_simulink, ':', 'LineWidth', 1, 'color', [84, 54, 107]/255);

hold on;
plot([min(time), max(time)+1], [3.7, 3.7], '-', 'color', [0.4 0.4 0.4])
hold on;
plot([min(time), max(time)+1], [2, 2], '-', 'color', [0.4 0.4 0.4])

set(gca, 'FontSize', 15)
xlim([56, max(time)+0.1]);
ylim([1.4,3])
set(gca, 'Xtick', []);
set(gca, 'Xticklabel', []);
set(gca, 'Ytick', []);
set(gca, 'Yticklabel', []);

%plot figure-in-figure
[h_m h_i]=inset(fig1,fig2);
