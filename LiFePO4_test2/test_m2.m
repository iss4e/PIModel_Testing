% test the model

addpath('/u1/fkazhamiaka/Models/PIModel_Testing');

model = IPmodel_noiteration;

model.v_charging_filename='v_curve_charging_lifepo4.csv'; % CSV File containing constant-current charging voltage curves
model.v_discharging_filename='v_curve_discharging_lifepo4.csv'; % CSV File containing constant-current discharging voltage curves
model.nominal_capacity=1.1; % Nominal capacity (Ampere-hours)
model.R_i=0.05; % Internal impedance (Ohms)
model.max_charging_current=4.4; % Maximum charging current (Amperes)
model.max_discharging_current=11; % Maximum discharging current (Amperes)
model.initial_energy_content=1; % Initial energy content (Wh)
model.time_step = 1/360; % simulation time-step duration (hours)

load('test_current_input.mat') % current

load('voltage_measured.mat') % voltage


power = current.*voltage;
apply_power = zeros(1,numel(power));

v_modelled = zeros(1,numel(power));
i_modelled = zeros(1,numel(power));
b_modelled = zeros(1,numel(power));


for i=1:numel(power)
    
    VIB = zeros(1,3);
    apply_power(i) = power(i);
    
    % if discharging, ensure we discharge fully so that we can test
    % charging from empty battery.
    if (apply_power(i) < 0)
        %apply_power(i) = apply_power(i)*1.05;
    end
    success = 0;
    % try applying power. If exception is raised, adjust power accordingly
    % to prevent violating constraints of the battery model.
    if (i==17466)
        i = i;
    end
    while(~success)
        try

            VIB = model.step(apply_power(i));
            
            % if we get to this point, the model has not thrown an
            % exception, so we call it a success and exit the loop
            success = 1;
        catch exception

            if (strcmp(exception.identifier, 'Battery:Overcharge'))
                
                causeException = exception.cause{1};
                % exception tells us approximately what the maximum power that can be used to
                % charge the battery without violating the energy content constraint.
                max_apply_power = causeException.message;
                apply_power(i) = str2num(max_apply_power);
            
            elseif (strcmp(exception.identifier, 'Battery:Undercharge'))
                
                causeException = exception.cause{1};
                % exception tells us approximately what the maximum power that can be used to
                % discharge the battery without violating the energy content constraint.
                max_apply_power = causeException.message;
                apply_power(i) = str2num(max_apply_power);
                
            elseif (strcmp(exception.identifier, 'Battery:ChargingOrDischargingRate'))
                
                causeException = exception.cause{1};
                % exception tells us approximately what the maximum power that can be used to
                % charge or discharge the battery (whichever is being attempted) without violating the rate constraint.
                max_apply_power = causeException.message;
                apply_power(i) = str2num(max_apply_power);
            elseif (strcmp(exception.identifier, 'Battery:NoIntersection'))
                
                causeException = exception.cause{1};
                % exception tells us approximately what the maximum power that can be used to
                % discharge the battery without violating the energy content constraint.
                max_apply_power = causeException.message;
                apply_power(i) = str2num(max_apply_power);
            else 
                rethrow(exception);
            end
        end
    end
    
    v_modelled(i) = VIB(1);
    i_modelled(i) = VIB(2);
    b_modelled(i) = VIB(3);
    
end

figure;
plot(voltage)
hold on;
plot(v_modelled)

xlabel('Time')
ylabel('Voltage')
set(gca, 'FontSize', 15)

figure;
plot(current)
hold on;
plot(i_modelled)

xlabel('Time')
ylabel('Current')
set(gca, 'FontSize', 15)

save('v_PI_noiterate.mat', 'v_modelled');