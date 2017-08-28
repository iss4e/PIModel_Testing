classdef IPmodel < matlab.System & matlab.system.mixin.Propagates
    % Integrated Power-based battery model: Matlab System block implementation. Note that voltage curves should be specified in three columns representing the following: C-rate, Ampere-hour content (Ah in Ah), Battery/Cell Voltage (in Volts). Rows corresponding to the same C-rate should appear in contiguous chronological order.

    % Public, tunable properties
    properties
        
    end

    properties(Nontunable)
        v_charging_filename='v_curve_charging_LTO.csv'; % CSV File containing constant-current charging voltage curves
        v_discharging_filename='v_curve_discharging_LTO.csv'; % CSV File containing constant-current discharging voltage curves
        nominal_capacity=30; % Nominal capacity (Ampere-hours)
        R_i=0.002; % Internal impedance (Ohms)
        max_charging_current=150; % Maximum charging current (Amperes)
        max_discharging_current=150; % Maximum discharging current (Amperes)
        initial_energy_content=41.3; % Initial energy content (Wh)
        time_step=1; % simulation time-step duration (hours)
        StringSetProperty='Throw BMS exceptions'; % BMS option
    end

    properties(Hidden, Constant)
        StringSetPropertySet = matlab.system.StringSet({'Throw BMS exceptions', 'No exceptions, apply BMS power suggestion automatically'});
    end
        
    properties(DiscreteState)
        V; % voltage
        I; % current
        b; % energy content
    end
    
    % Pre-computed constants
    properties(Access = public)
        a_1;
        a_2;
        M_function;
        eff_c;
        eff_d;
        alpha_c;
        alpha_d;
        Vnom_c;
        Vnom_d;
        v_prev;
        b_prev;
        tolerance;
        V_min;
        V_max;
        throw_exception;
        b_grid;
        c_grid;
        v_grid;
        
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            
            % Read the voltage files, and split them according to C-rate.
%            coder.extrinsic('-sync:on','csvread');

            charging_voltages = csvread(obj.v_charging_filename);
            num_charging_points = double(size(charging_voltages, 1));

            discharging_voltages = csvread(obj.v_discharging_filename);
            num_discharging_points = double(size(discharging_voltages, 1));

            % Step 1: clean charging and discharging data by removing duplicates and
            % sorting.

            unique_charging_rates = [];
            unique_charging_rate_indices = [];

            unique_discharging_rates = [];
            unique_discharging_rate_indices = [];

            for i=1:num_charging_points
                % identify unique c-rates
                if ~any(charging_voltages(i,1)==unique_charging_rates)
                    unique_charging_rates = [unique_charging_rates, charging_voltages(i)];
                    unique_charging_rate_indices = [unique_charging_rate_indices, i];
                end
            end
            num_unique_charging_rates = length(unique_charging_rates);

            % sort the values and remove duplicates in the Ah dimension
            new_charging_voltages = [];
            for i=1:num_unique_charging_rates

                start_index = unique_charging_rate_indices(i);

                end_index = [];
                if (i+1 > num_unique_charging_rates)
                    end_index = length(charging_voltages);
                else
                    end_index = unique_charging_rate_indices(i+1)-1;
                end

                voltage_curve = charging_voltages(start_index:end_index,3);
                [ah_curve, order1, order2] = unique(charging_voltages(start_index:end_index,2));
                voltage_curve = voltage_curve(order1);

                new_charging_voltages = [new_charging_voltages;[ones(length(ah_curve),1)*unique_charging_rates(i),ah_curve, voltage_curve]];

            end

            charging_voltages = new_charging_voltages;
            num_charging_points = size(charging_voltages, 1);

            for i=1:num_discharging_points
                % identify unique c-rates
                if ~any(discharging_voltages(i,1)==unique_discharging_rates)
                    unique_discharging_rates = [unique_discharging_rates, discharging_voltages(i)];
                    unique_discharging_rate_indices = [unique_discharging_rate_indices, i];
                end
            end
            num_unique_discharging_rates = length(unique_discharging_rates);

            % sort the values and remove duplicates in the Ah dimension
            new_discharging_voltages = [];
            for i=1:num_unique_discharging_rates

                start_index = unique_discharging_rate_indices(i);

                end_index = [];
                if (i+1 > num_unique_discharging_rates)
                    end_index = length(discharging_voltages);
                else
                    end_index = unique_discharging_rate_indices(i+1)-1;
                end

                voltage_curve = discharging_voltages(start_index:end_index,3);
                [ah_curve, order1, order2] = unique(discharging_voltages(start_index:end_index,2));
                voltage_curve = voltage_curve(order1);

                new_discharging_voltages = [new_discharging_voltages;[ones(length(ah_curve),1)*unique_discharging_rates(i),ah_curve, voltage_curve]];

            end

            discharging_voltages = new_discharging_voltages;
            num_discharging_points = size(discharging_voltages, 1);

            % Step 2: redo the calculation of unique indices and rates on the cleaned charging
            % and discharging data

            unique_charging_rates = [];
            unique_charging_rate_indices = [];

            unique_discharging_rates = [];
            unique_discharging_rate_indices = [];

            for i=1:num_charging_points
                % identify unique c-rates
                if ~any(charging_voltages(i,1)==unique_charging_rates)
                    unique_charging_rates = [unique_charging_rates, charging_voltages(i)];
                    unique_charging_rate_indices = [unique_charging_rate_indices, i];
                end
            end
            num_unique_charging_rates = length(unique_charging_rates);

            for i=1:num_discharging_points
                % identify unique c-rates
                if ~any(discharging_voltages(i,1)==unique_discharging_rates)
                    unique_discharging_rates = [unique_discharging_rates, discharging_voltages(i)];
                    unique_discharging_rate_indices = [unique_discharging_rate_indices, i];
                end
            end
            num_unique_discharging_rates = length(unique_discharging_rates);


            % Step 3: get nominal voltages corresponding to each charging and
            % discharging rate. This is done by calculating the mean of each
            % voltage curve. The voltage curve is interpolated evenly
            % before calculating the mean, to avoid sample bias.

            obj.Vnom_c = zeros(2, num_unique_charging_rates);
            obj.Vnom_d = zeros(2, num_unique_discharging_rates);

            for i=1:num_unique_charging_rates

                start_index = unique_charging_rate_indices(i);

                end_index = [];
                if (i+1 > num_unique_charging_rates)
                    end_index = length(charging_voltages);
                else
                    end_index = unique_charging_rate_indices(i+1)-1;
                end

                voltage_curve = charging_voltages(start_index:end_index,3);
                ah_curve = charging_voltages(start_index:end_index,2);

                min_ah = min(ah_curve);
                max_ah = max(ah_curve);

                interpolation_distance = (max_ah-min_ah)/100;

                obj.Vnom_c(1,i) = mean(interp1(ah_curve, voltage_curve, min_ah:interpolation_distance:max_ah));
                obj.Vnom_c(2,i) = unique_charging_rates(i)*obj.nominal_capacity;
            end

            for i=1:num_unique_discharging_rates

                start_index = unique_discharging_rate_indices(i);

                end_index = [];
                if (i+1 > num_unique_discharging_rates)
                    end_index = length(discharging_voltages);
                else
                    end_index = unique_discharging_rate_indices(i+1)-1;
                end

                % get the voltage and Ah curves, removing any duplicates in the Ah
                % dimension.
                voltage_curve = discharging_voltages(start_index:end_index,3);
                ah_curve = discharging_voltages(start_index:end_index,2);

                min_ah = min(ah_curve);
                max_ah = max(ah_curve);

                interpolation_distance = (max_ah-min_ah)/100;

                obj.Vnom_d(1,i) = mean(interp1(ah_curve, voltage_curve, min_ah:interpolation_distance:max_ah));
                obj.Vnom_d(2,i) = -unique_discharging_rates(i)*obj.nominal_capacity;
            end

            % Step 4: calculate the efficiency

            obj.eff_c = zeros(2, num_unique_charging_rates);
            obj.eff_d = zeros(2, num_unique_discharging_rates);

            for i=1:num_unique_charging_rates
                obj.eff_c(1,i) = 1 - ((obj.R_i*(unique_charging_rates(i)*obj.nominal_capacity))/obj.Vnom_c(1,i));
                obj.eff_c(2,i) = unique_charging_rates(i)*obj.nominal_capacity;
            end

            for i=1:num_unique_discharging_rates
                obj.eff_d(1,i) = 1 - ((obj.R_i*(unique_discharging_rates(i)*obj.nominal_capacity))/obj.Vnom_d(1,i));
                obj.eff_d(2,i) = -unique_discharging_rates(i)*obj.nominal_capacity;
            end


            % Step 5: calculate energy content. This is done by taking the riemann
            % sum to approximate the area under the voltage curve, and
            % taking the efficiency into account. At the same time,
            % calculate the energy content limits, which is the maximum
            % energy we could charge with or the energy left in the battery
            % when discharging.

            E_c = zeros(num_charging_points,1);
            E_d = zeros(num_discharging_points,1);

            obj.a_2 = zeros(2, num_unique_charging_rates);
            obj.a_1 = zeros(2, num_unique_discharging_rates);

            figure;
            
            for i=1:num_unique_charging_rates

                start_index = unique_charging_rate_indices(i);

                end_index = [];
                if (i+1 > num_unique_charging_rates)
                    end_index = length(charging_voltages);
                else
                    end_index = unique_charging_rate_indices(i+1)-1;
                end

                E_c(start_index) = 0;
                
                % "stretch" the ah_difference to make up for voltage
                % curves that start at a non-zero Ah count. 
                stretch_factor = 1 + (charging_voltages(start_index,2)/charging_voltages(end_index,2));
                
                for j=(start_index+1):end_index;
                                      
                    ah_difference = (charging_voltages(j,2)-charging_voltages(j-1,2))*stretch_factor;
                    E_c(j) = E_c(j-1) + charging_voltages(j,3)*ah_difference*obj.eff_c(1,i);

                end
                
%                 plot([0, E_c(start_index:end_index)'], [1.5, charging_voltages(start_index:end_index,3)'], 'LineWidth', 2);
%                 if (i < num_unique_charging_rates)
%                     hold on;
%                 else
%                     hold on;
%                     plot([0,72.5], [1.5, 1.5], '--k')
%                     hold on;
%                     plot([0,72.5], [2.7, 2.7], '--k')
%                     legend('0.1C', '0.5C', '1C', '2C', '3C', '4C', '5C')
%                     xlabel('Energy Content (Wh)')
%                     ylabel('Voltage (V)')
%                     set(gca, 'FontSize', 15)
%                     xlim([0,72.5])
%                 end
                
                obj.a_2(1,i) = E_c(end_index);
                obj.a_2(2,i) = unique_charging_rates(i)*obj.nominal_capacity;

            end

            for i=1:num_unique_discharging_rates

                start_index = unique_discharging_rate_indices(i);

                end_index = [];
                if (i+1 > num_unique_discharging_rates)
                    end_index = length(discharging_voltages);
                else
                    end_index = unique_discharging_rate_indices(i+1)-1;
                end

                E_d(start_index) = 0;

                for j=(start_index+1):end_index;

                    ah_difference = discharging_voltages(j,2)-discharging_voltages(j-1,2);
                    E_d(j) = E_d(j-1) + discharging_voltages(j,3)*ah_difference/obj.eff_d(1,i);

                end

            end

            % since the battery was being discharged, need to shift
            % energy content so that it represents how much energy is
            % left in the battery. We use the largest difference in energy
            % content to undertand the maximum energy that can be obtained.

            max_content = max(E_d);
            for i=1:num_unique_discharging_rates

                start_index = unique_discharging_rate_indices(i);

                end_index = [];
                if (i+1 > num_unique_discharging_rates)
                    end_index = length(discharging_voltages);
                else
                    end_index = unique_discharging_rate_indices(i+1)-1;
                end

                max_rate_content = max(E_d(start_index:end_index));

                for j=start_index:end_index
                    E_d(j) = E_d(j) + (max_content - max_rate_content);
                end

                obj.a_1(1,i) = max_content - max_rate_content;
                obj.a_1(2,i) = -unique_discharging_rates(i)*obj.nominal_capacity;

            end
            
            
            % Step 6: set maximum charging and discharging current
            
            obj.alpha_c = obj.max_charging_current;
            obj.alpha_d = -obj.max_discharging_current;
            
            
            % Step 7: set up the scattered interpolants (M function)

            all_currents = vertcat(charging_voltages(:,1), -1*discharging_voltages(:,1))*obj.nominal_capacity;
            all_bks = vertcat(E_c, E_d);
            all_vs = vertcat(charging_voltages(:,3), discharging_voltages(:,3));

            M_data = scatteredInterpolant(all_currents, all_bks, all_vs, 'linear', 'linear');
            %obj.M_function = scatteredInterpolant(all_currents, all_bks, all_vs, 'linear', 'linear');

            % V_min and V_max can be extracted from the voltage curves
            obj.V_min = min(all_vs);
            obj.V_max = max(all_vs);
            
            max_b = max(all_bks);
            min_b = min(all_bks);
            
            max_curr = max(all_currents);
            min_curr = min(all_currents);
            
            % calculate the INTerpolation RESolution
            b_int_res = (max_b - min_b)/100;
            curr_int_res = (max_curr - min_curr)/100;
            
            b_range = min_b:b_int_res:(max_b);
            curr_range = min_curr:curr_int_res:max_curr;
            
            
            [curr_grid,b_grid] = ndgrid(curr_range,b_range);
            
            
            v_grid = zeros(length(curr_range), length(b_range));
            
            for i=1:length(curr_range)
                for j=1:length(b_range)
                    % use M_function to get voltage estimates for the grid
                    % defined by the current and energy content. Limit the
                    % voltage estimates to be between V_min and V_max.
                    v_grid(i,j) = M_data(curr_grid(i,j), b_grid(i,j));
                end
            end
            
            obj.b_grid = b_grid;
            obj.v_grid = v_grid;
            obj.c_grid = curr_grid;
            
            obj.M_function = griddedInterpolant(curr_grid,b_grid,v_grid);
%             
%             test_bs = [72:-1:20];
%             
%             test_curs = ones(1,53)*(-100);
%             
%             test_vs = obj.M_function(test_curs,test_bs);
            

            % plot the M function;
%             figure;
%             scatter3(all_currents, all_bks, all_vs);
%             xlabel('Current (A)')
%             ylabel('Energy Content (Wh)')
%             zlabel('Voltage (V)')
%             
%             hold on;
%             scatter3(reshape(curr_grid, [1,numel(curr_grid)]), reshape(b_grid, [1,numel(b_grid)]), reshape(v_grid, [1,numel(v_grid)]));
            
%             hold on;
%             scatter3(test_curs,test_bs,test_vs)

            % Step 8: specify remaining constants
            
            % initialize the previous voltage value to a reasonable
            % estimate.
            obj.v_prev = obj.Vnom_c(1,1);
            
            % set the initial energy content
            obj.b_prev = obj.initial_energy_content;
            
            % tolerance that determines voltage conversion. Fixed at the
            % following value:
            obj.tolerance = 0.005; 
            
            % check if we want to throw exceptions or not.
            obj.throw_exception = false;
            if (strcmp(obj.StringSetProperty, 'Throw BMS exceptions'))
                obj.throw_exception = true;
            end
            
        end

        function y = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states. u is the input power.
            
            [curr_estimate, v_estimate, b_estimate, a_1_val, a_2_val, eff_c_val, eff_d_val] = iterateVoltage(obj,u);
 
            % check that current does not exceed limits
            if (curr_estimate > obj.alpha_c) || (curr_estimate < obj.alpha_d)
                
                msgID = 'Battery:ChargingOrDischargingRate';
                msg = 'Charging or discharging power is too high';
                baseException = MException(msgID,msg);
                
                if (curr_estimate > 0)
                    
                    % use binary search to find power whose current estimate is
                    % lower than alpha_c.
                    L = 0;
                    R = u;
                    power_resolution = (R-L)/50;
                    power_range = L:power_resolution:R;
                    max_pow = 0;
                    max_curr = 0;
                    curr_estimate = 0;
                    b_estimate = obj.b_prev;
                    v_estimate = obj.v_prev;

                    for p_index = 1:numel(power_range)

                        power_prime = power_range(p_index);
                        
                        [curr_prime, v_prime, b_prime, a_1_prime, a_2_prime, eff_c_prime, eff_d_prime] = iterateVoltage(obj,power_prime);

                        if ((curr_prime < obj.alpha_c) && (abs(curr_prime) > abs(max_curr)))
                            max_pow = power_prime;
                            max_curr = curr_prime;
                            curr_estimate = curr_prime;
                            b_estimate = b_prime;
                            v_estimate = v_prime;
                        end
                        
                    end

                    msgID2 = 'Battery:maxChargingPower';
                    msg2 = num2str(max_pow);
                    causeException = MException(msgID2, msg2);
                    baseException = addCause(baseException,causeException);

                else
                    
                    % Search to find the power corresponding to the max current
                    % estimate that is lower (in magnitude) than alpha_d.
                    L = u;
                    R = 0;
                    power_resolution = (R-L)/50;
                    power_range = L:power_resolution:R;
                    max_pow = 0;
                    max_curr = 0;
                    curr_estimate = 0;
                    b_estimate = obj.b_prev;
                    v_estimate = obj.v_prev;

                    for p_index = 1:numel(power_range)

                        power_prime = power_range(p_index);
                        
                        [curr_prime, v_prime, b_prime, a_1_prime, a_2_prime, eff_c_prime, eff_d_prime] = iterateVoltage(obj,power_prime);

                        if ((curr_prime > obj.alpha_d) && (abs(curr_prime) > abs(max_curr)))
                            max_pow = power_prime;
                            max_curr = curr_prime;
                            curr_estimate = curr_prime;
                            b_estimate = b_prime;
                            v_estimate = v_prime;
                        end
                        
                    end
                    
                    msgID2 = 'Battery:maxDishargingPower';
                    % use v_estimate instead of v_prev because it is lower, to give some
                    % room for imprecise current estimates.
                    msg2 = num2str(max_pow);
                    causeException = MException(msgID2, msg2);
                    baseException = addCause(baseException,causeException);
                end    
                    
                if (obj.throw_exception)
                    throw(baseException);
                end
            
            % check that battery energy content does not exceed limits.
            % Throw an exception if it occurs.
            elseif curr_estimate > 0 && b_estimate > a_2_val + obj.tolerance
                msgID = 'Battery:Overcharge';
                msg = 'Battery energy content exceeds upper limit';
                baseException = MException(msgID,msg);
                
                % Search to find power corresponding to max current that
                % doesnt violate a_2 constraint.
                L = 0;
                R = u;
                power_resolution = (R-L)/50;
                power_range = L:power_resolution:R;
                max_pow = 0;
                max_curr = 0;
                curr_estimate = 0;
                b_estimate = obj.b_prev;
                v_estimate = obj.v_prev;

                for p_index = 1:numel(power_range)

                    power_prime = power_range(p_index);

                    [curr_prime, v_prime, b_prime, a_1_prime, a_2_prime, eff_c_prime, eff_d_prime] = iterateVoltage(obj,power_prime);

                    if ((b_prime < a_2_prime) && (abs(curr_prime) > abs(max_curr)))
                        max_pow = power_prime;
                        max_curr = curr_prime;
                        curr_estimate = curr_prime;
                        b_estimate = b_prime;
                        v_estimate = v_prime;
                    end

                end
                
                msgID2 = 'Battery:maxPowerBeforeOvercharge';
                msg2 = num2str(max_pow);
                causeException = MException(msgID2, msg2);
                baseException = addCause(baseException,causeException);
                
                if (obj.throw_exception)
                    throw(baseException);
                end
            
            elseif (curr_estimate < 0) && (b_estimate < (a_1_val - obj.tolerance))
                msgID = 'Battery:Undercharge';
                msg = 'Battery energy content below lower limit';
                baseException = MException(msgID,msg);
                
                % search to find power corresponding to maximum current
                % that doesn't break a_1 constraint.
                L = u;
                R = 0;
                power_resolution = (R-L)/50;
                power_range = L:power_resolution:R;
                max_pow = 0;
                max_curr = 0;
                curr_estimate = 0;
                b_estimate = obj.b_prev;
                v_estimate = obj.v_prev;
                for p_index = 1:numel(power_range)

                    power_prime = power_range(p_index);

                    [curr_prime, v_prime, b_prime, a_1_prime, a_2_prime, eff_c_prime, eff_d_prime] = iterateVoltage(obj,power_prime);

                    if ((b_prime > a_1_prime) && (abs(curr_prime) > abs(max_curr)))
                        max_pow = power_prime;
                        max_curr = curr_prime;
                        curr_estimate = curr_prime;
                        b_estimate = b_prime;
                        v_estimate = v_prime;
                    end

                end
                
                msgID2 = 'Battery:maxPowerBeforeUndercharge';
                msg2 = num2str(max_pow);
                causeException = MException(msgID2, msg2);
                baseException = addCause(baseException,causeException);
                if (obj.throw_exception)
                    throw(baseException);
                end
            end
            
            obj.V = v_estimate;
            obj.I = curr_estimate;
            obj.b = b_estimate;
            
            obj.v_prev = obj.V;
            obj.b_prev = obj.b;
            
            y = [obj.V, obj.I, obj.b];
            
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end

        function name1 = getInputNamesImpl(obj)
            % Return input port names for System block
            name1 = 'Power';
        end

        function name1 = getOutputNamesImpl(obj)
            % Return output port names for System block
            name1 = 'V, I, b';
        end

        function out = getOutputSizeImpl(obj)
            % Return size for each output port
            out = [1 3];
        end

        function out = getOutputDataTypeImpl(obj)
            % Return data type for each output port
            out = 'double';
        end

        function out = isOutputComplexImpl(obj)
            % Return true for each output port with complex data
            out = false;
        end

        function out = isOutputFixedSizeImpl(obj)
            % Return true for each output port with fixed size
            out = true;
        end

        function [sz,dt,cp] = getDiscreteStateSpecificationImpl(obj,name)
            % Return size, data type, and complexity of discrete-state
            % specified in name
            sz = [1 1];
            dt = 'double';
            cp = false;
        end
       
        function [curr_best, v_best, b_best, a_1_best, a_2_best, eff_c_best, eff_d_best] = iterateVoltage(obj,u)
            
            max_iterations = 100;
            
            v_best = 0;
            v_diff = inf;
            curr_best = 0;
            b_best = 0;
            a_1_best = 0;
            a_2_best = inf;
            eff_c_best = 1;
            eff_d_best = 1;      
            
            v_res = (obj.V_max - obj.V_min)/4;
            
            for v_0 = [obj.v_prev, obj.V_min:v_res:obj.V_max]
            
                v_estimate = v_0;
                curr_estimate = 0;
                b_estimate = obj.b_prev;

                % initialize change in energy in this time step
                delta_e = 0;

                % initialize a1,a2,eff_c and eff_d values
                a_1_val = 0;
                a_2_val = inf;
                eff_c_val = 1;
                eff_d_val = 1;      

                num_iterations = 0;

                while(num_iterations < max_iterations)

                    % estimate the current
                    curr_estimate = u/v_estimate;

                    % get a2 and eff_c values corresponding to the
                    % current, if it is positive (charging)

                    if (curr_estimate > 0)

                        a_2_val = interp1(obj.a_2(2,:), obj.a_2(1,:), curr_estimate, 'linear', 'extrap');
                        eff_c_val = interp1(obj.eff_c(2,:), obj.eff_c(1,:), curr_estimate, 'linear', 'extrap');

                        delta_e = u*eff_c_val*obj.time_step;

                    end

                    % get a1 and eff_d values corresponding to the
                    % current, if it is negative (discharging)

                   if (curr_estimate < 0)

                        a_1_val = interp1(obj.a_1(2,:), obj.a_1(1,:), curr_estimate, 'linear', 'extrap');
                        eff_d_val = interp1(obj.eff_d(2,:), obj.eff_d(1,:), curr_estimate, 'linear', 'extrap');

                        delta_e = (u/eff_d_val)*obj.time_step;

                   end

                    % estimate battery energy content
                    b_estimate = obj.b_prev + delta_e;

                    % estimate voltage
                    % Bound by V_min and V_max because we shouldnt be going
                    % out of this range.
                    v_new_estimate = max(obj.V_min, min(obj.V_max,(obj.M_function(curr_estimate, b_estimate))));
                    %v_new_estimate = obj.M_function(curr_estimate, b_estimate);

                    num_iterations = num_iterations + 1;
                    if (num_iterations == 100)
                        num_iterations = num_iterations;
                    end
                    

                    % stop iterating if the voltage value has converged
                    if (abs(v_new_estimate - v_estimate) < obj.tolerance)
                        v_estimate = v_new_estimate;
                        break;
                    else
                        v_estimate = v_new_estimate;
                    end

                end
                
                % check to see if estimate is closest to v_prev. If yes,
                % update the return values.
                if (abs(v_estimate - obj.v_prev) < (v_diff - obj.tolerance))
%                     if (v_diff < inf)
%                         disp(strcat('v_0 = ', num2str(v_0), ', v_prev = ', num2str(obj.v_prev), ', b_prev = ', num2str(obj.b_prev), ', b = ', num2str(b_estimate), ', u = ', num2str(u), ', v_best = ', num2str(v_estimate)));
%                     end
                    v_best = v_estimate;
                    curr_best = curr_estimate;
                    b_best = b_estimate;
                    a_1_best = a_1_val;
                    a_2_best = a_2_val;
                    eff_c_best = eff_c_val;
                    eff_d_best = eff_d_val;

                    v_diff = abs(v_estimate - obj.v_prev);
                end
                
            end
            
        end
        
    end
    
end
