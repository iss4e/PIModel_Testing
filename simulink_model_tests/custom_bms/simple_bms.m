classdef simple_bms < matlab.System
    % untitled Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties

    end

    properties(Nontunable)
        Vmin = 1.5; % minimum voltage
        Vmax = 2.7;  % maximum voltage
    end
        
    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)
        
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
        end

        function y = stepImpl(obj,voltage,current)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            y = 1; % default case, dont do anything

            % if we are discharging and the voltage falls below the
            % minimum, open the switch
            if (current < 0 && voltage < obj.Vmin)
                y = 0;
            end
            % if we are charging and the voltage rises above the maximum,
            % open the switch
            if (current > 0 && voltage > obj.Vmax)
                y = 0;
            end

        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end

        function num = getNumInputsImpl(obj)
            % Define total number of inputs for system with optional inputs
            num = 2;
            % if obj.UseOptionalInput
            %     num = 2;
            % end
        end
    end
end
