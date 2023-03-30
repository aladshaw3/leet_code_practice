classdef obj1
    %obj1 -

    properties
        %   Data objects
        x
        varMap
        t
    end

    methods
        % Default Constructor
        function obj = obj1()
            %obj1 Construct an instance of this class
            %   Set default values for all properties
            obj.x = zeros(1,1);
            obj.varMap = containers.Map({'x1'},{obj.x(1)})
            obj.t = 0
        end

        function F = Residual(obj, x)
            %Residual -
            F(1) = 0;
        end

        function Edit(obj)
            %Residual -
            disp('here')
            obj.x(1) = 0;
            obj.t = 1
        end

    end

end
