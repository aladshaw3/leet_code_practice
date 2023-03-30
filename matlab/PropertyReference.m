classdef PropertyReference < handle
    %PropertyReference Reference to a property in another object
    properties
        sourceHandle
        sourceFieldName
    end

    properties (Dependent = true)
         Value
    end

    methods
        function obj = PropertyReference (source, fieldName)
            obj.sourceHandle = source;
            obj.sourceFieldName = fieldName
        end
        function value = get.Value( obj )
            value = obj.sourceHandle.(obj.sourceFieldName);
        end

        function set.Value( obj, value )
            obj.sourceHandle.(obj.sourceFieldName) = value;
        end
        function disp( obj )
            disp(obj.Value);
        end
    end
end
