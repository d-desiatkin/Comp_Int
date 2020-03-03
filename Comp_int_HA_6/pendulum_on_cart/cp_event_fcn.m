function [position,isterminal,direction] = cp_event_fcn(t,x)
    position = cos(x(2)); % The value that we want to be zero
    isterminal = 1;  % Halt integration 
    direction = 0; 
end