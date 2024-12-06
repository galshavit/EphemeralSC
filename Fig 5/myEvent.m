function [value, isterminal, direction] = myEvent(~, y)
value      = y<1e-3;
isterminal = 1;   % Stop the integration
direction  = 0;
end