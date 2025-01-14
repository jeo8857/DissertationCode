function [position,isterminal,direction] = MultipleBreathsEvent(t,s)
s
pause
% Tind = find(abs((breaths-1)*res.T-t)==min(abs((breaths-1)*res.T-t)));
  position = s(12)-2.5*10^4; % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end