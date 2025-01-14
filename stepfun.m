function stepy = stepfun(t,nsteps,maxp,delay, T)
%Creates stair-step function used to measure 16 s pressure-volume loop
%   Uses tanh functions for smoothness

difp = maxp/nsteps;
tw = 0.05; % make narrow tanh
step_up = @(t, t0, tw) 0.5*difp + 0.5*difp*tanh((t - t0)/tw);
step_down = @(t, t0, tw) -0.5*difp - 0.5*difp*tanh((t - t0)/tw);

t0 = (0:(T/2/nsteps):t) + delay;
sms = zeros(length(t), length(t0));
for i = 1:length(t0)
    if t0(i) <= T/2
       sms(:,i) = step_up(t, t0(i), tw);
    elseif t0(i) > T/2
       sms(:,i) = step_down(t, t0(i), tw);
    end
end

stepy = sum(sms, 2);

end
