function result = simps(t,f)
%     for i = ind1:ind2
%         dt = t(i)-t(i-1);
%         rec(i) = dt*func(i);
%     end
%     area = sum(rec);
%     avg = sum(rec)/(ind2-ind1);

% Simpsons rule for irregularly spaced data
N = length(t)-1;
for i = 1:N-1
    h(i) = t(i+1)-t(i);
end

result = 0;
for i = 2:2:N-1
    h0 = h(i-1); h1 = h(i);
    hph = h1+h0; hdh = h1/h0; hmh = h1*h0;
    result = result+(hph/6)*( ...
        (2-hdh)*f(i-1)+(hph^2/hmh)*f(i)+(2-1/hdh)*f(i+1) ...
        );
end
    
end