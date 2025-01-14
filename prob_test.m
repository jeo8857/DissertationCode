r = [0.0058 0.0042 .0031 .0024]';
L = [0.0436 .0178 .0096 .01]';
Ba = [.5760 .5934 .3840 .3491]';
par = F_Parameters();
TV = 0.5/1e3;
par.dp = 10*10^(-6);
p = G_impaction(r,L,Ba,par,TV);

% [p1_in,p1_out] = G_impaction(r(1),L(1),Ba(1),par,TV)
% [p2,~] = G_impaction(r(2),L(2),Ba(2),par,TV)
% p3 = G_impaction(r(3),L(3),Ba(3),par,TV)
% p4 = G_impaction(r(4),L(4),Ba(4),par,TV)
% pb = p(1)+p(2)+p(3)+p(4)...
%     -(p(1)*p(2)+p(1)*p(3)+p(1)*p(4)+p(2)*p(3)+p(2)*p(4)+p(3)*p(4))...
%     +(p(1)*p(2)*p(3)+p(1)*p(2)*p(4)+p(1)*p(3)*p(4)+p(2)*p(3)*p(4))...
%     -p(1)*p(2)*p(3)*p(4)

% p(1) = 0.047637683104654; %0.0159;
% p(2) = 0.0227;
% p(3) = 0.0200;
% p(4) = 0.0231;
% 
% p = [p(1) p(2) p(3) p(4)];
% 
% % k = length(p);
% % 
% % % when k = 1
% % A = sum(p);
% % 
% % for i = 1:length(p)
% %     if i==1
% %         p(i)*p(2:end);
% %     elseif i==length(p)
% %         p(i)*p(1:end-1);
% %     else
% %         p(i)*[p(1:i-1) p(i+1:end)];
% %     end
% % end
% 

