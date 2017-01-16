%% load parameters
file = fopen('test6.csv');
data = textscan(file,'%f%s','delimiter',',');
fclose(file);
%%
p = data{1};
barrierCall = BarrierOption(p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),...
                             p(10),p(11),p(12));
tic,[price,error] = barrierCall.computePrice();toc
%% 
tic,[delta1,delta2] = barrierCall.computeDelta(0.01);toc

%% 
tic,[gamma1,gamma2] = barrierCall.computeGamma(0.01);toc

%%
tic,[vega1,vega2] = barrierCall.computeVega(0.01);toc
%% 
tic,theta = barrierCall.computeTheta(0.01);toc


%% ouput csv file
c = {price 'option price';error '95% error';delta1 'delta1';delta2 'delta2';...
    gamma1 'gamma1';gamma2 'gamma2';vega1 'vega1';vega2 'vega2';theta 'theta'};
f = fopen('output.csv', 'w') ;
for i = 1:length(c)
fprintf(f, '%f,%s\n', c{i,:});
end
fclose(f) ;



