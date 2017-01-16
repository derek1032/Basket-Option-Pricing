classdef BarrierOption
   
    properties
        S10;
        S20;
        Sb; % barrier
        K;
        r;
        T;
        sigma1;
        sigma2;
        rho;
        method;
        NStep;
        NRepl;    
    end
    
    methods
        function ob = BarrierOption(S10,S20,Sb,K,r,T,sigma1,sigma2,rho,method,NStep,NRepl)
            ob.S10 = S10;
            ob.S20 = S20;
            ob.Sb = Sb;
            ob.K = K;
            ob.r = r;
            ob.T = T;
            ob.sigma1 = sigma1;
            ob.sigma2 = sigma2;
            ob.rho = rho;
            ob.method = method;
            ob.NStep = NStep;
            ob.NRepl = NRepl;
        end
       
        %%%%%%%%%%%%%%%%%%%%%%% Finite Difference %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% Explicit Method   %%%%%%%%%%%%%%%%%%%%%%%%
        function p = FiniteDiff(ob,M,N)
            %Define the grid
            s1_max = 500;
            s2_max = 500;
            dt = ob.T/M;
            dx = s1_max/(N);
            dy = s2_max/(N);
            sol_prev = zeros(N+1,N+1);
            %Initialize the sol_post
            grid = repmat((0:1:N)*dx,N+1,1)+repmat((0:1:N)'*dy,1,N+1);
            sol_post = max(grid-ob.K,0);
            for k = 1:M
                sol_post = sol_post.*(grid<ob.Sb);
                for i = 2 : N
                    for j = 2: N
                        
                        sol_prev(i,j) = (1-dt*ob.sigma1^2*(i-1)^2-dt*ob.sigma2^2*(j-1)^2-dt*ob.r)*sol_post(i,j)...
                                        +(0.5*dt*ob.sigma1^2*(i-1)^2+0.5*dt*(i-1)*ob.r)*sol_post(i+1,j)...
                                        +(0.5*dt*ob.sigma1^2*(i-1)^2-0.5*dt*(i-1)*ob.r)*sol_post(i-1,j)...
                                        +(0.5*dt*ob.sigma2^2*(j-1)^2+0.5*dt*(j-1)*ob.r)*sol_post(i,j+1)...
                                        +(0.5*dt*ob.sigma2^2*(j-1)^2-0.5*dt*(j-1)*ob.r)*sol_post(i,j-1)...
                                        + 0.25*ob.rho*ob.sigma1*ob.sigma2*(i-1)*(j-1)*dt*sol_post(i+1,j+1)...
                                        - 0.25*ob.rho*ob.sigma1*ob.sigma2*(i-1)*(j-1)*dt*sol_post(i+1,j-1)...
                                        - 0.25*ob.rho*ob.sigma1*ob.sigma2*(i-1)*(j-1)*dt*sol_post(i-1,j+1)...
                                        + 0.25*ob.rho*ob.sigma1*ob.sigma2*(i-1)*(j-1)*dt*sol_post(i-1,j-1);
                                    
                    end 
                end
                sol_prev(1,:) = 2*sol_prev(2,:)-sol_prev(3,:);
                sol_prev(N+1,:) = 2*sol_prev(N,:)-sol_prev(N-1,:);
                sol_prev(:,1) = 2*sol_prev(:,2)-sol_prev(:,3);
                sol_prev(:,N+1) = 2*sol_prev(:,N)-sol_prev(:,N-1);
                sol_post = sol_prev;
                sol_prev = zeros(N+1,N+1);
            end
            p = sol_post;
            %p = interp2(0:dx:s1_max,0:dx:s2_max, sol_post, ob.S10, ob.S20);
         
        end
        %%%%%%%%%%%%%%%%implicit method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function p = Finite_diff_implicit(ob,M,N)
            %Define the grid
            s1_max = 500;
            s2_max = 500;
            dt = ob.T/ob.NStep;
            dx = s1_max/(N);
            dy = s2_max/(N);
            I1 = ones(N,N);
            I1(1,1:2) = 0;
            I1(N,(N-1):N) = 0;
            
            I2 = diag([1,zeros(1,(N-1))]);
            I3 = diag([1,zeros(1,N-2)],1);
            I4 = diag([zeros(1,N-1),1]);
            I5 = diag([zeros(1,N-2),1],-1);
            
            a1 = @(i,j) (1+dt*ob.sigma1^2*(i).^2+dt*ob.sigma2^2*(j).^2+dt*ob.r);
            a2 = @(i) (-0.5*dt*ob.sigma1^2*(i).^2-0.5*dt*(i).*ob.r);
            a3 = @(i) (-0.5*dt*ob.sigma1^2*(i).^2+0.5*dt*(i).*ob.r);
            a4 = @(j) (-0.5*dt*ob.sigma2^2*(j).^2-0.5*dt*(j).*ob.r);
            a5 = @(j) (-0.5*dt*ob.sigma2^2*(j).^2+0.5*dt*(j).*ob.r);
            a6 = @(i,j) (-0.25*ob.rho*ob.sigma1*ob.sigma2*(i).*(j).*dt);
            a7 = @(i,j) (0.25*ob.rho*ob.sigma1*ob.sigma2*(i).*(j).*dt);
            a8 = @(i,j) (0.25*ob.rho*ob.sigma1*ob.sigma2*(i).*(j).*dt);
            a9 = @(i,j) (-0.25*ob.rho*ob.sigma1*ob.sigma2*(i).*(j).*dt);
            %left j-1
            
            A = @(j) (a5(j)*diag(ones(N,1))+diag(a7(1:N-1,j),1)+diag(a9(2:N,j),-1))...
                .*I1 + (2*a9(1,j)+a5(j)).*I2+(a7(1,j)-a9(1,j)).*I3+(a5(j)+2*a9(N,j)).*I4+(a9(N,j)-a7(N,j)).*I5;
            %center j
            B = @(j) (diag(a1(1:N,j))+diag(a2(1:N-1),1)+diag(a3(2:N),-1))...
                .*I1 + (2*a3(1)+a1(1,j)).*I2+(a2(1)-a3(1)).*I3+(a1(N,j)+2*a3(N)).*I4+(a3(N)-a2(N)).*I5;
            %right j+1
            C = @(j) (a4(j)*diag(ones(N,1))+diag(a6(1:N-1,j),1)+diag(a8(2:N,j),-1))...
                .*I1 + (2*a8(1,j)+a4(j)).*I2+(a6(1,j)-a8(1,j)).*I3+(a4(j)+2*a8(N,j)).*I4+(a8(N,j)-a6(N,j)).*I5;
            
            
            
            grid = repmat((1:1:N)*dx,N,1)+repmat((1:1:N)'*dy,1,N);
            sol_post = max(grid-ob.K,0);
            I = (grid)<ob.Sb;
            sol_post_reshape = reshape(sol_post,N*N,1);
            I_reshape = reshape(I,N*N,1);
            
            Z = cell(N,N);
            %initilize
            for i = 1:N
                for j = 1:N
                    Z{i,j} = zeros(N,N);
                end
            end 
            Z{1,1} = 2*A(1)+B(1);
            Z{1,2} = C(1)-A(1);
            Z{N,N} = B(N)+2*C(N);
            Z{N,N-1} = A(N)-C(N); 
            
            for j = 2:N-1
 
                   Z{j,j-1} = A(j);
                   Z{j,j} = B(j);
                   Z{j,j+1} = C(j);     
            end 
            
            U = cell2mat(Z);
            D = U^(-1);
             for k = 1:ob.NStep
                sol_post_reshape = sol_post_reshape.*I_reshape;
                sol_prev_reshape = D*sol_post_reshape;
                sol_post_reshape = sol_prev_reshape; 
                k
             end 
                
            sol = reshape(sol_post_reshape,N,N);
            p = sol;
            %p = interp2(dx:dx:s1_max,dy:dy:s2_max, sol,ob.S10,ob.S20);
            
        end 
        
        function p = Finite_diff_price(ob,S1,S2)
            s1_max = 500;
            s2_max = 500;
            dx = s1_max/100;
            dy = s2_max/100;
            sol = ob.Finite_diff_implicit(1000,100);
            p = interp2(dx:dx:s1_max,dy:dy:s2_max, sol, S1,S2,'spline');
        end 
           
        function [delta_x,delta_y] = Finite_diff_delta(ob,S1,S2)
            s1_max = 500;
            s2_max = 500;
            dx = s1_max/100;
            dy = s2_max/100;
            delta_x = 0.1;
            delta_y = 0.1;
            %sol = ob.FiniteDiff(1000,100);
            sol = ob.Finite_diff_implicit(100,100);
            p1 = interp2(dx:dx:s1_max,dy:dy:s2_max, sol, S1,S2,'spline');
            p2 = interp2(dx:dx:s1_max,dy:dy:s2_max, sol, S1+delta_x,S2,'spline');
            p3 = interp2(dx:dx:s1_max,dy:dy:s2_max, sol, S1,S2+delta_y,'spline');
            delta_x = (p2-p1)/delta_x;
            delta_y = (p3-p1)/delta_y;
            
        end
        
        function [g1,g2] = Finite_diff_gamma(ob,S1,S2)
            s1_max = 500;
            s2_max = 500;
            dx = s1_max/5000;
            dy = s2_max/5000;
            d1 = ob.Finite_diff_delta(S1,S2);
            d2 = ob.Finite_diff_delta(S1+dx,S2);
            d3 = ob.Finite_diff_delta(S1,S2+dy);
            g1 = (d2-d1)/dx;
            g2 = (d3-d1)/dy;
            
        end 
        
        
        function [v1,v2] = Finite_diff_vega(ob)
            p1 = ob.Finite_diff_price(ob.S10,ob.S20);
            delta_sigma1 = 0.01;
            delta_sigma2 = 0.01;
            ob_temp1 = BarrierOption(ob.S10,ob.S20,ob.Sb,ob.K,...
            ob.r,ob.T,ob.sigma1+delta_sigma1,ob.sigma2,ob.rho,ob.method,ob.NStep,ob.NRepl);
            ob_temp2 = BarrierOption(ob.S10,ob.S20,ob.Sb,ob.K,...
            ob.r,ob.T,ob.sigma1,ob.sigma2+delta_sigma2,ob.rho,ob.method,ob.NStep,ob.NRepl);
            p2 = ob_temp1.Finite_diff_price(ob.S10,ob.S20);
            p3 = ob_temp2.Finite_diff_price(ob.S10,ob.S20);
            v1 = (p2-p1)/delta_sigma1;
            v2 = (p3-p1)/delta_sigma2;
      
        end 
        
         function t = Finite_diff_theta(ob)
            p1 = ob.Finite_diff_price(ob.S10,ob.S20);
            
            delta_t = 0.1;
            ob_temp1 = BarrierOption(ob.S10,ob.S20,ob.Sb,ob.K,...
            ob.r,ob.T+delta_t,ob.sigma1,ob.sigma2,ob.rho,ob.method,ob.NStep,ob.NRepl);
            p2 = ob_temp1.Finite_diff_price(ob.S10,ob.S20);
            t = -(p2-p1)/delta_t;
            
         end
        
         
         %%%%%%%%%%%%%%%%%%%%%%%%% Monte Carlo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [price,error] = MonteCarlo(ob)
            rng(1);
            dt = ob.T/ob.NStep;
            % pilot simulation
            NPilot = 4000;
            R  = randn(NPilot,ob.NStep);
            multi1 = exp((ob.r-0.5*ob.sigma1^2)*dt+ob.sigma1*sqrt(dt)*R);
            multi2 = exp((ob.r-0.5*ob.sigma2^2)*dt+ob.sigma2*sqrt(dt)*(ob.rho*R+sqrt(1-ob.rho^2)*randn(NPilot,ob.NStep)));
            S1Path = cumprod([ob.S10*ones(NPilot,1),multi1],2);
            S2Path = cumprod([ob.S20*ones(NPilot,1),multi2],2);
            crossIdx = abs(any(S1Path+S2Path>=ob.Sb,2)-1);
            callprice = exp(-ob.r*ob.T)*max(0,S1Path(:,ob.NStep+1)+S2Path(:,ob.NStep+1)-ob.K).*crossIdx;
            covmat = cov(callprice,S1Path(:,ob.NStep+1));
            c = -covmat(1,2)/(ob.S10^2*exp(2*ob.r*ob.T)*(exp(ob.sigma1^2*ob.T)-1));
            % simulation
            R  = randn(ob.NRepl,ob.NStep);
            multi1 = exp((ob.r-0.5*ob.sigma1^2)*dt+ob.sigma1*sqrt(dt)*R);
            multi2 = exp((ob.r-0.5*ob.sigma2^2)*dt+ob.sigma2*sqrt(dt)*(ob.rho*R+sqrt(1-ob.rho^2)*randn(ob.NRepl,ob.NStep)));
            S1Path = cumprod([ob.S10*ones(ob.NRepl,1),multi1],2);
            S2Path = cumprod([ob.S20*ones(ob.NRepl,1),multi2],2);
            crossIdx = abs(any(S1Path+S2Path>=ob.Sb,2)-1);
            callprice = exp(-ob.r*ob.T)*max(0,S1Path(:,ob.NStep+1)+S2Path(:,ob.NStep+1)-ob.K).*crossIdx;
            price = mean(callprice + c*(S1Path(:,ob.NStep+1)-ob.S10*exp(ob.r*ob.T)));
            error = std(callprice + c*(S1Path(:,ob.NStep+1)-ob.S10*exp(ob.r*ob.T)))*norminv(0.975)/sqrt(ob.NRepl);

        end
         %%%%%%%%%%%%%%%%%%%%%%%% Tree Method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [price,error] = Tree(ob)
            % invariant quantities
            N = ob.NStep;
            deltaT = ob.T/N;
            nu1 = ob.r - 0.5*ob.sigma1^2;
            nu2 = ob.r - 0.5*ob.sigma2^2;
            u1 = exp(ob.sigma1*sqrt(deltaT));
            d1 = 1/u1;
            u2 = exp(ob.sigma2*sqrt(deltaT));
            d2 = 1/u2;
            discount = exp(-ob.r*deltaT);
            p_uu = discount*0.25*(1+sqrt(deltaT)*(nu1/ob.sigma1+nu2/ob.sigma2)+ob.rho);
            p_ud = discount*0.25*(1+sqrt(deltaT)*(nu1/ob.sigma1-nu2/ob.sigma2)-ob.rho);
            p_du = discount*0.25*(1+sqrt(deltaT)*(-nu1/ob.sigma1+nu2/ob.sigma2)-ob.rho);
            p_dd = discount*0.25*(1+sqrt(deltaT)*(-nu1/ob.sigma1-nu2/ob.sigma2)+ob.rho);
            
            % set up stock values
            S1Vals = zeros(2*N+1,1);
            S2Vals = zeros(2*N+1,1);
            S1Vals(1) = ob.S10*d1^N;
            S2Vals(1) = ob.S20*d2^N;
            for i = 2:2*N+1
               S1Vals(i) = u1*S1Vals(i-1);
               S2Vals(i) = u2*S2Vals(i-1);
            end
            
            % set up terminal values
            CVals = zeros(2*N+1,2*N+1);
            for i=1:2:2*N+1
                for j=1:2:2*N+1
                    if (S1Vals(i) + S2Vals(j))>= ob.Sb;
                        CVals(i,j) = 0;
                    else
                        CVals(i,j) = max(0,S1Vals(i)+S2Vals(j)-ob.K);
                    end
                end
            end
            % roll back
            for tau=1:N
                for i = (1+tau):2:(2*N+1-tau)
                    for j = (1+tau):2:(2*N+1-tau)
                        if (S1Vals(i) + S2Vals(j))>= ob.Sb
                            CVals(i,j) = 0;
                        else
                            CVals(i,j) = p_uu*CVals(i+1,j+1)+ p_ud*CVals(i+1,j-1)...
                                +p_du*CVals(i-1,j+1)+p_dd*CVals(i-1,j-1);
                        end
                    end
                end
            end
            price = CVals(N+1,N+1);
            error = NaN;
        end
        
        function [price,error] = computePrice(ob)
            switch ob.method
                case 1
                    [price,error] = ob.MonteCarlo();
                case 2
                    price = ob.Finite_diff_price(ob.S10,ob.S20);
                    error = NaN;
                case 3
                    [price,error] = ob.Tree();
            end
        end
        %%%%%%%%%%%%%%%%COMPUTE THETA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function theta = computeTheta(ob,incr)
             switch ob.method
                case 1
                    oldT = ob.T; % store original T
                    newT = ob.T * (1+incr);
                    deltaT = newT - oldT;
                    price1 = ob.MonteCarlo();
                    ob.T = newT;
                    price2 = ob.MonteCarlo();
                    theta = -(price2 - price1)/deltaT; % compute vega2
                    ob.T = oldT;
                case 2
                   theta = ob.Finite_diff_theta();
                case 3
                    oldT = ob.T; % store original T
                    newT = ob.T * (1+incr);
                    deltaT = newT - oldT;
                    price1 = ob.Tree();
                    ob.T = newT;
                    price2 = ob.Tree();
                    theta = -(price2 - price1)/deltaT; % compute vega2
                    ob.T = oldT;
            end
            
        end
        %%%%%%%%%%%%%%%%COMPUTE VEGA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [vega1,vega2] = computeVega(ob,incr)
             switch ob.method
                case 1
                    oldS2 = ob.sigma2; % store original sigma2
                    newS2 = ob.sigma2 * (1+incr);
                    deltaS = newS2 - oldS2;
                    price1 = ob.MonteCarlo();
                    ob.sigma2 = newS2;
                    price2 = ob.MonteCarlo();
                    vega2 = (price2 - price1)/deltaS; % compute vega2
                    ob.sigma2 = oldS2;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    oldS1 = ob.sigma1; % store original sigma1
                    newS1 = ob.sigma1 * (1+incr);
                    deltaS = newS1 - oldS1;
                    price1 = ob.MonteCarlo();
                    ob.sigma1 = newS1;
                    price2 = ob.MonteCarlo();
                    vega1 = (price2 - price1)/deltaS; % compute vega1
                    ob.sigma1 = oldS1;
                    
                    
                case 2
                    [vega1,vega2] = ob.Finite_diff_vega();
                    
                case 3
                    oldS2 = ob.sigma2; % store original sigma2
                    newS2 = ob.sigma2 * (1+incr);
                    deltaS = newS2 - oldS2;
                    price1 = ob.Tree();
                    ob.sigma2 = newS2;
                    price2 = ob.Tree();
                    vega2 = (price2 - price1)/deltaS; % compute vega2
                    ob.sigma2 = oldS2;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    oldS1 = ob.sigma1; % store original sigma1
                    newS1 = ob.sigma1 * (1+incr);
                    deltaS = newS1 - oldS1;
                    price1 = ob.Tree();
                    ob.sigma1 = newS1;
                    price2 = ob.Tree();
                    vega1 = (price2 - price1)/deltaS; % compute vega1
                    ob.sigma1 = oldS1;
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%COMPUTE GAMMA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [gamma1,gamma2] = computeGamma(ob,incr)
            switch ob.method
                case 1
                    [delta1_1,delta1_2] = ob.computeDelta(incr);
                    [delta2_1,delta2_2] = ob.computeDelta(-incr);
                    deltaS_1 = (ob.S10*(2+incr))/2 - (ob.S10*(2-incr))/2;
                    deltaS_2 = (ob.S20*(2+incr))/2 - (ob.S20*(2-incr))/2;
                    gamma1 = (delta1_1 - delta2_1)/deltaS_1;
                    gamma2 = (delta1_2 - delta2_2)/deltaS_2;
                    
                case 2
                    
                    [gamma1,gamma2] = ob.Finite_diff_gamma(ob.S10,ob.S20);
                    
                case 3
                    [delta1_1,delta1_2] = ob.computeDelta(incr);
                    [delta2_1,delta2_2] = ob.computeDelta(-incr);
                    deltaS_1 = (ob.S10*(2+incr))/2 - (ob.S10*(2-incr))/2;
                    deltaS_2 = (ob.S20*(2+incr))/2 - (ob.S20*(2-incr))/2;
                    gamma1 = (delta1_1 - delta2_1)/deltaS_1;
                    gamma2 = (delta1_2 - delta2_2)/deltaS_2;
            end
            
            
        end
        %%%%%%%%%%%%%%%%%%COMPUTE DELTA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [delta1,delta2] = computeDelta(ob,incr)
            switch ob.method
                case 1
                    oldS1 = ob.S10; % store original S10
                    newS1 = ob.S10 * (1+incr);
                    deltaS = newS1 - oldS1;
                    price1 = ob.MonteCarlo();
                    ob.S10 = newS1;
                    price2 = ob.MonteCarlo();
                    delta1 = (price2 - price1)/deltaS; % compute delta1
                    ob.S10 = oldS1;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    oldS2 = ob.S20; % store original S20
                    newS2 = ob.S20 * (1+incr);
                    deltaS = newS2 - oldS2;
                    price1 = ob.MonteCarlo();
                    ob.S20 = newS2;
                    price2 = ob.MonteCarlo();
                    delta2 = (price2 - price1)/deltaS; % compute delta2
                    ob.S20 = oldS2;
                case 2
                    
                    [delta1,delta2] = ob.Finite_diff_delta(ob.S10,ob.S20);
                    
                case 3
                    oldS1 = ob.S10; % store original S10
                    newS1 = ob.S10 * (1+incr);
                    deltaS = newS1 - oldS1;
                    price1 = ob.Tree();
                    ob.S10 = newS1;
                    price2 = ob.Tree();
                    delta1 = (price2 - price1)/deltaS; % compute delta1
                    ob.S10 = oldS1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    oldS2 = ob.S20; % store original S20
                    newS2 = ob.S20 * (1+incr);
                    deltaS = newS2 - oldS2;
                    price1 = ob.Tree();
                    ob.S20 = newS2;
                    price2 = ob.Tree();
                    delta2 = (price2 - price1)/deltaS; % compute delta2
                    ob.S20 = oldS2;
            end
            
            
         end 
         
         

         
    end 
        
        
end
