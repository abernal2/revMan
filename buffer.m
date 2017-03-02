clear all, close all, clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differences from previous versions:
%           1) penalty cost per class
%           2) cancellation probability per class
%           3) Have to check whether customer arrived before at-capacity or
%              after, which determines what the system generated in revenue
%           4) Variable "booking_profile" contains an extra column which
%              determines whether that customer arrived or canceled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Revenue
r = [3;5;1;5;2;4;1;2]; 
% Penalty
pen = [4.5;7;2;6.5;4;5;2;3];    
% No cancellation probability
q = [0.7477; 0.3447; 0.6428; 0.6815; 0.7227; 0.7407; 0.9612; 0.6434]; 
% Poisson arrival rates
lam = [2; 3; 2; 2; 1; 2; 3; 1];
% Capacity
C = 60;
% Advanced reservation possibilities
advRes = [2;0;2;6;1;3;4;0];
% Possible service times
serTime = [1 4 2 1 1 4 4 3];
% Buffer size
xbuffer = 0;

%% Simulating the e-CSP policy to get the e-CSP long-run average revenue

%----------------------------Begin Simulation------------------------------

% Matlab's exponential distribution has pdf f(x)=(1/mu)*exp(-x/mu) with
% E[X]=mu. Also, since the results are derived with the assumption that the
% system is in steady state. Therefore, the simulation has to be started
% for a brief period of T time units until steady state is reached.

Tp = 10500;
T = 1000;
sim_ELRAR = 0;
const = Tp - T;
total_lam = sum(lam);
class_prob = [0; cumsum(lam)]/total_lam;
len = length(xbuffer);
sim_runs = 10;
tic
for i = 1:sim_runs
    
    % This is used for when sim_runs > 0
    revenue = 0;
    
    current_time = 0;
    % Initializing the arrival process. First arrival is always serviced.
    booking_profile = [];
    rnd_class = sum((class_prob < rand(1)));
    % This conveys whether customer canceled or not
    arrived = (rand(1) < q(rnd_class));

    % Checking whether the first arrival is eligible using the e-CSP policy.
    % The booking profile consists of [class, current time, start_time, end_time]
    % I got rid off the cutoff class because for analyzing the revenue, you
    % can reduce that problem to just looking at the accepted classes with
    % a modified rate, namely lam_i*(total_lam(accepted
    % classes)/total_lam(all classes)) for all accepted classes i. That
    % doesn't mean that the accepted classes will not be rejected; they can
    % be rejected if the capacity is full.
    
    current_time = exprnd(1/total_lam); 
    rnd_d = advRes(rnd_class); % class support is {0,1,...,d}
    rnd_s = serTime(rnd_class);
    booking_profile = [arrived, rnd_class, current_time, current_time + rnd_d, current_time + rnd_d + rnd_s];

    % The rest of the simulation up to time T. We assume the system has reached 
    % steady state at time T and will neglect that part of simulation.
    while(current_time <= T)

        rnd_class = sum((class_prob < rand(1)));
        current_time = current_time + exprnd(1/total_lam);

        % Determining if the specified interval is at capacity or not. We dont need
        % the finish time of the requested service interval because it only matters
        % how many customers are/or will be using it at the beginning of the requested      
        % service interval.

        % Determining advance reservation and service time
        rnd_d = advRes(rnd_class); 
        rnd_s = serTime(rnd_class);

        % Time the current customer is requesting to start
        temp = current_time + rnd_d; 

        % Deleting customers in the booking profile that will not
        % factor into the maximum usage.
        if(~isempty(booking_profile))
            bool_array = (booking_profile(:,end) <= temp | booking_profile(:,end-1) >= temp + rnd_s);  
            dummy_profile = booking_profile(~bool_array,:);   
            
            % The reason I didn't have to change the "check_Max_Usage()"
            % function is because I am supplying a modified "dummy_profile"
            % variable (see next 2 lines). 
            bool_array = (dummy_profile(:,end-1) <= current_time & dummy_profile(:,1) == 0);
            dummy_profile(bool_array,:) = [];
            if(~isempty(dummy_profile))
                max_usage = check_Max_Usage(dummy_profile(:,end-1),dummy_profile(:,end));
            else
                max_usage = 0;
            end
        else
            max_usage = 0;
        end

        if (max_usage < C + xbuffer(i)) % Books when booking profile is NOT at capacity

            % Add current customer
            arrived = (rand(1) < q(rnd_class));
            booking_profile(end+1,:) = [arrived, rnd_class, current_time, temp, temp + rnd_s]; 

        end

    end

    % The rest of the simulation up to time Tp. We assume the system has reached 
    % steady state at time T and will neglect that part of simulation. Now, we
    % record the revenue after time T up to time Tp. These records will be used
    % to calculate the long-run average revenue over the time interval [T,Tp].
    % First, we need to calculate the revenue of the customers ALREADY in the
    % system at time T. Then, the simulation is resumed up to time Tp.

    bool_array = (booking_profile(:,end) <= T);
    booking_profile(bool_array,:) = [];

    % Checking the case when all customers are done before T, in which case
    % the line above will return: "Empty matrix: 4-by-0"
    if(isempty(booking_profile)) booking_profile = []; end
    rej = 0;
    
    while(current_time <= Tp)

        rnd_class = sum((class_prob < rand(1)));
        current_time = current_time + exprnd(1/total_lam);

        % Determining the advance reservation and service time
        rnd_d = advRes(rnd_class); 
        rnd_s = serTime(rnd_class);

        % Determining if the specified interval is at capacity or not. We dont need
        % the finish time of the requested service interval because it only matters
        % how many customers are/or will be using it at the beginning of the requested      
        % service interval.
        temp = current_time + rnd_d; % Time the current customer is requesting to start

        % Deleting customers in the booking profile that will not
        % factor into the maximum usage.
        if(~isempty(booking_profile))
            bool_array = (booking_profile(:,end) <= temp | booking_profile(:,end-1) >= temp + rnd_s);  
            dummy_profile = booking_profile(~bool_array,:);   
            
            % The reason I didn't have to change the "check_Max_Usage()"
            % function is because I am supplying a modified "dummy_profile"
            % variable (see next 2 lines). 
            bool_array = (dummy_profile(:,end-1) <= current_time & dummy_profile(:,1) == 0);
            dummy_profile(bool_array,:) = [];
            if(~isempty(dummy_profile))
                max_usage = check_Max_Usage(dummy_profile(:,end-1),dummy_profile(:,end));
            else
                max_usage = 0;
            end
        else
            max_usage = 0;
        end

        if (max_usage < C + xbuffer(i)) % Books when booking profile is NOT at capacity

            % Add current customer
            arrived = (rand(1) < q(rnd_class));
            booking_profile(end+1,:) = [arrived, rnd_class, current_time, temp, temp + rnd_s]; 

        else %<---- I put this 'else' statement to test if customers are being rejected
            rej = rej + 1;
        end

    end

    % Deleting customers who will START service after time Tp
     bool_array = (booking_profile(:,end-1) >= Tp);
     booking_profile(bool_array,:) = [];

    % Calculating the revenue for customers who received service in the
    % interval [T,Tp] for one simulation
    % revenue = sum(arrayfun(price,booking_profile(:,1)))/const; (OLD)   
    bool_array = (booking_profile(:,1) == 0);
    booking_profile(bool_array,:) = [];
    leng = size(booking_profile, 1);
    temp = [booking_profile(:,end-1) ones(leng,1); booking_profile(:,end) -ones(leng,1)];
    temp = sortrows(temp);
    temp = [temp(:,1) cumsum(temp(:,2))];
    
    for k = 1:leng
        idx = max(find(temp(:,1) <= booking_profile(i,end-1)-0.0001));
        if (temp(idx,2) < C)
            revenue = revenue + r(booking_profile(i,2));
        else
            revenue = revenue + r(booking_profile(i,2)) - p(booking_profile(i,2));
        end
    end

    % Updating the long-run average revenue. This computes a running
    % average so there is no need to store each simulated revenue
    sim_ELRAR = ((i-1)*sim_ELRAR + revenue)/i;

end


