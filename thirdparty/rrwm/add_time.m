function [time] = add_time(time1, time2)

time = time1;
time(6) = time(6) + time2;

%% Seconds
while time(6) >= 60
    time(5) = time(5) + 1;
    time(6) = time(6) - 60;
end

%% Minutes
while time(5) >= 60
    time(4) = time(4) + 1;
    time(5) = time(5) - 60;
end

%% Hours
while time(4) >= 24
    time(3) = time(3) + 1;
    time(4) = time(4) - 24;
end
