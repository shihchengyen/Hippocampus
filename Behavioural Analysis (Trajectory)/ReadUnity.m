function ReadUnity(data)

% col 3 & 4 are x & y
% col 5 is angle in degrees
% but to display correctly, we need to multiply by -1 and add 90 degrees
% then we will convert to radians as the sin and cos functions expect radians
dangle = deg2rad(-data(:,5) + 90);    

for i = 1:size(data,1);  

% plot using quiver function
% quiver(data(:,3),data(:,4),cos(dangle),sin(dangle));      
pause(0.02); %plot data after delta time 
    
%     if data(i,1) ~= 0,
%         disp(data(i,1));
%     else
%     end

    q = quiver(data(i,3),data(i,4),cos(dangle(i)),sin(dangle(i)));
    q.LineWidth = 1.5;
    hold on; %retains plots in current axes so that new plots added to the axes do not delete existing plots
    xlim([-15 15]);
    ylim([-10 10]);
end
