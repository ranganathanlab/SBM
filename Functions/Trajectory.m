function [J_traj,h_traj] = Trajectory(wt,q,verbose,options)
% record the trajectory of an inference through the mean of its parameters
% for the purpose of checking for convergence
    for i=1:(size(wt,1)-1)
        [J,h]=Jw(wt(i,:)',q);
        [J,h]=IsingGauge(J,h);
        NORMS=Frob(J);
        J_traj(i)=mean(mean(NORMS));
        h_traj(i)=mean(mean(h.^2));
    end

    if verbose==1
        figure()
        hold on
        plot(J_traj);
        title(' mean J trajectory')
        ylabel('mean J')
        xlabel(['iterations (in steps of ' num2str(options.skip) ')'])
        figure()
        hold on
        plot(h_traj);
        title(' mean h trajectory')
        ylabel('mean h')
        xlabel(['iterations (in steps of ' num2str(options.skip) ')'])
    end
end

