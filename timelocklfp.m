function [timelock ] = timelocklfp(cfg,data,figure_dir,figurename)
%written by Seth Konig 11/8/2015
%calculate time/event lock avearge LFP and plot and save it.

timelock = ft_timelockanalysis(cfg,data);
figure
plot(timelock.time,timelock.avg')
xlim([-0.75 1.0])
line([0 0],get(gca,'ylim'),'Color','k','LineStyle','--')
xlabel('Time (sec)')
ylabel('Voltage (uV)')
title(figurename)
save_and_close_fig(figure_dir,figurename)