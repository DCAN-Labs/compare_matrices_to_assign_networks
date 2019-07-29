
dyy=diff(yy);

subplot 211
plot(yy)


subplot 212
plot(dyy)
line(xlim,[0 0])

%%
subplot 211
%resp_trace=yy;

%x=filloutliers(resp_trace,'linear','movmedian',100);
smoothyy=smoothdata(yy,'sgolay',3000);
smoothdyy=diff(smoothyy);

%[foo ix_peak]=min(abs(smoothdyy(4000:end));
%limit the minimium to after 4000 bins.
[~,ix_peak]=min(abs(smoothdyy(4000:end)));
ix_peak = ix_peak +4000;

plot(smoothyy)
hold all
plot(yy)
plot(ix_peak,yy(ix_peak),'or')
hold off
legend('smoothed','Orig')



subplot 212
plot(dyy)
hold all
line(xlim,[0 0])
plot(smoothdyy)
plot(ix_peak,smoothdyy(ix_peak),'or')
hold off
