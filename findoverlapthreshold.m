function [MuI_threshhold_all_networks] = findoverlapthreshold(eta_to_template_vox,network_names,Zscore_eta, method)
addpath(genpath('/home/faird/shared/code/internal/utilities/plotting-tools'));
% This function works by taking in a matrix of network X voxtel  of eta squared
% values and sets a threshold by finding the local minimum
% input is matrix with the voxels as rows and
% output is

%HARDCODE WARNING - but no other method is currently available.
%method = 'hist_localmin';
method = 'smooth_then_derivative';
plot_data = 0; % Originally was 1, changed by Cristian Morales to 0 (to run as job on MSI, in case is 1, requires an interactive session)
%network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'};


%etaZabove1 = etaZ(:,j) > thresholds(k);
%overlapetaZ(:,j,k) = etaZabove1*j;
if exist('plot_data','var') && plot_data == 1
    F = figure();
    G = figure();
else
end


for j=1:length(network_names)
    if j==4 || j ==6
        continue
    end
    
    if exist('Zscore_eta','var') && Zscore_eta == 1
        etaZ(:,j) = zscore(eta_to_template_vox(:,j));
        X =etaZ;
    else
        X = eta_to_template_vox(:,j);
    end
    
    
    %% Encapsulate data as cell
    
    Xcell=encapsulate_X(X);
    
    %% Count elements per bin
    
    
    options.shown_as='stairs';
    options.n_bins = 100;
    options.LineWidth=1.5;
    % options.n_bins=211;
    my_color=[27,158,119
        217,95,2
        117,112,179
        231,41,138]/255;
    
    ct{1}='box';
    ct{2}='curve';
    ct{3}='stairs';
    ct{4}='contour';
    options.n_bins=[];
    options.xlim=[];
    [x_hist, y_hist]=get_hist_data(Xcell,options);
    
    dx=mean(diff(x_hist));
    x=x_hist;
    x=[x(1)-dx; x(:); x(end)+dx];
    y=[0; y_hist(:); 0];
    
    n_points=1e4;
    xx=linspace(x(1),x(end),n_points);
    yy=spline(x,y,xx);
    yy(yy<0)=0;
    
    switch method
        
        case 'hist_localmin'
            %invert yy to find min signal.
            inyy = max(yy(:)) - yy;
            [peaks, locs] = findpeaks(inyy,'MinPeakProminence',100);
            [locs, idx] = max(locs); % take the max in case locs find to peaks.
            peaks = peaks(idx);
            MuI_threshhold = xx(locs);
            networkthreshold(j)=locs;
            MuI_threshhold_all_networks(j) = MuI_threshhold;
            
            %plot data
            if exist('plot_data','var') && plot_data == 1
                subplot(4,4,j)
                plot(xx,yy);hold on;
                plot(xx(locs),yy(locs),'o');
                axis([0 0.6 0 7500]);
                title([network_names{j}])
            else
            end
            
            
        case 'smooth_then_derivative'
            dyy=diff(yy);
            %resp_trace=yy;
            
            %x=filloutliers(resp_trace,'linear','movmedian',100);
            %smoothyy=smoothdata(yy,'sgolay',2000)
            smoothyy=smooth(yy,2000,'sgolay');
            smoothdyy=diff(smoothyy);
            
            %[foo ix_peak]=min(abs(smoothdyy(4000:end));
            %limit the minimium to after 4000 bins and before 7000 bins.
            [~,locs]=min(abs(smoothdyy(4000:7000)));
            locs = locs +4000;
            
            
            MuI_threshhold = xx(locs);
            networkthreshold(j)=locs;
            MuI_threshhold_all_networks(j) = MuI_threshhold;
            
            if exist('plot_data','var') && plot_data == 1
                figure(1)
                subplot(4,4,j)
                plot(xx,smoothyy)
                hold all
                plot(xx,yy)
                plot(xx(locs),yy(locs),'or')
                axis([0 0.65 0 6000]);
                hold off
                title([network_names{j}])
                %legend('smoothed','Orig')
                
                figure(2)
                subplot(4,4,j)
                plot(dyy)
                hold all
                line(xlim,[0 0])
                plot(smoothdyy)
                plot(locs,smoothdyy(locs),'or')
                hold off
                title([network_names{j}])
            else
            end
            
        case 'etaZ'
            thresholds = [1:0.25:3.5];
            for k=1: size(thresholds,2)
                for j=1:length(network_names)
                    if j==4 || j ==6
                        continue
                    end
                    
                    %etaZ(:,j) = zscore(eta_to_template_vox(:,j));
                    etaZabove1 = etaZ(:,j) > thresholds(k);
                    overlapetaZ(:,j,k) = etaZabove1*j;
                end
            end
            
        otherwise
            
    end %switch
end %go through every network
end %function