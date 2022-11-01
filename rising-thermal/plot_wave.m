function plot_wave()
    for i = 0:2970
        filename = ['Theta_' num2str(i) '.csv'];
        if exist(filename) == 2
            data = dlmread(filename);
            name = num2str(i);
            High.Theta.(genvarname(name)) = data;
        end
    end
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    filename = '202x202.gif';
    for i = 1:2970
        name = num2str(i);
        if isfield(High.Theta,(genvarname(name)))
            contourf(High.Theta.(genvarname(name)));
            colorbar;
            caxis([300 300.5]);
            %hold on
            %quiver(Coriolis.u.(genvarname(name)),Normal.v.(genvarname(name)));
            drawnow;
            % Capture the plot as an image 
            frame = getframe(h); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            % Write to the GIF File 
            if i == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
            else 
                imwrite(imind,cm,filename,'gif','WriteMode','append'); 
            end 
        end
    end
end

plot(cell2mat(struct2cell(Low.max_T))(3:end),'LineWidth',1);
hold on
plot(cell2mat(struct2cell(Mid.max_T))(3:end),'LineWidth',1);
plot(cell2mat(struct2cell(High.max_T))(3:end),'LineWidth',1);
legend('101 x 101','202 x 202','404 x 404');