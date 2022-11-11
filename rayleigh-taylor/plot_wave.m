function plot_wave()
    for i = 0:50:5000
        filename = ['Theta_' num2str(i) '.csv'];
        if exist(filename) == 2
            data = dlmread(filename);
            name = num2str(i);
            RK4_HighRes.Theta.(genvarname(name)) = data;
        end
    end
    for i = 0:50:5000
        filename = ['Omega_' num2str(i) '.csv'];
        if exist(filename) == 2
            data = dlmread(filename);
            name = num2str(i);
            RK4_HighRes.Omega.(genvarname(name)) = data;
        end
    end
    for i = 0:50:5000
        filename = ['Streamfunction_' num2str(i) '.csv'];
        if exist(filename) == 2
            data = dlmread(filename);
            name = num2str(i);
            RK4_HighRes.Streamfunction.(genvarname(name)) = data;
        end
    end
    for k = 0:10:5000
        name = num2str(k);
        if isfield(Test5.Streamfunction,(genvarname(name)))
            for j = 2:400
                for i = 2:100
                    Test5.v.(genvarname(name))(j,i) = (Test5.Streamfunction.(genvarname(name))(j+1,i) - Test5.Streamfunction.(genvarname(name))(j-1,i))/0.0025;
                    Test5.u.(genvarname(name))(j,i) = (Test5.Streamfunction.(genvarname(name))(j,i+1) - Test5.Streamfunction.(genvarname(name))(j,i-1))/0.0025;
                    Test5.v.(genvarname(name))(401,i) = 0;
                    Test5.u.(genvarname(name))(401,i) = 0;
                end
                Test5.v.(genvarname(name))(j,101) = 0;
                Test5.u.(genvarname(name))(j,101) = 0;
            end
        end
    end
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    gifname = 'RK4_resolutions.gif';
    for i = 0:10:1000
        name1 = num2str(i);
        name2 = num2str(i*5);
        if isfield(AB3_LowRes.Theta,(genvarname(name1)))
            subplot(1,2,1);
            contourf(AB3_LowRes.Theta.(genvarname(name1)));
            title('AB3: 101x401');
            colorbar;
            caxis([299 300]);
            subplot(1,2,2);
            contourf(AB3_HighRes.Theta.(genvarname(name2)));
            title('AB3: 201x801');
            colorbar;
            caxis([299 300]);
            drawnow;
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if i == 0
                imwrite(imind,cm,gifname,'gif','Loopcount',inf);
            else
               imwrite(imind,cm,gifname,'gif','WriteMode','append'); 
        
            end
        end
    end
    subplot(1,2,2);
    for i = 0:10:1000
        name = num2str(i);
        if (isfield(AB3_LowRes.Theta,(genvarname(name))))
            contourf(AB3_LowRes.Theta.(genvarname(name)));
            title('AB3: 101x401');
            colorbar;
            caxis([299 300]);
            drawnow;
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,gifname,'gif','WriteMode','append');
        end
    end
    
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    gifname = 'low_resolution.gif';
    subplot(1,2,1);
    for i = 0:20:990
        name = num2str(i);
        if isfield(Test4.Theta,(genvarname(name)))
            contour(Test4.Theta.(genvarname(name)));
            title('101x401, AB3');
            colorbar;
            caxis([299 300]);
            drawnow;
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if i == 0
                imwrite(imind,cm,gifname,'gif','Loopcount',inf);
            else
               imwrite(imind,cm,gifname,'gif','WriteMode','append'); 
        
            end
        end
    end
    subplot(1,2,2);
    for i = 0:20:990
        name = num2str(i);
        if (isfield(Test6.Theta,(genvarname(name))))
            contour(Test6.Theta.(genvarname(name)));
            title('101x401, RK4');
            colorbar;
            caxis([299 300]);
            drawnow;
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,gifname,'gif','WriteMode','append');
        end
    end
end