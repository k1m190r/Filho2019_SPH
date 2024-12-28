%-----------------------------------------------------------------
% Instructions:
% 1.Within the working directory containing the FORTRAN program,
% create and open a new file in the MATLAB Editor˜
% 2.Copy and paste this script and save as ‘‘name.m''
% 3.Create a subfolder called figures within the output folder
% 4.Update the path name of the working directory in this MATLAB script
% 5.Run
% 6.Figures will be save within the figures subfolder
%-----------------------------------------------------------------
%-----------------------------------------------------------------
%----------------- PLOTTING GRAPH (SERIES SOLUTION)---------------
%-----------------------------------------------------------------
clc
clear all
close all
a = load('working_directory/output/GEOMETRY.DAT');
Xm = a(1, 1);
Xmx = a(1, 2);
Ym = a(2, 1);
Ymx = a(2, 2);
b = load('working_directory/N_PART_SERIES.DAT');
part_side = sqrt(b(1, 1));
dt = b(2, 1);
c = load('working_directory/NUMBER_OF_ITERATIONS.DAT');
iterations = c(1, 1);
d = load('working_directory/STEP_OUT.DAT');
step_out = d(1, 1);
f = load ('working_directory/output/VERIFYING_TEMP_SERIES.DAT');
fr = sortrows(f);
x_plot = fr(:, 2);
y_plot = fr(:, 3);
T1 = fr(:, 4);
x_series = reshape(x_plot, part_side, part_side);
y_series = reshape(y_plot, part_side, part_side);
T_series = reshape(T1, part_side, part_side);
figure(1);
contourf(x_series, y_series, T_series, 20);
grid on
daspect ([1 1 1]);
colormap jet;
cH = colorbar;
set(gcf, 'Position', get (0, 'Screensize'));
set(cH, 'FontSize', 12);
set(get(cH, 'title'), 'string', 'T(ˆ{o}C)', 'FontSize', 12);
caxis([0 100]);
axis ([Xm Xmx Ym Ymx]);
set(gca, 'xticklabel', num2str(get(gca, 'ytick')', '%.2f'), 'fontsize', 10);
set(gca, 'yticklabel', num2str(get(gca, 'ytick')', '%.2f'), 'fontsize', 10);
title('SOLUTION BY SERIES', 'fontsize', 14);
xlabel('x (m)', 'fontsize', 12);
ylabel('y (m)', 'fontsize', 12);
disp('Figures will be saved in working_directory/output/figures')
disp('Press any key to continue...');
pause
dir_out = (['cd working_directory/output/figures']);
eval(dir_out);
img = getframe(gcf);
imwrite(img.cdata, ['SERIES_SOLUTION.png']);
close;
%-----------------------------------------------------------------
%----------------PLOTTING GRAPHS (SPH SOLUTIONS)------------------
%-----------------------------------------------------------------
dir_in = (['cd working_directory/output/temperature']);
%----------------------INITIAL DISPOSITION------------------------
FileName = (['0000001']);
t = load ('working_directory/output/temperature/0000001.dat');
fr = sortrows(t);
x_plot = fr(:, 2);
y_plot = fr(:, 3);
T1 = fr(:, 4);
x_SPH = reshape(x_plot, part_side, part_side);
y_SPH = reshape(y_plot, part_side, part_side);
T_SPH = reshape(T1, part_side, part_side);
eval(dir_out);
figure(2);
contourf(x_SPH, y_SPH, T_SPH, 20);
grid on
daspect ([1 1 1]);
colormap jet;
cH = colorbar;
set(gcf, 'Position', get (0, 'Screensize'));
set(cH, 'FontSize', 12);
set(get(cH, 'title'), 'string', 'T(ˆ{o}C)', 'FontSize', 12);
caxis([0 100]);
axis ([Xm Xmx Ym Ymx]);
set(gca, 'xticklabel', num2str(get(gca, 'ytick')', '%.2f'), 'fontsize', 10);
set(gca, 'yticklabel', num2str(get(gca, 'ytick')', '%.2f'), 'fontsize', 10);
title1 = (['SPH SOLUTION - Cubic Spline Kernel - t = 0.00 s']);
title(title1, 'Units', 'Normalized', 'fontsize', 12)
title(title1, 'fontsize', 12, 'fontweight', 'b')
xlabel('x (m)', 'fontsize', 12);
ylabel('y (m)', 'fontsize', 12);
img = getframe(gcf);
imwrite(img.cdata, [FileName, '.png']);
clear FileName;

for i = step_out:step_out:iterations,
    num = int2str(i);
    char1 = num2str(i);
    char2 = num2str(iterations);
    eval(dir_in);

    if (i $ > $ = 0) & (i $ < $ 10)
        FileName = (['000000' num]);
    else

        if (i $ > $ = 10) & (i $ < $ 100)
            FileName = (['00000' num]);
        else

            if (i $ > $ = 100) & (i $ < $ 1000)
                FileName = (['0000' num]);
            else

                if (i $ > $ = 1000) & (i $ < $ 10000)
                    FileName = (['000' num]);
                else

                    if (i $ > $ = 10000) & (i $ < $ 100000)
                        FileName = (['00' num]);
                    else

                        if (i $ > $ = 100000) & (i $ < $ 1000000)
                            FileName = (['0' num]);
                        else
                            FileName = ([num]);
                        end

                    end

                end

            end

        end

    end

    disp(['opening' FileName '.dat']);
    t = load ([FileName '.dat']);
    fr = sortrows(t);
    x_plot = fr(:, 2);
    y_plot = fr(:, 3);
    T1 = fr(:, 4);
    clear FileName;
    x_SPH = reshape(x_plot, part_side, part_side);
    y_SPH = reshape(y_plot, part_side, part_side);
    T_SPH = reshape(T1, part_side, part_side);
    figure(i);
    contourf(x_SPH, y_SPH, T_SPH, 20);
    grid on
    daspect ([1 1 1]);
    colormap jet;
    cH = colorbar;
    set(gcf, 'Position', get (0, 'Screensize'));
    set(cH, 'FontSize', 12);
    set(get(cH, 'title'), 'string', 'T(ˆ{o}C)', 'FontSize', 12);
    caxis([0 100]);
    axis ([Xm Xmx Ym Ymx]);
    set(gca, 'xticklabel', num2str(get(gca, 'ytick')', '%.2f'), 'fontsize', 10);
    set(gca, 'yticklabel', num2str(get(gca, 'ytick')', '%.2f'), 'fontsize', 10);
    aux2 = i * dt;
    num1 = num2str(aux2);
    title1 = (['SPH SOLUTION - t = ' num1 ' s']);
    title(title1, 'Units', 'Normalized', 'fontsize', 12)
    title(title1, 'fontsize', 12, 'fontweight', 'b')
    xlabel('x (m)', 'fontsize', 12);
    ylabel('y (m)', 'fontsize', 12);
    eval(dir_out);
    img = getframe(gcf);

    if i $ < $ 10
        nome = (['000000' num '.png']);
        saveas(gcf, nome);
    end

    if (i $ > $ = 10) & (i $ < $ 100)
        nome = (['00000' num '.png']);
        imwrite(img.cdata, nome);
    end

    if (i $ > $ = 100) & (i $ < $ 1000)
        nome = (['0000' num '.png']);
        imwrite(img.cdata, nome);
    end

    if (i $ > $ = 1000) & (i $ < $ 10000)
        nome = (['000' num '.png']);
        imwrite(img.cdata, nome);
    end

    if (i $ > $ = 10000) & (i $ < $ 100000)
        nome = (['00' num '.png']);
        imwrite(img.cdata, nome);
    end

    if (i $ > $ = 100000) & (i $ < $ 1000000)
        nome = (['0' num '.png']);
        imwrite(img.cdata, nome);
    end

    if (i $ > $ = 1000000)
        nome = ([num '.png']);
        imwrite(img.cdata, nome);
    end

    close all
end
