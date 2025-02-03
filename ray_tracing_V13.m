clear;
clc;
close all;
warning('off');



% changes made for x intercept = 0 are not done in this code as V9.
% update: dont know what changes coz code seems to work for x intercept=0


% note: checked if mu_medium1_medium2 use should it be * or \ in calculations?




%--------------------------------------------------------------------------
% wave properties
fluid_depth = 10; % fluid depth in mm
amp = 0.02; % amplitude of wave in mm
lambda = 7; % wavelength of wave in mm
tank_dim = 60; % tank length and width in mm
delta1 = 0.25;%01; % surface discretization
delta2 = 0.25; % light source discretization
screen_level = -550; % vertical distance between bottom and screen
mu_air_fluid = 1.33;    % refraction index for air-fluid
base_theta_i = 0; % base angle of source with flat fluid surface
% (considering depth and x_intercept)
x_intercept_list = -10:50:10; 
exp_decay_const = -0.05;%-0.05;
glass_th = 2; % thickness of glass
mu_fluid_glass = 1.13;
mu_air_glass = mu_air_fluid*mu_fluid_glass;

light_source_side_skip = 0.2;

x = -tank_dim/2:delta1:tank_dim/2;
y = -tank_dim/2:delta1:tank_dim/2;
z = fluid_depth*ones(size(x, 2), size(y, 2));




% loop over x intercepts for plotting purpose
for jjj=1:length(x_intercept_list)
    x_intercept = x_intercept_list(jjj);
    disp(x_intercept);
    close all;
        
    % loop to plot base level with 0 amplitude

    for gg=1:1
%         disp(gg);

        for i=1:length(x)
            for j=1:length(y)
                r = sqrt(x(i)^2+y(j)^2);
        %         if r<=tank_dim/2
        %             z(i, j) = z(i, j)+amp*cos(2*pi*r/lambda);
        %         else
        %             z(i, j) = 0;
        %         end
                % make sure the function definition is same in the cost function at
                % the bottom of the code
                if gg==1
                    z(i, j) = z(i, j)+amp*cos(2*pi*r/lambda)*exp(exp_decay_const*r);
                end  
            end
        end


        if gg==1
            surface(x, y, z, 'FaceAlpha', 1.0, 'EdgeColor', 'none');
%             surface(x, y, z, 'FaceColor', 'cyan', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            hold on;
        end


        % plot bottom of fluid
        bottom_level = 0; % keep this 0 so that depth defines the surface and bottom is alyways at 0
        glass_level = bottom_level-glass_th;

        if gg==1
            z_bottom = bottom_level*ones(size(x, 2), size(y, 2));
            z_glass = glass_level*ones(size(x, 2), size(y, 2));
            surface(x, y, z_bottom, 'FaceColor', 'cyan', 'FaceAlpha', 0.5, 'EdgeColor', 'cyan');
            surface(x, y, z_glass, 'FaceColor', 'cyan', 'FaceAlpha', 0.5, 'EdgeColor', 'cyan');

%             surface(x, y, z_bottom, 'FaceColor', 'black', 'FaceAlpha', 1.0, 'EdgeColor', 'black');
        end
        
        
        
        
        %--------------------------------------------------------------------------
        % light source

        % angle of incidence not explicitely defined for diverging light source
        
        
%         light_center_x = -174;%500;
%         light_center_y = 0;
        light_center_z_list = 650;%(fluid_depth+10):2:(fluid_depth+10)+0.5;
        % length of list and delta z_list in above and below needs to be same
        light_spread_z_list = fluid_depth:2:fluid_depth+0.5;%10:1:10;
        light_spread_x = x_intercept;%-10;%0;
        light_spread_y_list = -tank_dim/2+light_source_side_skip:delta2:tank_dim/2-light_source_side_skip;
       
        light_center_x = -1*tand(base_theta_i)*abs(light_center_z_list(1) - fluid_depth) + x_intercept;

        % old code for calculating base_theta_i, now we calculate
        % light_center_x
%         base_theta_i = atand(abs((light_center_x - x_intercept) / (light_center_z_list(1) - fluid_depth);

        % for storing angles of incideces and refractions for wave and
        % plane surface of fluid, in degrees
        if gg==1
            theta_i_wave = zeros(1, length(light_spread_y_list));
            theta_r_wave = zeros(1,length(light_spread_y_list));
        else
            theta_i_plane = zeros(1, length(light_spread_y_list));
            theta_i_plane = zeros(1, length(light_spread_y_list));
        end

        for jj=1:length(light_center_z_list)
            light_center_z = light_center_z_list(jj);
%             disp(light_center_z);


            %--------------------------------------------------------------------------
            % ray tracing in air

            surface_intersect_x = zeros(1, length(light_spread_y_list));
            surface_intersect_y = zeros(1, length(light_spread_y_list));
            surface_intersect_z = zeros(1, length(light_spread_y_list));

            glass_top_intersect_x = zeros(1, length(light_spread_y_list));     
            glass_top_intersect_y = zeros(1, length(light_spread_y_list));     
            glass_top_intersect_z = zeros(1, length(light_spread_y_list));
            
            glass_bottom_intersect_x = zeros(1, length(light_spread_y_list));     
            glass_bottom_intersect_y = zeros(1, length(light_spread_y_list));     
            glass_bottom_intersect_z = zeros(1, length(light_spread_y_list));

            screen_intersect_x = zeros(1, length(light_spread_y_list));     
            screen_intersect_y = zeros(1, length(light_spread_y_list));     
            screen_intersect_z = zeros(1, length(light_spread_y_list));     

            init_guess = 0;
            t_list = zeros(1, length(light_spread_y_list));

            for i=1:length(light_spread_y_list)
%                 disp(i);
%                 disp('y-coordinate');
%                 disp(light_spread_y_list(i));
                light_center_y = 0;%light_spread_y_list(i);
                light_source = [light_center_x, light_center_y, light_center_z];
                if gg==1
                    scatter3(light_center_x, light_center_y, light_center_z, '.');
                    hold on;
                end
        
                ray_direction_x = light_spread_x-light_center_x;
                ray_direction_z = light_spread_z_list(jj)-light_center_z_list(jj);

                
                ray_direction_y = light_spread_y_list(i)-light_center_y;
                ray_direction = [ray_direction_x, ray_direction_y, ray_direction_z];
                ray_direction = ray_direction/norm(ray_direction); % unit ray direction
%                 disp('ray direction');
%                 disp(ray_direction);
                if gg==1
                    wave_flag = true;
                else
                    wave_flag = false;
                end

                surface_intersect_para = @(t) cost_function(t, light_source, ray_direction, fluid_depth, amp, lambda, wave_flag, exp_decay_const);

                t = fzero(surface_intersect_para, init_guess);
%                 disp('t parameter');
%                 disp(t);
                t_list(i) = t;
                intersect = light_source+t*ray_direction;
%                 disp('intersect coordinates');
%                 disp(intersect);
%                 disp('------------------');

                surface_intersect_x(i) = intersect(1);
                surface_intersect_y(i) = intersect(2);
                surface_intersect_z(i) = intersect(3);
                if gg==1
                    plot3([light_center_x, surface_intersect_x(i)],...
                        [light_center_y, surface_intersect_y(i)],...
                        [light_center_z, surface_intersect_z(i)]);
                end
                init_guess = min(t_list(t_list>0));%[1.1*t 0.9*t];
            %     disp(init_guess);
                % for a fixed initial guess, the code finds intersect near to that
                % point even if the ray has a intersection before that, hopefully this
                % solves the problem since now the code searches for t value near to
                % the previously found value
            end

%             if gg==1
%                 scatter3(surface_intersect_x, surface_intersect_y, surface_intersect_z, '.');
%             end

            %--------------------------------------------------------------------------
            % refracted ray tracing in fluid

            unit_surface_normal_list = 1e6*ones(3, length(light_spread_y_list));

            for i=1:length(light_spread_y_list)
                light_center_y = light_spread_y_list(i);
%                 disp(i);
                % surface normal obtained from partial derivatives of parametric
                if gg==1
                    % surface eq F(x, y) = (x, y, depth+amp*cos(sqrt(x^2+y^2)*2pi/lambda))
                    r = sqrt(surface_intersect_x(i)^2+surface_intersect_y(i)^2);
                    kkk = -amp*(surface_intersect_x(i)/r)*exp(exp_decay_const*r)*(sin(2*pi*r/lambda)*2*pi/lambda...
                        -cos(2*pi*r/lambda)*exp_decay_const);
                    dF_dx = [1; 0; kkk];
                    
                    kkk = -amp*(surface_intersect_y(i)/r)*exp(exp_decay_const*r)*(sin(2*pi*r/lambda)*2*pi/lambda...
                        -cos(2*pi*r/lambda)*exp_decay_const);                   
                    dF_dy = [0; 1; kkk];
                    surface_normal = cross(dF_dx, dF_dy) / norm(cross(dF_dx, dF_dy));
                else
                    surface_normal = [0;0;1];
                end

                unit_surface_normal_list(:, i) = surface_normal;

                % calculate angle of surface normal with the vertical line to know the peaks of waves 
                cos_surface_normal_vertical_angle = dot(surface_normal,[0; 0; 1])/norm(surface_normal);
                surface_normal_vertical_angle = acos(cos_surface_normal_vertical_angle);
            %     disp(surface_normal_vertical_angle*180/pi);

            %     finding angle of incidence using dot product
            %     incindent ray vector is actually used in opposite direction because 
            %     thats how the angle of incidence is defined
                incident_ray = -1*[surface_intersect_x(i)-light_center_x;...
                                surface_intersect_y(i)-light_center_y;...
                                surface_intersect_z(i)-light_center_z];
%                 disp(incident_ray);
                cos_theta_i = dot(surface_normal,incident_ray)/...
                                (norm(surface_normal)*norm(incident_ray));
                            
%                 disp('air to fluid refraction');
                
                theta_i = acos(cos_theta_i);
                theta_i_deg = theta_i*180/pi;
%                 disp(theta_i_deg);

                theta_r = asin(sin(theta_i)/mu_air_fluid);
                theta_r_deg = theta_r*180/pi;
%                 disp(theta_r_deg);
                
%                 disp('calculated refractive index is');
%                 disp(sin(theta_i)/sin(theta_r));

                if gg==1
    %                 disp(gg);
                    theta_i_wave(i) = theta_i_deg;
                    theta_r_wave(i) = theta_r_deg;
                else
    %                 disp(gg);
                    theta_i_plane(i) = theta_i_deg;
                    theta_r_plane(i) = theta_r_deg;
                 end

            %     total internal reflection
            %     if theta_i_deg > 47
            %         disp('skipped');
            %         continue;
            %     end

                unit_incident_ray = incident_ray/norm(incident_ray); 
%                 disp(unit_incident_ray);
                unit_surface_normal = surface_normal;

                
                % translate back normal to origin for calculations
                A1 = [cos(theta_i), 1;1,cos(theta_i)];
                B1 = [-cos(theta_r); -cos(theta_i-theta_r)];
                coeff = A1\B1;
%                 disp(coeff);
%                 if i==ceil(length(light_spread_y_list)/2)
%                     disp(A1);
%                     disp(B1);
%                     disp(coeff);
%                 end

                % this is just to take care of theta_i=theta_r=0 case which
                % results in coeff = NaN as it cant find a unique solution
                if isnan(coeff)
                    coeff = [0; -1];
                end
                refracted_ray = coeff(1)*unit_incident_ray+coeff(2)*unit_surface_normal;
%                 disp(refracted_ray);

                % plot normals and refracted ray
%                 if gg==1
%                      surface_normal = surface_normal*5+...
%                         [surface_intersect_x(i); surface_intersect_y(i) ; surface_intersect_z(i)];
%     
%                     plot3([surface_normal(1);surface_intersect_x(i)],...
%                         [surface_normal(2);surface_intersect_y(i)],...
%                         [surface_normal(3);surface_intersect_z(i)]);   
%                 end

                scale = abs((bottom_level-surface_intersect_z(i))/refracted_ray(3));
                % scale the refracted ray to make all of them end at screen level
                refracted_ray_mod = refracted_ray*scale...
                    +  [surface_intersect_x(i); surface_intersect_y(i) ; surface_intersect_z(i)];
              %     scatter3(refracted_ray_mod(1),refracted_ray_mod(2),refracted_ray_mod(3), '.', 'r');

                glass_top_intersect_x(i) = refracted_ray_mod(1);
                glass_top_intersect_y(i) = refracted_ray_mod(2);
                glass_top_intersect_z(i) = refracted_ray_mod(3);

    %             if abs(surface_normal_vertical_angle*180/pi) < 0.5
                    if gg==1
                        plot3([refracted_ray_mod(1);surface_intersect_x(i)],...
                        [refracted_ray_mod(2);surface_intersect_y(i)],...
                        [refracted_ray_mod(3);surface_intersect_z(i)]);
            %         disp('air-fluid refraction');
            %         disp('angle of incidence = '); disp(theta_i_deg);
            %         disp('angle of refraction = '); disp(theta_r_deg);
                    end
    %             end

                %--------------------------------------------------------------------------
                
                % refraction from fluid to glass
                
                 % surface normal at bottom
                unit_surface_normal_glass_top = [0; 0; 1];

%                 disp('fluid to glass refraction');

                % finding angle of incidence
                cos_theta_i_glass_top = dot(unit_surface_normal_glass_top,-1*refracted_ray);%/...
            %                     (norm(unit_surface_normal_glass_top)*norm(-1*refracted_ray));
            %     disp(-1*refracted_ray); disp(norm(-1*refracted_ray));
                

                theta_i_glass_top = acos(cos_theta_i_glass_top);
%                 theta_i_glass_top_deg = theta_i_glass_top*180/pi;
%                 disp(theta_i_glass_top_deg);
                
                theta_r_glass_top = asin(sin(theta_i_glass_top)/mu_fluid_glass);
%                 theta_r_glass_top_deg = theta_r_glass_top*180/pi;
%                 disp(theta_r_glass_top_deg);
                
%                 disp('calculated refractive index is');
%                 disp(sin(theta_i_glass_top)/sin(theta_r_glass_top));
                
                A2 = [cos(theta_i_glass_top), 1;1,cos(theta_i_glass_top)];
                B2 = [-cos(theta_r_glass_top); -cos(theta_i_glass_top-theta_r_glass_top)];
                coeff = A2\B2;
                % this is just to take care of theta_i=theta_r=0 case which
                % results in coeff = NaN as it cant find a unique solution
                if isnan(coeff)
                    coeff = [0; -1];
                end
                
                % refracted_ray for fluid is incident ray for glass (in
                % opposite direction)
                double_refracted_ray = coeff(1)*(-1*refracted_ray)+coeff(2)*unit_surface_normal_glass_top;
            %     disp(double_refracted_ray);    
                
                scale = abs(glass_level/double_refracted_ray(3));
                % scale the refracted ray to make all of them end at screen level
                double_refracted_ray_mod = double_refracted_ray*scale...
                    +  [refracted_ray_mod(1); refracted_ray_mod(2); refracted_ray_mod(3)];

                if gg==1
                    scatter3(double_refracted_ray_mod(1),double_refracted_ray_mod(2),double_refracted_ray_mod(3), '.');
                else
                    scatter3(double_refracted_ray_mod(1),double_refracted_ray_mod(2),double_refracted_ray_mod(3), 'r', '.');
                end   
                
                
                glass_bottom_intersect_x(i) = double_refracted_ray_mod(1);
                glass_bottom_intersect_y(i) = double_refracted_ray_mod(2);
                glass_bottom_intersect_z(i) = double_refracted_ray_mod(3);


                % if abs(surface_normal_vertical_angle*180/pi) < 0.5
                    if gg==1
                        plot3([double_refracted_ray_mod(1);refracted_ray_mod(1)],...
                            [double_refracted_ray_mod(2);refracted_ray_mod(2)],...
                            [double_refracted_ray_mod(3);refracted_ray_mod(3)]);
                %         disp('fluid-air refraction');
                %         disp('angle of incidence = '); disp(theta_i_bottom_deg);
                %         disp('angle of refraction = '); disp(theta_r_bottom_deg);
                    end
    %             end
                
                %--------------------------------------------------------------------------
                
                % refraction back in air
                
                % surface normal at bottom
                unit_surface_normal_glass_bottom = [0; 0; 1];   
                
%                 disp('glass to air refraction');

                
                %finding angle of incidence
                cos_theta_i_glass_bottom = dot(unit_surface_normal_glass_bottom,-1*double_refracted_ray);%/...
            %                     (norm(unit_surface_normal_glass_bottom)*norm(-1*double_refracted_ray));
            %     disp(-1*double_refracted_ray); disp(norm(-1*double_refracted_ray));

                theta_i_glass_bottom = acos(cos_theta_i_glass_bottom);
%                 theta_i_glass_bottom_deg = theta_i_glass_bottom*180/pi;
%                 disp(theta_i_glass_bottom_deg);
                
                theta_r_glass_bottom = asin(sin(theta_i_glass_bottom)*mu_air_glass);
%                 theta_r_glass_bottom_deg = theta_r_glass_bottom*180/pi;
%                 disp(theta_r_glass_bottom_deg);                
                
%                 disp('calculated refractive index is');
%                 disp(sin(theta_i_glass_bottom)/sin(theta_r_glass_bottom));
                
                
                % translate back normal to origin for calculations
                A3 = [cos(theta_i_glass_bottom), 1;1,cos(theta_i_glass_bottom)];
                B3 = [-cos(theta_r_glass_bottom); -cos(theta_i_glass_bottom-theta_r_glass_bottom)];
                coeff = A3\B3;
                % this is just to take care of theta_i=theta_r=0 case which
                % results in coeff = NaN as it cant find a unique solution
                if isnan(coeff)
                    coeff = [0; -1];
                end
                triple_refracted_ray = coeff(1)*(-1*double_refracted_ray)+coeff(2)*unit_surface_normal_glass_bottom;
            %     disp(double_refracted_ray);                
            
                scale = abs(screen_level/triple_refracted_ray(3));
                % scale the refracted ray to make all of them end at screen level
                triple_refracted_ray_mod = triple_refracted_ray*scale...
                    +  [double_refracted_ray_mod(1); double_refracted_ray_mod(2); double_refracted_ray_mod(3)];

                if gg==1
                    scatter3(triple_refracted_ray_mod(1),triple_refracted_ray_mod(2),triple_refracted_ray_mod(3), '.');
                else
                    scatter3(triple_refracted_ray_mod(1),triple_refracted_ray_mod(2),triple_refracted_ray_mod(3), 'r', '.');
                end

                screen_intersect_x(i) = triple_refracted_ray_mod(1);
                screen_intersect_y(i) = triple_refracted_ray_mod(2);
                screen_intersect_z(i) = triple_refracted_ray_mod(3);            
            
    %             if abs(surface_normal_vertical_angle*180/pi) < 0.5
                    if gg==1
                        plot3([triple_refracted_ray_mod(1);double_refracted_ray_mod(1)],...
                            [triple_refracted_ray_mod(2);double_refracted_ray_mod(2)],...
                            [triple_refracted_ray_mod(3);double_refracted_ray_mod(3)]);
                %         disp('fluid-air refraction');
                %         disp('angle of incidence = '); disp(theta_i_bottom_deg);
                %         disp('angle of refraction = '); disp(theta_r_bottom_deg);
                    end
    %             end            

%             disp('--------')
            end
        end
    end

    view(90, 0);
%    set(gca,'visible','off');
%     view(0,0);
    % random plotting

    % figure;
    % plot(theta_i_wave);
    % hold on;
    % if gg==2
    %     plot(theta_i_plane);
    %     legend('wave', 'plane');
    % end
    % title('theta incident');
    % 
    % 
    % figure;
    % plot(theta_r_wave);
    % hold on;
    % if gg==2
    %     plot(theta_r_plane);
    %     legend('wave', 'plane');
    % end
    %     title('theta refracted');


    
    %--------------------------------------------------------------------------
    % INVERSE PROBLEM SOLUTION
    %--------------------------------------------------------------------------
    % solving inverse problem where direction of incident ray, location on
    % interface and screen are known, have to find surface normals
    % TODO if this works for thick light source

    tStart = tic;    
    
    disp('inverse algorithm starts');
    
    mean_surface_intersect_x = light_spread_x;%mean(surface_intersect_x);% light_spread_x
    mean_surface_intersect_z = fluid_depth; %mean(surface_intersect_z);
    mean_glass_top_intersect_z = bottom_level;
    mean_glass_bottom_intersect_z = glass_level;
    mean_screen_intersect_z = screen_level;

    y_base = light_spread_y_list;%-tank_dim/2:delta2:tank_dim/2;
        
    unit_surface_normal_calc_final_list = 1e6*ones(3, length(y_base));
    snell_check_1 = zeros(length(y_base));
    snell_check_2 = zeros(length(y_base));
    snell_check_3 = zeros(length(y_base));
    
    surface_normal_coplanar_error = 1e6*ones(length(y_base));
    
    glass_top_intersect_inv_temp = 1e6*ones(3, length(y_base));
    glass_bottom_intersect_inv_temp = 1e6*ones(3, length(y_base));
    
    glass_top_intersect_inv = 1e6*ones(3, length(y_base));
    glass_bottom_intersect_inv = 1e6*ones(3, length(y_base));



    for i=1:length(y_base)  % loop over point on screen
        disp(i);
        unit_surface_normal_calc_temp_list = 1e6*ones(3, length(y_base));
        for j=1:length(y_base)  % loop over point on surface
%             disp(j);
            screen_vec = [screen_intersect_x(i); screen_intersect_y(i); mean_screen_intersect_z];
            surface_vec = [mean_surface_intersect_x; y_base(j); mean_surface_intersect_z];

            guess = [0; 0; 0; 0];  % initial guess can be intersection of a straight line
%             fermat_function(guess, mean_glass_top_intersect_z, screen_vec, surface_vec, mu_air_fluid)
            fermat_function_instance = @(glass_intersect_guess) fermat_function(glass_intersect_guess,...
                mean_glass_top_intersect_z, mean_glass_bottom_intersect_z, screen_vec, surface_vec,...
                mu_air_fluid, mu_air_glass, mu_fluid_glass);
            
            if i==10 && j==10
                options = optimset('Display', 'off', 'MaxFunEvals', 5e3, 'MaxIter', 5e3, 'TolFun', 1e-12, 'TolX', 1e-12);
                glass_intersect = fminsearch(fermat_function_instance, guess, options);
            else
                options = optimset('Display', 'off', 'MaxFunEvals', 5e3, 'MaxIter', 5e3, 'TolFun', 1e-12, 'TolX', 1e-12);
                glass_intersect = fminsearch(fermat_function_instance, guess, options);
            end
            
%             options = optimoptions('fsolve','Display','off', 'FinDiffType', 'central',...
%                 'MaxIter', 4000, 'MaxFunEvals', 2000000, 'FunctionTolerance', 1e-12);
%             glass_intersect = fsolve(fermat_function_instance, guess, options);
            glass_bottom_vec = [glass_intersect(1:2); mean_glass_bottom_intersect_z];
            glass_top_vec = [glass_intersect(3:4); mean_glass_top_intersect_z];
            
            glass_bottom_intersect_inv_temp(:, j) = glass_bottom_vec';
            glass_top_intersect_inv_temp(:, j) = glass_top_vec';

%             disp(glass_bottom_vec);
%             disp(glass_top_vec);
            
            % verifying Snell's law
%             disp('air-glass interface check');
            incident_ray = screen_vec - glass_bottom_vec; % vector direction reversed to measure angle correct
            refracted_ray = glass_top_vec - glass_bottom_vec;
            theta_i_test = acosd(dot(incident_ray,[0; 0; -1])/norm(incident_ray));
%             disp(theta_i_test);
            theta_r_test = acosd(dot(refracted_ray,[0; 0; 1])/norm(refracted_ray));
%             disp(theta_r_test);
            snell_check_1(i, j) = sind(theta_i_test) / sind(theta_r_test);

%             disp('glass-fluid interface check');
            incident_ray = glass_bottom_vec - glass_top_vec; % vector direction reversed to measure angle correct
            refracted_ray = surface_vec - glass_top_vec;
            theta_i_test = acosd(dot(incident_ray,[0; 0; -1])/norm(incident_ray));
%             disp(theta_i_test);
            theta_r_test = acosd(dot(refracted_ray,[0; 0; 1])/norm(refracted_ray));
%             disp(theta_r_test);
            snell_check_2(i, j) = sind(theta_i_test) / sind(theta_r_test);

            %------------------------------------------
            
            
            
            
            
            
            
            
            
            
            
            
            
            light_center_y = light_spread_y_list(j);
            light_source = [light_center_x, light_center_y, light_center_z];
            
            source = transpose(light_source);
            unit_incident_vec = -1*(surface_vec-glass_top_vec)/norm(surface_vec-glass_top_vec);
            unit_refracted_vec = (source-surface_vec)/norm(source-surface_vec);
%             disp('unit incident and refracted vec'); disp(unit_incident_vec'); disp(unit_refracted_vec');
            
            coeff_guess = [1; 1]; 
            snell_function_instance = @(coeff_guess) snell_function(coeff_guess, unit_incident_vec, unit_refracted_vec, mu_air_fluid);
            options = optimoptions('fsolve','Display','off', 'FinDiffType', 'central',...
                'MaxIter', 4000, 'MaxFunEvals', 2000000, 'FunctionTolerance', 1e-12);
            coeff = fsolve(snell_function_instance, coeff_guess, options);
%             disp(coeff);
            surface_normal = coeff(1)*unit_incident_vec + coeff(2)*unit_refracted_vec;
            % surface normal obtained above is directed downwards, towards
            % incident ray, thus we are looking for opoosite of that
            unit_surface_normal = -1*surface_normal/norm(surface_normal);
    %         disp(unit_surface_normal);
            unit_surface_normal_calc_temp_list(:,j) = unit_surface_normal;

            % verifying Snell's law   
    %         disp('fluid-air check');
            theta_i_test = acosd(dot(unit_incident_vec, unit_surface_normal)); % vector direction reversed to measure angle correct
%             disp(theta_i_test);
            theta_r_test = acosd(dot(unit_refracted_vec, unit_surface_normal));
%             disp(theta_r_test);
            snell_check_3(i, j) = sind(theta_r_test) / sind(theta_i_test);


            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            % verify if the surface normal is a valid surface normal for an axi-symmetric case
            % i.e. it lies in the span of 2 vectors viz. vertical and radius vector
            center_vec = [0;0;fluid_depth]; % for current surface this is zero
            vertical_vec = [0;0;1]; % always vertical
            radius_vec = surface_vec - center_vec;
%             disp(radius_vec')
            vert_rad_cross = cross(vertical_vec, radius_vec);
            unit_vert_rad_cross = vert_rad_cross / norm(vert_rad_cross);

            % this can be considered like autocorrelation calculations
            surface_normal_coplanar_error(i, j)= abs(dot(unit_surface_normal, unit_vert_rad_cross));
    %         disp(surface_normal_coplanar_error(i, j));

    %         scatter3(screen_vec(1), screen_vec(2), screen_vec(3));
    %         scatter3(bottom_vec(1), bottom_vec(2), bottom_vec(3));
    %         scatter3(surface_vec(1), surface_vec(2), surface_vec(3));
%             disp('-----');
        end  
        [val, id] = min(surface_normal_coplanar_error(i, :));
%         disp(id); disp(val);
%         disp(unit_surface_normal_calc_temp_list(1:3, id));
        unit_surface_normal_calc_final_list(:, i) = unit_surface_normal_calc_temp_list(1:3, id);
%         disp("----");


        glass_bottom_intersect_inv(:, i) = glass_bottom_intersect_inv_temp(1:3, id);
        glass_top_intersect_inv(:, i) = glass_top_intersect_inv_temp(1:3, id);
% 
    end

    surface_normal_error1 = unit_surface_normal_list - unit_surface_normal_calc_final_list;
  

    
    
    

    %--------------------------------------------------------------------------
    % reconstrucnt the surface using the surface normals calculated
    % in experimental data, first need to SORT the y-coordinates of surface
    % intersects obtained from inverse problem along with corresponding unit
    % surface normals because for construcntion of continous surface, we need
    % unit surface data of continuous poits
       
    z_base = 1e6*ones(1, length(y_base));
    z_ref =  1e6*ones(1, length(y));
    surface_normals_calculated_points = 1e6*ones(3, length(y_base));
    
    if gg==1
        for jjj=1:length(y_base)
            r = sqrt(mean_surface_intersect_x^2 + y_base(jjj)^2);
            % make sure the function definition is same in the cost function at
            % the bottom of the code
            
            z_base(jjj) = fluid_depth + amp*cos(2*pi*r/lambda)*exp(exp_decay_const*r);
            
            kkk = -amp*(surface_intersect_x(i)/r)*exp(exp_decay_const*r)*(sin(2*pi*r/lambda)*2*pi/lambda...
            -cos(2*pi*r/lambda)*exp_decay_const);

            dF_dx = [1; 0; kkk];

            kkk = -amp*(surface_intersect_y(i)/r)*exp(exp_decay_const*r)*(sin(2*pi*r/lambda)*2*pi/lambda...
            -cos(2*pi*r/lambda)*exp_decay_const);                   

            dF_dy = [0; 1; kkk];
                    
            surface_normals_calculated_points(:, jjj) = cross(dF_dx, dF_dy) / norm(cross(dF_dx, dF_dy));
        end
    
     
        for kkkk=1:length(y)
            r = sqrt(mean_surface_intersect_x^2 + y(kkkk)^2);
            z_ref(kkkk)= fluid_depth + amp*cos(2*pi*r/lambda)*exp(exp_decay_const*r);
        end
    end
    
    surface_intersect_calculated = ones(3, length(y_base));
    surface_intersect_calculated(1, :) = mean_surface_intersect_x;%mean(surface_intersect_x);
    surface_intersect_calculated(3, :) = z_base(1);

    

%     for i=1:length(y_base)-1
% %         disp(i);
%     %     disp(unit_surface_normal_calc_final_list(i));
%         surface_normal_scale = 1 / unit_surface_normal_calc_final_list(3, i);
%     %     disp(unit_surface_normal_calc_final_list(3, i));
% %         dF_dy = -1*unit_surface_normal_calc_final_list(2, i)*surface_normal_scale;

%         this step of averaging dFdy is not valid adn any averaging must be done in numerical integral
%           evaluation step (like Simpsons one third rule)  
%         dF_dy = -1*(unit_surface_normal_calc_final_list(2, i)+unit_surface_normal_calc_final_list(2, i+1))*surface_normal_scale/2;
%           
%         surface_intersect_calculated(2, i) = y_base(i);
%         surface_intersect_calculated(3, i+1) = surface_intersect_calculated(3, i) + delta2 * dF_dy;
%     %     disp(surface_normal_scale);
%     %     disp(delta2*dF_dy);
%     %     disp('---')
    
%     end
  
    dF_dy_list = 1e6*ones(1, length(y_base));
    % simpsons integration method
    for i=1:length(y_base)
        surface_normal_scale = 1 / unit_surface_normal_calc_final_list(3, i);    
%         dF_dy = -1*unit_surface_normal_calc_final_list(2, i)*surface_normal_scale;
        dF_dy_list(i) = -1*unit_surface_normal_calc_final_list(2, i)*surface_normal_scale; 
        surface_intersect_calculated(2, i) = y_base(i);

    end

    for i=1:2:length(y_base)-2 
        simpons_derivative = dF_dy_list(i)+4*dF_dy_list(i+1)+dF_dy_list(i+2);
        surface_intersect_calculated(3, i+2) = surface_intersect_calculated(3, i) + delta2 * simpons_derivative / 3; 
        % linear interpolation for intermediate points
%         surface_intersect_calculated(3, i+1) = (surface_intersect_calculated(3, i+2) + surface_intersect_calculated(3, i))/2;
    end 
    
    % spline interpolation for intermediate points
    surface_intersect_calculated(3, 2:2:end-1) = spline(surface_intersect_calculated(2, 1:2:end), surface_intersect_calculated(3, 1:2:end), surface_intersect_calculated(2, 2:2:end-1));
    tEnd = toc(tStart);

    %----------------------------------------------------------------------
    
    % plotting results
 
    surface_normal_error2 = surface_normals_calculated_points - unit_surface_normal_calc_final_list;
      
    error_list = abs(surface_intersect_calculated(3,:) - z_base);
    error_avg = mean(error_list);
    error_avg_pc = mean(error_list)*100 / amp;
    error_max = max(error_list);
    error_max_pc = max(error_list)*100 / amp;
    
    
    figure;
    plot(y_base, surface_intersect_calculated(3,:));
    hold on;
%     plot(y_base, z_base);
    plot(y, z_ref);
    ylim([fluid_depth-1.3*amp fluid_depth+2*amp]);
    legend('inverse solution','true shape');%, 'Location','northeastoutside');
    ylabel('Wave Amplitude (mm)');
    xlabel('Y coordinate (mm)');
%     suptitle('Amplitude vs. Y-coordinate');
    msg1 = ['l=', num2str(lambda), ', amp=', num2str(amp), ', psi=', num2str(base_theta_i), ', dy=', num2str(delta2)];
    msg2 = ['x_int=', num2str(light_spread_x), ', err(avg)=' num2str(error_avg_pc) ', err(max)=', num2str(error_max_pc) ', tEnd=' num2str(tEnd) ', k=' num2str(exp_decay_const)];
%     title({'Amplitude vs. Y-coordinate'; msg1; msg2});
%     filename = ['./inverse_solution_plots/lambda = ', num2str(lambda), ' mm, amplitude = ', num2str(amp), ' mm, base_theta_i = ', num2str(base_theta_i),', x intercept = ', num2str(light_spread_x), ', delta y = ', num2str(delta2), '_spline.png'];
    
    filename = ['./inverse_solution_plots/' msg1 ',' msg2 '.eps'];
%     filename2 = ['./inverse_solution_plots/' msg1 ',' msg2 '2.eps'];

    
    set(gca,'FontSize',10);%,'FontName','Times');
% %     print(filename,'-depsc2');
% %     saveas(gcf, filename, 'epsc');
%     print(gcf, filename,'-depsc2','-r300');
 
end

surface_intersect_forward = [surface_intersect_x; surface_intersect_y; surface_intersect_z];
glass_top_intersect_forward = [glass_top_intersect_x; glass_top_intersect_y; glass_top_intersect_z];
glass_bottom_intersect_forward = [glass_bottom_intersect_x; glass_bottom_intersect_y; glass_bottom_intersect_z];
screen_intersect_forward = [screen_intersect_x; screen_intersect_y; screen_intersect_z];


err1 = abs(glass_top_intersect_inv - glass_top_intersect_forward);
err2 = abs(glass_bottom_intersect_inv - glass_bottom_intersect_forward);

% disp(max(max(err1)));
% disp(max(max(err2)));




















%--------------------------------------------------------------------------

% function definitions

 function y = cost_function(t, light_source, ray_direction, tank_depth, amp, lambda, wave_flag, exp_decay_const)
    light_source_x = light_source(1);
    light_source_y = light_source(2);
    light_source_z = light_source(3);
    ray_direction_x = ray_direction(1);
    ray_direction_y = ray_direction(2); 
    ray_direction_z = ray_direction(3);
    r = sqrt((light_source_x+ray_direction_x*t)^2+(light_source_y+ray_direction_y*t)^2);
    if wave_flag==true
%         disp(wave_flag);
        y = light_source_z + ray_direction_z*t - tank_depth - amp*cos(2*pi*r/lambda)*exp(exp_decay_const*r);
    else
%         disp(wave_flag);
        y = light_source_z+ray_direction_z*t - tank_depth;
    end
 end
 
 %--------------------------------------------------------------------------
 
 % time required by ray to travel from 1 point to other, Fermat's theorem
 % used as minimization of time, guess vector is a 2d, z coordinate added
 % inside function as it is known parameter
 
  
 % while writing paper I noticed that I used wrong cost function, I used
 % partial derivatives in vector form, but I could have used just travel
 % time as scalar valued just like snells function, still dont know how
 % these two are equivalent, maybe fzero currently used gets replaced by
 % argmin. so this could be just trick to solve numerically
 
 function y = fermat_function(glass_intersect_guess, glass_top_intersect_z, glass_bottom_intersect_z, screen_vec, surface_vec,  mu_air_fluid, mu_air_glass, mu_fluid_glass)
    guess_vec_bottom = [glass_intersect_guess(1:2); glass_bottom_intersect_z];
    guess_vec_top = [glass_intersect_guess(3:4); glass_top_intersect_z];
    dist1 = norm(screen_vec - guess_vec_bottom);
    dist2 = norm(guess_vec_bottom - guess_vec_top);
    dist3 = norm(guess_vec_top - surface_vec);
%     disp(guess_vec);
%     disp(dist1);
%     disp(dist2);
    y = (dist1/mu_air_fluid)+(dist2*mu_fluid_glass)+dist3;
%         disp(y);
 end
 
 %--------------------------------------------------------------------------
 
 % Using Snell's law to define a cost function for obtaining the surface normal
 % here incident ray is from fluid to air (reverse direction of forward
 % problem) and both rays are unit length
 
 function y = snell_function(coeff_guess, unit_incident_vec, unit_refracted_vec, mu_air_fluid)
%     disp(coeff_guess);
%     disp(incident_vec);
%     disp(refracted_vec);
    surface_normal = coeff_guess(1)*unit_incident_vec + coeff_guess(2)*unit_refracted_vec;
%     disp(surface_normal);
    unit_surface_normal = surface_normal / norm(surface_normal);
    theta_i = acos(dot(unit_surface_normal, unit_incident_vec));
    theta_r = acos(dot(-1*unit_surface_normal, unit_refracted_vec));
    y = mu_air_fluid * sin(theta_i) - sin(theta_r);
    y=abs(y);
%     disp(theta_i*180/pi);
%     disp(theta_r*180/pi);
%     disp(y);
 end
 
 %--------------------------------------------------------------------------