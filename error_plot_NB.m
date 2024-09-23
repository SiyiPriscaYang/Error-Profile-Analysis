%% this script reads error vector files and draw them out
clear all;
clc;

th=100;%% the maximum number of error vectors you want to plot
% Specify the file name
%% this file reads the nonbinary parity check matrix 
% (0,1,2,3 for I,X,Y,Z, respectively)
H=load('H_3_7_m3.txt');
% you can also load Hx and Hz directly if you have the parity check matrix
Hz=(H>1);
Hx=mod(H-Hz,2);
%% change the name to the file that contains your error vector
% this contains nonbinary error vectors
% each row: a error vector appended with the error type (1,2 for
% convergent/non-convergent errors, respectively. 
% Comment line 36 out if your error vector file does not have this information)
% (0,1,2,3 for I,X,Y,Z, respectively)
% format: 
E = load('nb_deg_3_7_m3_0.04_500.txt');
Sizex = size(Hx);
Sizez = size(Hz);
%% change line 31-34, line 51-53 if you have other formats of error vector,
% just need to make sure ex and ez in 52-53 corresponds to the binary
% vector for X/Z errors.

nCx = Sizex(1);
nCz = Sizez(1);
nV = Sizex(2);

SizeE = size(E);
nE = floor(SizeE(1));
err_type=E(:,nV+1);
%% comment out this line if your error vector does not have anythig appended to it
E=E(:,1:(nV));

% Define the radius of the circles
R = 10;
r = 0.3;
ro=3;
% pauli[0][0]=I, pauli[1][0]=X, pauli[0][1]=Z, pauli[1][1]=Y;  
pauli={'I','Z';'X','Y'};
pcolor={'k','b';'r','g'};
w=[0.4;0.6];
lwidth=1;

for n=1:min(th,nE)
    if err_type(n)==2
    %% Create a figure
    figure;

    %% VNs
    e=E(n,:);
    ez=(e>1);% z error vector, checked by Hx
    ex=mod(e-ez,2);% x error vector, checked by Hz
%i,x,y,z, 0,1,2,3
%ez       0,0,1,1
%e-ez     0,1,1,2
%ex       0,1,1,0
    ind_ex=find(ex);% qbits with X errors
    ind_ez=find(ez);% qbits with Z errors
    ind=[ind_ex,ind_ez];% errorneous qbits
    ext=ex(ind)';
    ezt=ez(ind)';
    k=length(ind);% number of errorneous qbits (VNs)
    
    % Calculate the angles for equally spaced points
    theta = linspace(0, 2*pi, k + 1);
    % Calculate the coordinates of the small circle
    pos_vn_x = R * cos(theta);
    pos_vn_y = R * sin(theta);

    % Plot VNs (small circles) on the large circle with labels
    for i = 1:k
        idx=ind(i);
        % Plot the small circle with a label
        
        %smallCircle = viscircles([pos_vn_x(i), pos_vn_y(i)], r, 'EdgeColor', pcolor{ex(idx),ez(idx)}, 'Label', strcat(num2str(idx),pauli{ex(idx),ez(idx)}));
        rectangle('Position', [pos_vn_x(i)-r, pos_vn_y(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor', pcolor{ex(idx)+1, ez(idx)+1},'FaceColor', pcolor{ex(idx)+1, ez(idx)+1},'LineWidth', lwidth);
        hold on;
        % Add label within the circle
        %text(pos_vn_x(i), pos_vn_y(i), num2str(idx), 'Color', pcolor{ex(idx)+1, ez(idx)+1}, 'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                
        %hold on;
    end

    %% CNs
    Hz_t=Hz(:,ind);
    Hx_t=Hx(:,ind);

    deg_cz=sum(Hz_t,2);
    deg_cx=sum(Hx_t,2);

    % degree>=2 neighbors of errorneous qubits
    ind_cz_in=find(deg_cz>1);
    ind_cx_in=find(deg_cx>1);
    n_cz_in=length(ind_cz_in);
    n_cx_in=length(ind_cx_in);

    % calculate position of degree>=2 CNs
    Hz_d2=Hz_t(ind_cz_in,:);
    Hx_d2=Hx_t(ind_cx_in,:);
    %pos_cz_d2=zeros(n_cz_in,2);
    %pos_cx_d2=zeros(n_cx_in,2);
    sz=mod(Hz_d2*ext,2);
    sx=mod(Hx_d2*ezt,2);
    for i=1:n_cz_in
        cz_adj_vn=find(Hz_d2(i,:));
        if length(cz_adj_vn)>2
            posx=mean(pos_vn_x(cz_adj_vn));
            posy=mean(pos_vn_y(cz_adj_vn));
            if sz(i)
                rectangle('Position', [posx-r, posy-r, 2*r, 2*r], 'EdgeColor', 'b','FaceColor', 'b','LineWidth', lwidth);
                hold on;
            else
                rectangle('Position', [posx-r, posy-r, 2*r, 2*r], 'EdgeColor', 'b','LineWidth', lwidth);
                hold on;
            end
        else
            posx=pos_vn_x(cz_adj_vn)*w;
            posy=pos_vn_y(cz_adj_vn)*w;
            rectangle('Position', [posx-r, posy-r, 2*r, 2*r], 'EdgeColor', 'b','LineWidth', lwidth);
            hold on;
        end
        % Plot the rectangle
        
        % Save label for legend
        %legendStrings{i} = num2str(ind_cz_in(i));
        % Plot edges connecting rectangle to adjacent circles
        for j = cz_adj_vn
            plot([posx, pos_vn_x(j)], [posy, pos_vn_y(j)], 'Color', 'b', 'LineWidth', lwidth);
            hold on;
        end
    end

    for i=1:n_cx_in
        cx_adj_vn=find(Hx_d2(i,:));
        if length(cx_adj_vn)>2
            posx=mean(pos_vn_x(cx_adj_vn));
            posy=mean(pos_vn_y(cx_adj_vn));
            if sx(i)
                rectangle('Position', [posx-r, posy-r, 2*r, 2*r], 'EdgeColor', 'r','FaceColor', 'r','LineWidth', lwidth);
                hold on;
            else
                rectangle('Position', [posx-r, posy-r, 2*r, 2*r], 'EdgeColor', 'r','LineWidth', lwidth);
                hold on;
            end
        else
            posx=pos_vn_x(cx_adj_vn)*w;
            posy=pos_vn_y(cx_adj_vn)*w;
            rectangle('Position', [posx-r, posy-r, 2*r, 2*r], 'EdgeColor', 'r','LineWidth', lwidth);
            hold on;
        end
        % Save label for legend
        %legendStrings{i} = num2str(ind_cx_in(i));
        % Plot edges connecting rectangle to adjacent circles
        for j = cx_adj_vn
            plot([posx, pos_vn_x(j)], [posy, pos_vn_y(j)], 'Color', 'r', 'LineWidth', lwidth);
            hold on;
        end
    end

    % degree 1 neighbors of errorneous qubits
    ind_cz_out=find(deg_cz==1);
    ind_cx_out=find(deg_cx==1);
 
    % calculate position of degree 1 CNs
    Hz_d1=Hz_t(ind_cz_out,:);
    Hx_d1=Hx_t(ind_cx_out,:);

    for i=1:k
        adj_cz=find(Hz_d1(:,i));
        adj_cx=find(Hx_d1(:,i));
        nz=length(adj_cz);
        nx=length(adj_cx);
        nt=nz+nx;
        if nt<5
            gamma = linspace((1-nt)*pi/12, (nt-1)*pi/12, nt);
        else
            gamma = linspace(-pi/3, pi/3, nt);
        end
        posx=pos_vn_x(i)+ro*cos(gamma+theta(i));
        posy=pos_vn_y(i)+ro*sin(gamma+theta(i));
        for j=1:nt
            if j>nz
                tcolor='r';
                idx=ind_cx_out(j-nz);
                err=ezt(i);
            else
                tcolor='b';
                idx=ind_cz_out(j);
                err=ext(i);
            end
                    
            % Display points
            %scatter(posx(j), posy(j), 50, tcolor, 'filled', 'DisplayName', num2str(idx));
            %hold on;
            
            if err
                rectangle('Position', [posx(j)-r, posy(j)-r, 2*r, 2*r], 'EdgeColor', tcolor,'FaceColor', tcolor, 'LineWidth', lwidth);
                hold on;
            else
                rectangle('Position', [posx(j)-r, posy(j)-r, 2*r, 2*r], 'EdgeColor', tcolor,'LineWidth', lwidth);
                hold on;
            end
            % Display lines
            plot([posx(j), pos_vn_x(i)], [posy(j), pos_vn_y(i)], 'Color', tcolor, 'LineWidth', lwidth);
            hold on;
        end
    end
    
    % Set axis limits for better visualization
    axis equal;
    %set(gca,'XTick',[], 'YTick', []);
    axis off;
    %legend(legendStrings, 'Location', 'Best');
    end
end
