function [t,EA,EAd] = fig_animate(tout,Xout,Uout,Xdout,Udout,Uext,p)

flag_movie = p.flag_movie;

if flag_movie
    
    try
        name = ['test.mp4'];
        vidfile = VideoWriter(name,'MPEG-4');
    catch ME
        name = ['test'];
        vidfile = VideoWriter(name,'Motion JPEG AVI');
    end
    open(vidfile);
end

%% smoothen for animation
t = (tout(1):p.simTimeStep:tout(end));
X = interp1(tout,Xout,t);
U = interp1(tout,Uout,t);
Xd = interp1(tout,Xdout,t);
Ud = interp1(tout,Udout,t);
Ue = interp1(tout,Uext,t);

%% loop through frames
figure('Position',[200 100 1000 600]);
set(0, 'DefaultFigureRenderer', 'opengl');
set(gcf, 'Color', 'white')
N = 3;M = 3;    % subplot size

nt = length(t);
EA = zeros(nt,3);
EAd = zeros(nt,3);
for ii = 1:nt
    EA(ii,:) = fcn_X2EA(X(ii,:));
    EAd(ii,:) = fcn_X2EA(Xd(ii,:));
end

% for ZOH force
t2 = repelem(t,2);
t2(1) = []; t2(end+1) = t2(end);
U2 = repelem(U,2,1);

%%%%%% subfigure handles %%%%%%%
h_x = subplot(N,M,3);
h_dx = subplot(N,M,6);
h_w = subplot(N,M,9);
h_u = subplot(N,M,[7,8]);

for ii = 1:p.playSpeed:nt
    %% The main animation
    % plot setting
    pcom = X(ii,1:3)';
    h_main = subplot(N,M,[1,2,4,5]);
    hold on; grid on;axis square;axis equal;
    h_main.XLim = [pcom(1)-0.5 pcom(1)+0.5];
    h_main.YLim = [pcom(2)-0.5 pcom(2)+0.5];
    h_main.ZLim = [-0.2 0.6];
    
    viewPt = [0.2,0.5,0.2];
    view(viewPt);
    
    % plot robot & GRF
    % real
    fig_plot_robot(X(ii,:)',U(ii,:)',Ue(ii,:)',p)
    % desired
    fig_plot_robot_d(Xd(ii,:)',0*Ud(ii,:)',p)
    
    % text
    txt_time = ['t = ',num2str(t(ii),2),'s'];
    text(pcom(1),pcom(2),0.4,txt_time)
    txt_vd = ['vd = ',num2str(Xd(ii,4),2),'m/s'];
    text(pcom(1),pcom(2),0.5,txt_vd)
    txt_v = ['v = ',num2str(X(ii,4),2),'m/s'];
    text(pcom(1),pcom(2),0.45,txt_v)

    %% states
    %%%%%%%%% position %%%%%%%%%
    plot(h_x,t(1:ii),X(1:ii,1),'r',...
         t(1:ii),X(1:ii,2),'g',...
         t(1:ii),X(1:ii,3),'b',...
         t(1:ii),Xd(1:ii,1),'r--',...
         t(1:ii),Xd(1:ii,2),'g--',...
         t(1:ii),Xd(1:ii,3),'b--','linewidth',1)
    h_x.XLim = [t(1) t(end)];
    set( get(h_x,'Title'), 'String', 'Position [m]');

    %%%%%%%%% velocity %%%%%%%%%
    plot(h_dx,t(1:ii),X(1:ii,4),'r',...
         t(1:ii),X(1:ii,5),'g',...
         t(1:ii),X(1:ii,6),'b',...
         t(1:ii),Xd(1:ii,4),'r--',...
         t(1:ii),Xd(1:ii,5),'g--',...
         t(1:ii),Xd(1:ii,6),'b--','linewidth',1)
    h_dx.XLim = [t(1) t(end)];
    set( get(h_dx,'Title'), 'String', 'Velocity [m/s]');

    %%%%%%%%% Angular velocity %%%%%%%%% 
    plot(h_w,t(1:ii),X(1:ii,16),'r',...
         t(1:ii),X(1:ii,17),'g',...
         t(1:ii),X(1:ii,18),'b',...
         t(1:ii),Xd(1:ii,16),'r--',...
         t(1:ii),Xd(1:ii,17),'g--',...
         t(1:ii),Xd(1:ii,18),'b--','linewidth',1)
    h_w.XLim = [t(1) t(end)];
    set( get(h_w,'Title'), 'String', 'Angular velocity [rad/s]' );
    
    
    %% control
    plot(h_u,t2(1:2*ii),U2(1:2*ii,3),'r',...
         t2(1:2*ii),U2(1:2*ii,6),'g',...
         t2(1:2*ii),U2(1:2*ii,9),'b',...
         t2(1:2*ii),U2(1:2*ii,12),'k',...
         t(1:ii),Ud(1:ii,3),'r--',...
         t(1:ii),Ud(1:ii,6),'g--',...
         t(1:ii),Ud(1:ii,9),'b--',...
         t(1:ii),Ud(1:ii,12),'k--','linewidth',1)
    h_u.XLim = [t(1) t(end)];
    set( get(h_u,'Title'), 'String', 'Fz [N]' );
    
    %% make movie
    if flag_movie
        writeVideo(vidfile, getframe(gcf));
    end
    
    drawnow
    if ii < nt
        cla(h_main);
    end
    
end

if flag_movie
    close(vidfile);
end




