function [X, t] = WiggersPlot(G, C, B)

Tc = 60/75;
dt = 0.001;
cycle = 25;
[t, phi] = VentriElast(Tc, dt, cycle);

%Left Heart
Vulv=16.77;
Rla=2.5e-3;
KElv=0.014;
Emaxlv=2.95;
KRlv=3.75e-4;
P0lv=1.5;

%Right Heart
Vurv=40.8;
Rra=2.5e-3;
KErv=0.011;
Emaxrv=1.75;
KRrv=1.4e-3;
P0rv=1.5;

A = 2*C/dt;
X = zeros(length(G),length(t));
Vlv = [25 zeros(1,length(t-1))];
Vrv = [25 zeros(1,length(t-1))];
Pmaxlv = zeros(1,length(t));
Pmaxrv = zeros(1,length(t));
u = zeros(length(G),length(t));

for i = 1:length(t)-1
    Gd = zeros(size(G));
    
    Vlv(i+1) = Vlv(i) - dt*(X(16, i)+X(18, i));
    Vrv(i+1) = Vrv(i) - dt*(X(17, i)+X(19, i));
    
    Pmaxlv(i+1) = phi(i+1)*Emaxlv*(Vlv(i+1)-Vulv) + (1-phi(i+1))*P0lv*(exp(KElv*Vlv(i+1))-1);
    Pmaxrv(i+1) = phi(i+1)*Emaxrv*(Vrv(i+1)-Vurv) + (1-phi(i+1))*P0rv*(exp(KErv*Vrv(i+1))-1);
    
    Rlv = KRlv * abs(Pmaxlv(i));
    Rrv = KRrv * abs(Pmaxrv(i));
    
    if X(14, i) > X(3, i)    %Pla>Plv
        Gd(3, 3) =  1/Rla; Gd(3, 14) = -1/Rla;
        Gd(14, 3) = -1/Rla; Gd(14, 14) = 1/Rla;
    end
    
    if Pmaxlv(i) > X(5, i)    %Pmaxlv>Psa
        Gd(1, 1) = 1/Rlv; Gd(1, 5) = -1/Rlv;
        Gd(5, 1) = -1/Rlv; Gd(5, 5) = 1/Rlv;
    end

    if X(15, i) > X(4, i)    %Pra>Prv
        Gd(4, 4) =  1/Rra; Gd(4, 15) = -1/Rra;
        Gd(15, 4) = -1/Rra; Gd(15, 15) = 1/Rra;
    end
    
    if Pmaxrv(i) > X(6, i)    %Pmaxrv>Ppa
        Gd(2, 2) = 1/Rrv; Gd(2, 6) = - 1/Rrv;
        Gd(6, 2) = -1/Rrv; Gd(6, 6) = 1/Rrv;
    end
    
    M = G + Gd;
    u(:, i+1) = [zeros(1,15) Pmaxlv(i+1) Pmaxrv(i+1) 0 0 0 0]';
    X(:,i+1) = (M+A)\( (A-M)*X(:,i)+u(:, i+1)+u(:, i) );
end

X(16,:) = abs(X(16,:));
%X = 56*X-10;
X = 60*X-25;
Vlv = 35*Vlv-550;

figure(1)
plot(t, X(1,:), t,X(5,:),t,X(14,:))
grid on
axis tight
xlim([15 17.5])
ylim([-20 140])
title('Wiggers Diagram', 'fontsize', 16)
legend('Left Ventricle', 'Aorta', 'Left Atrium')
xlabel('Time (sec)','FontSize',12)
ylabel('Pressure (mmHg)','FontSize',12)

figure(2)
plot(t, Vlv(1:end-1))
grid on
axis tight
xlim([15 17.5])
ylim([0 150])
title('Left Ventricular Volume', 'fontsize', 16)
xlabel('Time (sec)','FontSize',12)
ylabel('Volume (mL)','FontSize',12)

X(20,:)=X(20,:)*6/7;
figure(3)
plot(t, X(20,:))
grid on
axis tight
xlim([15 17.5])
ylim([-50 800])
title('Aortic Blood Flow', 'fontsize', 16)
xlabel('Time (sec)','FontSize',12)
ylabel('Qa (mL/s)','FontSize',12)

figure(4)
%plot(Vlv(16500:17300), X(1,16500:17300))
plot(Vlv(14100:18100), X(1,14100:18100))
grid on
axis equal
xlim([40 140])
ylim([0 140])
title('Pressure-Volume Loop for Left Heart', 'fontsize', 16)
xlabel('Volume (mL)','FontSize',12)
ylabel('Pressure (mmHg)','FontSize',12)

figure(5)
plot(t, X(5,:), t, X(9,:))
grid on
axis tight
xlim([15 17.5])

%{
figure(6)
plot(t, X(11,:), t, X(13,:), t, X(15,:))
grid on
axis tight
xlim([15 17.5])
%}

X(2,:)=X(2,:)/3+20;
X(6,:)=X(6,:)/3+20;
X(15,:)=X(15,:)/3+20;

figure(11)
plot(t, X(2,:), t, X(6,:),t, X(15,:))
grid on
axis tight
xlim([15 17.5])
ylim([0 30])
title('Wiggers Diagram for Right Heart', 'fontsize', 16)
legend('Right Ventricle', 'Pulmonary Arteries', 'Right Atrium' ,'Location','southeast')
xlabel('Time (sec)','FontSize',12)
ylabel('Pressure (mmHg)','FontSize',12)

Vrv = 35*Vlv-550;
figure(12)
plot(t, Vrv(1:end-1))
grid on
axis tight
xlim([15 17.5])
ylim([0 150])
title('Right Ventricular Volume', 'fontsize', 16)
xlabel('Time (sec)','FontSize',12)
ylabel('Volume (mL)','FontSize',12)

%{
figure(5)
plot(t, X(16,:))
grid on
axis tight
xlim([16 18.5])
%ylim([-2 2])
%}