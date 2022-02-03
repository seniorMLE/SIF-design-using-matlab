clc
clear all
close all
data = dlmread('test2.txt', ' ');  % read data
endp=length(data);
% set_default_properties
format long;
%gyro measurements [cm deg/s]
Wb(:,1)=data(:,4).*((pi/180)/100);
Wb(:,2)=data(:,5).*((pi/180)/100);
Wb(:,3)=data(:,6).*((pi/180)/100);
%accelerometrs measurements [mg]
Fb(:,1)=data(:,7).*(9.81/1000);
Fb(:,2)=data(:,8).*(9.81/1000);
Fb(:,3)=data(:,9).*(9.81/1000);
%navigation angles MIDG solution [cdeg]
yawM=data(:,10).*((pi/180)/100);
pitchM=data(:,11).*((pi/180)/100);
rollM=data(:,12).*((pi/180)/100);
%Psition MIDG Solution ecef [cm]
PM(:,1)= data(:,15)./100;
PM(:,2)= data(:,16)./100;
PM(:,3)= data(:,17)./100;
[Plat_M, Plon_M, Palt_M] = ecef2lla(PM(:,1),PM(:,2),PM(:,3));
%Velocity MIDG Solution ecef[cm/s]
VM(:,1) = data(:,18)./100;
VM(:,2) = data(:,19)./100;
VM(:,3)= data(:,20)./100;
%Psition and Velocity GPS ecef [cm] [cm/s]
GPS_flag=data(:,30);
PG(:,1) = data(:,34)./100;
PG(:,2) = data(:,35)./100;
PG(:,3) = data(:,36)./100;
VG(:,1) = data(:,37)./100;
VG(:,2) = data(:,38)./100;
VG(:,3)= data(:,39)./100;
[Plat_G, Plon_G, Palt_G] = ecef2lla(PG(:,1),PG(:,2),PG(:,3));
% constants
omega = 7.292115 * 10^-5; % angualar velocity of earth in rad/s
We_ie=[0;0;omega];
J_2 = 1.082629989051944 * 10^-3; % coefficient of second zonal harmonics of earths potential function
meu = 3986005 * 10^8; % earth's gravitational constant in m^3/s^2
R_e = 6378137.0000; % Radius of earth (semimajor axis) in m
dt = 1/50; % Sample period
I=eye(3);
O=zeros(3);

%%% INS Initialization %%%
%Intial position and velocity, our solution [cm deg]
PS(1,1)= PM(1,1);
PS(2,1)= PM(1,2);
PS(3,1)= PM(1,3);
VS(1,1)= VM(1,1);
VS(2,1)= VM(1,2);
VS(3,1)= VM(1,3);
% Transformation from ECEF to Geodetic form
[Plat_S(1,1), Plon_S(1,1), Palt_S(1,1)] = ecef2lla(PS(1,1),PS(2,1),PS(3,1));

CN_E(:,:,1) = [ -sin(Plat_S(1,1))*cos(Plon_S(1,1)) -sin(Plat_S(1,1))*sin(Plon_S(1,1)) cos(Plat_S(1,1));
    -sin(Plon_S(1,1))                     cos(Plon_S(1,1))               0;
    -cos(Plat_S(1,1))*cos(Plon_S(1,1))   -cos(Plat_S(1,1))*sin(Plon_S(1,1)) -sin(Plat_S(1,1))];
%Intial navigation angles our solution [cm deg]
yawS(1,1)=yawM(1,1);
pitchS(1,1)=pitchM(1,1);
rollS(1,1)=rollM(1,1);
% Transformation from Body frame to Navigation frame NED
CN_B (:,:,1)= [cos(pitchS(1,1))*cos(yawS(1,1)) sin(rollS(1,1))*sin(pitchS(1,1))*cos(yawS(1,1))-cos(rollS(1,1))*sin(yawS(1,1)) cos(rollS(1,1))*sin(pitchS(1,1))*cos(yawS(1,1))+sin(rollS(1,1))*sin(yawS(1,1));
    cos(pitchS(1,1))*sin(yawS(1,1)) sin(rollS(1,1))*sin(pitchS(1,1))*sin(yawS(1,1))+cos(rollS(1,1))*cos(yawS(1,1)) cos(rollS(1,1))*sin(pitchS(1,1))*sin(yawS(1,1))-sin(rollS(1,1))*cos(yawS(1,1));
    -sin(pitchS(1,1))                      sin(rollS(1,1))*cos(pitchS(1,1))          cos(rollS(1,1))*cos(pitchS(1,1))];
%Transformation from Body frame to ECEF
CE_B(:,:,1) =CN_E(:,:,1)\CN_B(:,:,1);
% Intialization of the quaterions
qoS=1/2*sqrt(1+CE_B(1,1,1)+CE_B(2,2,1)+CE_B(3,3,1));
q1S=(CE_B(3,2,1)-CE_B(2,3,1))/(4*qoS);
q2S=(CE_B(1,3,1)-CE_B(3,1,1))/(4*qoS);
q3S=(CE_B(2,1,1)-CE_B(1,2,1))/(4*qoS);
qS=[];
qS=[qS;[qoS q1S q2S q3S]];
qS(1,:)= qS(1,:)'/norm(qS(1,:)');
% Intialization of the gravitional vector
% Range frome the earth center to vehicle center
Rec = sqrt(PS(1,1)^2+PS(2,1)^2+PS(3,1)^2);
% Gravitational force vector in ECEF
Gec(:,1) = [ (-meu/Rec^2)*(PS(1,1)/Rec)*(1-((3/2)*J_2*(R_e/Rec)^2*((5*PS(3,1)^2/Rec^2)-1)));
    (-meu/Rec^2)*(PS(2,1)/Rec)*(1-((3/2)*J_2*(R_e/Rec)^2*((5*PS(3,1)^2/Rec^2)-1)));
    (-meu/Rec^2)*(PS(3,1)/Rec)*(1-((3/2)*J_2*(R_e/Rec)^2*((5*PS(3,1)^2/Rec^2)-3)))];
Wb_eb(:,1)= Wb(1,:)' - CE_B(:,:,1)'* We_ie ;

%%% Estimation Initialization %%%
iold=1;
%process uncertainty
Q_c=[1*I        O             O           O
    O          O          O           O
    O           O          0*I          O
    O           O             O         O];

q_kk=0.001.*ones(12,1);
%measurements covariance
R=[0.01        0       0         0       0       0
    0          0.01   0       0       0       0
    0           0       0.01  0       0       0
    0           0       0       0.3^2  0        0
    0           0       0       0       0.1^2  0
    0           0       0       0       0       0.1^2];
%Measurment matrix
H=[I O O O O
    O I O O O];

%Initial process covariance
Phat(:,:,1)=[1^2*I         O         O          O        O;
    O           0.01*I       O          O        O;
    O             O      0.00001*I       O        O;
    O             O         O        0.00001*I    O;
    O             O         O          O      0.000001*I];

ba_b(:,1) = [0;0;0];
bg_b(:,1) = [0;0;0];
noise2= 0.0005.*randn(6,1);
noise31=0.000001.*randn(1,1);
noise32=0.000001.*randn(1,1);
noise33= 0.001.*randn(1,1);
noise3=[noise31;noise32;noise33];
xhat(:,1)=0.00001*ones(15,1);
Corterms(:,1)=[-((2*VS(2,1))+(omega*PS(1,1)))*omega;((2*VS(1,1))-(omega*PS(2,1)))*omega;0];

for i = 1:(endp-1)
    
    noise2= 0.0005.*randn(6,1);
    noise31=0.000001.*randn(1,1);
    noise32=0.000001.*randn(1,1);
    noise33= 0.001.*randn(1,1);
    noise3=[noise31;noise32;noise33];
    %%%%%%%%%%INS Mechanization%%%%%%%%%%%
    %Position and Velocity
    PS(:,i+1) = PS(:,i)+ VS(:,i).*dt;
    VS(:,i+1) = VS(:,i)+( CE_B(:,:,i)*Fb(i,:)' - Corterms(:,i) +Gec(:,i))*dt;
    %Angular Motion
    qSS = qS(i,:)'+.5.*skewq(Wb_eb(:,i))*qS(i,:)'.*dt;
    
    % dq(:,i)=[cos((sqrt((Wb_eb(1,i)^2)+(Wb_eb(2,i)^2)+(Wb_eb(3,i)^2))*dt)/2)
    %     (sin((sqrt((Wb_eb(1,i)^2)+(Wb_eb(2,i)^2)+(Wb_eb(3,i)^2))*dt)/2))*(Wb_eb(1,i)/(sqrt((Wb_eb(1,i)^2)+(Wb_eb(2,i)^2)+(Wb_eb(3,i)^2))))
    %     (sin((sqrt((Wb_eb(1,i)^2)+(Wb_eb(2,i)^2)+(Wb_eb(3,i)^2))*dt)/2))*(Wb_eb(2)/(sqrt((Wb_eb(1,i)^2)+(Wb_eb(2,i)^2)+(Wb_eb(3,i)^2))))
    %     (sin((sqrt((Wb_eb(1,i)^2)+(Wb_eb(2,i)^2)+(Wb_eb(3,i)^2))*dt)/2))*(Wb_eb(3,i)/(sqrt((Wb_eb(1,i)^2)+(Wb_eb(2,i)^2)+(Wb_eb(3,i)^2))))];
    %
    %
    %     qSS=[dq(1,i) -dq(2,i) -dq(3,i) -dq(4,i)
    %        [dq(2,i);dq(3,i);dq(4,i)] dq(1,i)*(I-skew([dq(2,i);dq(3,i);dq(4,i)]))]*qS(i,:)';
    %
    qSS= qSS/norm(qSS);
    qS =[qS;qSS'];
    
    
    CE_B(:,:,i+1)=[ qS(i+1,1)^2+qS(i+1,2)^2-qS(i+1,3)^2-qS(i+1,4)^2 , 2*(qS(i+1,2)*qS(i+1,3)-qS(i+1,4)*qS(i+1,1)), 2*(qS(i+1,2)*qS(i+1,4)+qS(i+1,1)*qS(i+1,3));
        2*(qS(i+1,2)*qS(i+1,3)+qS(i+1,1)*qS(i+1,4)), qS(i+1,1)^2-qS(i+1,2)^2+qS(i+1,3)^2-qS(i+1,4)^2  , 2*(qS(i+1,3)*qS(i+1,4)-qS(i+1,1)*qS(i+1,2));
        2*(qS(i+1,2)*qS(i+1,4)-qS(i+1,1)*qS(i+1,3)) , 2*(qS(i+1,3)*qS(i+1,4)+qS(i+1,1)*qS(i+1,2)) , qS(i+1,1)^2-qS(i+1,2)^2-qS(i+1,3)^2+qS(i+1,4)^2];
    % Transformation from ECEF to Geodetic form
    [Plat_S(1,i+1), Plon_S(1,i+1), Palt_S(1,i+1)] = ecef2lla(PS(1,i+1),PS(2,i+1),PS(3,i+1));
    
    CN_E(:,:,i+1) = [ -sin(Plat_S(1,i+1))*cos(Plon_S(1,i+1)) -sin(Plat_S(1,i+1))*sin(Plon_S(1,i+1)) cos(Plat_S(1,i+1))
        -sin(Plon_S(1,i+1))                     cos(Plon_S(1,i+1))               0;
        -cos(Plat_S(1,i+1))*cos(Plon_S(1,i+1))   -cos(Plat_S(1,i+1))*sin(Plon_S(1,i+1)) -sin(Plat_S(1,i+1))];
    %Transformation from Body frame to Navigation fram
    CN_B(:,:,i+1)= CN_E(:,:,i+1)* CE_B(:,:,i+1) ;
    
    rollS(1,i+1)=atan2(CN_B(3,2,i+1),CN_B(3,3,i+1));
    pitchS(1,i+1)=asin(-CN_B(3,1,i+1));
    yawS(1,i+1)=atan2(CN_B(2,1,i+1),CN_B(1,1,i+1));
    
    Wb_eb(:,i+1)= Wb(i+1,:)'-CE_B(:,:,i+1)'*We_ie ;
    % Range frome the earth center to vehicle center
    Rec = sqrt(PS(1,i+1)^2+PS(2,i+1)^2+PS(3,i+1)^2);
    Gec(:,i+1) = [ (-meu/Rec^2)*(PS(1,i+1)/Rec)*(1-((3/2)*J_2*(R_e/Rec)^2*((5*PS(3,i+1)^2/Rec^2)-1)));
        (-meu/Rec^2)*(PS(2,i+1)/Rec)*(1-((3/2)*J_2*(R_e/Rec)^2*((5*PS(3,i+1)^2/Rec^2)-1)));
        (-meu/Rec^2)*(PS(3,i+1)/Rec)*(1-((3/2)*J_2*(R_e/Rec)^2*((5*PS(3,i+1)^2/Rec^2)-3)))];
    Corterms(:,i+1)=[-((2*VS(2,i+1))+(omega*PS(1,i+1)))*omega;((2*VS(1,i+1))-(omega*PS(2,i+1)))*omega;0];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if GPS_flag(i+1) == 1
        dtg=((i+1)-iold)*(1/50);
        iold=i+1;
        %gravity tensor
        gG=[ (-meu/Rec^3)*(1-3*(PS(1,i+1)/Rec)^2)       (3*meu/Rec^3)*(PS(1,i+1)*PS(2,i+1)/Rec^2)         3*(meu/Rec^3)*(PS(1,i+1)*PS(3,i+1)/Rec^2);
            3*meu/Rec^3*(PS(1,i+1)*PS(2,i+1)/Rec^2)      -meu/Rec^3*(1-3*(PS(2,i+1)/Rec)^2)              3*(meu/Rec^3)*(PS(2,i+1)*PS(3,i+1)/Rec^2);
            (3*meu/Rec^3)*(PS(1,i+1)*PS(3,i+1)/Rec^2)          3*meu/Rec^3*(PS(2,i+1)*PS(3,i+1)/Rec^2)      -meu/Rec^3*(1-3*(PS(3,i+1)/Rec)^2)];
        
        Somg=skew(We_ie);
        
        A=[     O             I              O              O                             O;
            gG-Somg^2    -2*Somg   -2*skew(CE_B(:,:,i+1)*Fb(i+1))     CE_B(:,:,i+1)       O;
            O             O           -Somg                             O      (0.5)*CE_B(:,:,i+1);
            O             O              O                               O             O;
            O             O              O                                O             O];
        F=eye(15)+A*dtg;
        
        B=[   O                     O       O    O;
            CE_B(:,:,i+1)         O         O    O;
            O       (0.5)*CE_B(:,:,i+1)   O    O;
            O            O                I    O;
            O            O               O    I];
        Gamma=B*dtg;
        z(:,i+1)=[PS(:,i+1)-PG(i+1,:)';VS(:,i+1)-VG(i+1,:)'];
        
        % %     %obtain the a priori estimate of the state
        %        xbar(:,i+1) = F * xhat(:,i);
        %
        % %     %obtain the a priori measurement
        %       zbar(:,i+1) = H * xbar(:,i+1);
        % %
        %    %obtain the innovation
        %     neu(:,i+1) = z(:,i+1) - zbar(:,i+1);
        %
        %Propagate the covariance
        Pbar(:,:,i+1) = F * Phat(:,:,i) * F' + Gamma * Q_c * Gamma';
        
        %The innovation covariance
        S(:,:,i+1) = H * Pbar(:,:,i+1) * H' + R;
        
        %The Filter gain
        %W(:,:,i+1) = Pbar(:,:,i+1) * H' * S(:,:,i+1)^(-1);
        %
        %delta = [0.1;0.1;0.1;0.1;0.1;0.1];
       zz = diag(abs(z(:,i+1))./(10*diag(R)));
       sat_z = max(-1, min(1, zz));
        W(:,:,i+1) = pinv(H)*sat_z;
       
        %Obtain the updated state
        xhat(:,i+1) =  W(:,:,i+1) * z(:,i+1);
        
        %Obtain the updated state covariance
        Phat(:,:,i+1) = Pbar(:,:,i+1) - W(:,:,i+1) *S(:,:,i+1) * W(:,:,i+1)';
        
        %estimate states update
        PS(:,i+1)=PS(:,i+1)-xhat(1:3,i+1);
        VS(:,i+1)=VS(:,i+1)-xhat(4:6,i+1);
        
        d_q=[1;xhat(7,i+1);xhat(8,i+1);xhat(9,i+1)];
        qSShat=[qSS(1) -qSS(2) -qSS(3) -qSS(4);
            qSS(2)  qSS(1)  qSS(4) -qSS(3);
            qSS(3) -qSS(4)  qSS(1)  qSS(2);
            qSS(4)  qSS(3) -qSS(2)  qSS(1)]* d_q;
        qSShat=qSShat/norm(qSShat);
        qS(i+1,:) = qSShat';
        ba_b(:,i+1) = ba_b(:,i)- xhat(10:12,i+1);
        bg_b(:,i+1) = bg_b(:,i) - xhat(13:15,i+1);
        %removing bias from gyro and accelerommeter  readingis
        Fb(i+1,:)= Fb(i+1,:) - ba_b(:,i+1)';
        Wb(i+1,:)= Wb(i+1,:) - bg_b(:,i+1)'; %%% kamal - noise3
        
        %%% re-calculate below valuse after apply ekf and we got the updated value
        CE_B(:,:,i+1)=[ qS(i+1,1)^2+qS(i+1,2)^2-qS(i+1,3)^2-qS(i+1,4)^2 , 2*(qS(i+1,2)*qS(i+1,3)-qS(i+1,4)*qS(i+1,1)), 2*(qS(i+1,2)*qS(i+1,4)+qS(i+1,1)*qS(i+1,3));
            2*(qS(i+1,2)*qS(i+1,3)+qS(i+1,1)*qS(i+1,4)), qS(i+1,1)^2-qS(i+1,2)^2+qS(i+1,3)^2-qS(i+1,4)^2  , 2*(qS(i+1,3)*qS(i+1,4)-qS(i+1,1)*qS(i+1,2));
            2*(qS(i+1,2)*qS(i+1,4)-qS(i+1,1)*qS(i+1,3)) , 2*(qS(i+1,3)*qS(i+1,4)+qS(i+1,1)*qS(i+1,2)) , qS(i+1,1)^2-qS(i+1,2)^2-qS(i+1,3)^2+qS(i+1,4)^2];
        % Transformation from ECEF to Geodetic form
        [Plat_S(1,i+1), Plon_S(1,i+1), Palt_S(1,i+1)] = ecef2lla(PS(1,i+1),PS(2,i+1),PS(3,i+1));
        
        CN_E(:,:,i+1) = [ -sin(Plat_S(1,i+1))*cos(Plon_S(1,i+1)) -sin(Plat_S(1,i+1))*sin(Plon_S(1,i+1)) cos(Plat_S(1,i+1))
            -sin(Plon_S(1,i+1))                     cos(Plon_S(1,i+1))               0;
            -cos(Plat_S(1,i+1))*cos(Plon_S(1,i+1))   -cos(Plat_S(1,i+1))*sin(Plon_S(1,i+1)) -sin(Plat_S(1,i+1))];
        %Transformation from Body frame to Navigation fram
        CN_B(:,:,i+1)= CN_E(:,:,i+1)* CE_B(:,:,i+1) ;
        
        rollS(1,i+1)=atan2(CN_B(3,2,i+1),CN_B(3,3,i+1));
        pitchS(1,i+1)=asin(-CN_B(3,1,i+1));
        yawS(1,i+1)=atan2(CN_B(2,1,i+1),CN_B(1,1,i+1));
        
        Wb_eb(:,i+1)= Wb(i+1,:)'-CE_B(:,:,i+1)'*We_ie ;
        Rec = sqrt(PS(1,i+1)^2+PS(2,i+1)^2+PS(3,i+1)^2);
        Gec(:,i+1) = -meu/Rec^2.*[ 1-3/2*J_2*(R_e/Rec)^2*(5*PS(3,i+1)^2/Rec^2-1);
            1-3/2*J_2*(R_e/Rec)^2*(5*PS(3,i+1)^2/Rec^2-1);
            1-3/2*J_2*(R_e/Rec)^2*(5*PS(3,i+1)^2/Rec^2-3)]./Rec.*PS(:,i+1);
        Corterms(:,i+1)=[-((2*VS(2,i+1))+(omega*PS(1,i+1)))*omega;((2*VS(1,i+1))-(omega*PS(2,i+1)))*omega;0];
        
    end
    
    if GPS_flag(i+1) ~= 1
        xhat(:,i+1)=xhat(:,i);
        Phat(:,:,i+1)=Phat(:,:,i);
        ba_b(:,i+1) = ba_b(:,i);
        bg_b(:,i+1) = bg_b(:,i);
        Fb(i+1,:)= Fb(i+1,:) - ba_b(:,i+1)' ;
        Wb(i+1,:)= Wb(i+1,:) - bg_b(:,i+1)'; %%% kamal - noise3
        
    end
    
    
end

screen=[0 50 1100 500];
legend1='MIDG';
legend2='solution';
time=(1:endp).*dt;
%Path
figure('Position',screen );
plot(Plon_S*180/pi,Plat_S*180/pi,'b','LineWidth',3);
hold on;
plot(Plon_M*180/pi,Plat_M*180/pi,'-.green','LineWidth',3);
hold on
plot(Plon_G*180/pi,Plat_G*180/pi,':r','LineWidth',3);
legend(legend2,legend1,'GPS');
xlabel('\Lambda [deg]');ylabel('\Phi [deg]');

% plot positions error histogram for car 1 (ECEF)
figure('Position',screen );
subplot(3,1,1);
hold on;
histfit(PS(1,:)'-PM(:,1));
xlabel('P_x^e error [m]');ylabel('relative frequency');
subplot(3,1,2);
hold on
histfit(PS(2,:)'-PM(:,2));
xlabel('P_y^e error [m]');ylabel('relative frequency');
subplot(3,1,3);
hold on
histfit(PS(3,:)'-PM(:,3));
xlabel('P_z^e error [m]');ylabel('relative frequency');

% plot velocities error histogram for car 1 (ECEF)
figure('Position',screen );
subplot(3,1,1);
hold on;
histfit(VS(1,:)'-VM(:,1));
xlabel('V_x^e error [m/s]');ylabel('relative frequency');
subplot(3,1,2);
hold on
histfit(VS(2,:)'-VM(:,2));
xlabel('V_y^e error [m/s]');ylabel('relative frequency');
subplot(3,1,3);
hold on
histfit(VS(3,:)'-VM(:,3));
xlabel('V_z^e error [m/s]');ylabel('relative frequency');

%plot accelerometer bias for car 1
figure('Position',screen );
subplot(3,1,1);
hold on;
plot(time,ba_b(1,:),'LineWidth',1.5);
xlabel('Time [s]');ylabel('x bias_ab');
subplot(3,1,2);
hold on
plot(time,ba_b(2,:),'LineWidth',1.5);
xlabel('Time [s]');ylabel('y bias_ab');
subplot(3,1,3);
hold on
plot(time,ba_b(3,:),'LineWidth',1.5);
xlabel('Time [s]');ylabel('z bias_ab');

%plot gyroscope bias for car 1
figure('Position',screen );
subplot(3,1,1);
hold on;
plot(time,bg_b(1,:),'LineWidth',1.5);
xlabel('Time [s]');ylabel('x bias_gb');
subplot(3,1,2);
hold on
plot(time,bg_b(2,:),'LineWidth',1.5);
xlabel('Time [s]');ylabel('y bias_gb');
subplot(3,1,3);
hold on
plot(time,bg_b(3,:),'LineWidth',1.5);
xlabel('Time [s]');ylabel('z bias_gb');
