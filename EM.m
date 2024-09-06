

% 主要的能源管理算法在function Action1=s_rule2(EnvConstants,LoggedSignals0,Action0,LoggedSignals1)中
% 提取上一步的信息LoggedSignals0和控制Action0以及当前信息LoggedSignals1得到下一步的Action1

clc;clear;close all;


% prepare calculation-------------------------------------------
% 载入数据，result_v_qP0为v的q_P=0的结果，result_v_qP1为v从0到50的q_P=max的结果，分别对应未满充和满充情况
% ？_qP0:满充情况，此时太阳能电池的功率输出会累积热效应，即q_P=0太阳能功率匹配负载
% ?_qP1:未满充情况，此时太阳能电池的功率输出不会累积热效应，即q_P=max
% 正向计算：V确定散热和功耗，t确定太阳能，E确定是否满充，得到P_batt
% [T,P_sol]=f(V,t)

load('Stratobus_v50.mat');

% 构建插值（从时间t和太阳能电池功率P_sol得到v）：取得不同速度v从0到50与时间t对应的锂电池功率和氦气变化率
data.t=linspace(0,24,(height(result_v(1).He)-1)/2+1);% 采用48h中后24h数据
data.v=0:30;
data.Hedm=zeros(length(data.t),length(data.v));
data.Psol=zeros(length(data.t),length(data.v));
for i=1:length(data.v)
    data.Hedm(:,i)=double(result_v(i).He.dm(length(data.t):height(result_v(1).He)));
    data.Psol(:,i)=double(result_v(i).NY.P_sol(length(data.t):height(result_v(1).He)));
    data.Preq(:,i)=double(result_v(i).NY.P_req(length(data.t):height(result_v(1).He)));
    data.T_He(:,i)=double(result_v(i).He.T(length(data.t):height(result_v(1).He)));
end

data.Hedm=data.Hedm'*3600;%kg/h
data.Psol=data.Psol'/1000;%kW
data.Preq=data.Preq'/1000;%kW
data.T_He=data.T_He';
%data.Hedm(:,1:8)=repmat(data.Hedm(:,9),1,8);%初始温度平衡慢
clearvars -except data Design

[data.T,data.V]=meshgrid(data.t,data.v);

% 构建插值：
Interp.tV2Psol= scatteredInterpolant(data.T(:),data.V(:),data.Psol(:), 'linear');%有外插
%Interp.tV2Psol.ExtrapolationMethod="none";%无外插
Interp.tV2Hedm= scatteredInterpolant(data.T(:),data.V(:),data.Hedm(:), 'linear');%有外插
Interp.tV2Hedm.ExtrapolationMethod="none";%无外插
Interp.tV2T_He= scatteredInterpolant(data.T(:),data.V(:),data.T_He(:), 'linear');%有外插
Interp.tV2T_He.ExtrapolationMethod="none";%无外插



Interp.tP2V=scatteredInterpolant(data.T(:),data.Preq(:)-data.Psol(:),data.V(:), 'linear');%有外插
Interp.tP2V.ExtrapolationMethod="none";%无外插
Interp.t2Pmin=griddedInterpolant(data.t,data.Preq(1,:)-data.Psol(1,:));% t到储能系统功率最小边界
Interp.t2Pmax=griddedInterpolant(data.t,data.Preq(end,:)-data.Psol(end,:));% t到储能系统功率最大边界
Interp.V2Preq=griddedInterpolant(data.V(:,1),data.Preq(:,1));

% plot wind vs t
figure(Name='wind data');
rng(0);wind=smooth(28*rand([73,1]));wind=wind/mean(wind)*15;mean(wind)
plot(wind);xlabel('t (h)');ylabel('wind velocity (m/s)');set(gca,'Fontname', 'Times New Roman','FontSize',12);
Interp.wind=griddedInterpolant(0:72,wind);%ave=15m/s,%m/s wind velocity
% plot Preq vs v
figure(Name='power requirement vs velocity');
plot(data.V(:,1),data.Preq(:,1));xlabel('velocity (m/s)');ylabel('power requirement (kW)');
annotation('textarrow',[0.171428571428571 0.133333333333333],[0.278571428571429 0.196031746031746],'String',{'fixed power consumption'},'Fontname','Times New Roman','FontSize',12);
set(gca,'Fontname', 'Times New Roman','FontSize',12);



% %% 绘图确定采样点都在插值面上
% t_rand=rand(100000,1)*12;
% p_rand=rand(100000,1)*300-200;
% v_rand=rand(100000,1)*30;
% v_interp=Interp.tP2V(t_rand,p_rand);
% v_interp(p_rand>Interp.t2Pmax(t_rand) | p_rand<Interp.t2Pmin(t_rand))=nan;
% Hedm_interp=Interp.tV2Hedm(t_rand,v_interp);
% figure(1);hold on;surf(data.T,data.V,data.Preq-data.Psol,"EdgeColor","none");scatter3(t_rand,v_rand,-Interp.tV2Psol(t_rand,v_rand)+Interp.V2Preq(v_rand));xlabel('t');ylabel('v')
% figure(2);hold on;surf(data.T,data.V,data.Hedm);surf(data.T,data.V,data.Hedm,"EdgeColor","none");xlabel('t (h)');ylabel('v (m/s)');zlabel('$\dot m_{\mathrm{He}}$','Interpreter','latex');set(gca,'Fontname', 'Times New Roman')
% %scatter3(t_rand,v_interp,Hedm_interp);
% figure(3);hold on;surf(data.T,data.Preq-data.Psol,data.V);xlabel('t (h)');ylabel('P (kW)');zlabel('v (m/s)');set(gca,'Fontname', 'Times New Roman');scatter3(t_rand,p_rand,v_interp);xlabel('$\dot t_2 $','Interpreter','latex')



% RFC数据：FC data from https://doi.org/10.1155/2019/7860214;elec from doi:10.3390/en13030612
load('EFC_eff.mat');

dEdt_RFC=(P_RFC<0).*-P_RFC.*eff_RFC+(P_RFC>0).*-P_RFC./eff_RFC;
dEdt_RFC(P_RFC==0)=-3;%考虑热管理1kW不用dEdt_RFC(P_RFC==0)=0;
Interp.P2dEdt_RFC=griddedInterpolant(P_RFC,dEdt_RFC);
Interp.dEdt2P_RFC=griddedInterpolant(flipud(dEdt_RFC),flipud(P_RFC));
% plot RFC
figure(Name='RFC data');
yyaxis left;plot(P_RFC,eff_RFC,".-");ylabel('charge/dischage efficiency');ylim([0 1]);
yyaxis right;plot(P_RFC,dEdt_RFC,".-");ylabel('Energy change rate (kWh/h)');
xlabel('Power (kW)');
set(gca,'Fontname', 'Times New Roman','FontSize',12);


% EnvConstants:72h 可放空，软驻留约束

% Design
Design.E_bat_max=150;%kWh,8kW14h
Design.E_RFC_max=900;%kWh
Design.P_bat_max=100;%100kW
Design.P_RFC_max=100;%100kW
Design.eff_dischrg_bat=0.97;
Design.eff_chrg_bat=0.97;
% Design.eff_dischrg_RFC=0.6;
% Design.eff_chrg_RFC=0.82;
Design.v_max=30;%m/s

% mission EnvConstants：
EnvConstants.H=20000;%m


%
EnvConstants.Ts=1;%h
EnvConstants.t=(0:EnvConstants.Ts:72)';%h, 时间离散
EnvConstants.t_24h=EnvConstants.t-24*floor(EnvConstants.t/24);
[EnvConstants.t_N,~]=size(EnvConstants.t);

EnvConstants.PenaltyForFalling=-100;%失败后每一步都惩罚该值
EnvConstants.r_max=50;%km，区域驻留半径

% J_temp(isnan(J_temp))=inf;%把nan转成inf

% discretization variables：
EnvConstants.x_N=100;% 离散个数，也许应当使步长小于每个时间步长的变量
EnvConstants.x1=linspace(0,Design.E_bat_max,EnvConstants.x_N);
EnvConstants.x2=linspace(0,Design.E_RFC_max,EnvConstants.x_N);

% initial condition：
init.E_bat=0.9*Design.E_bat_max;%kWh
init.E_RFC=0.9*Design.E_RFC_max;%kWh
init.R_back=0;%km
init.m_He=Design.m_He_initial;%kg

EnvConstants.init=init;%初始能量
EnvConstants.Interp=Interp;%获取插值
EnvConstants.Design=Design;%设计参数
clearvars -except EnvConstants

% J_temp=J_temp_bat+J_temp_RFC+J_temp_back+J_temp_leak,量级（0.1C约100kW，100kW，54km/h，1kg/d）都相近于1
% EnvConstants.J_temp_cal=@(t,x_cur,x_next) max(x_cur(:,1)-x_next(:,1),0)/200+max(x_cur(:,2)-x_next(:,2),0)/200+(max(0,x_cur(:,3)-EnvConstants.r_max)+max(0,x_next(:,3)-EnvConstants.r_max))/108+max(0,x_cur(:,4)-x_next(:,4))*24;

%% demo：
tic
u_demo=set2u(EnvConstants.t,15,15,8,17,1/7,1/7,9,18);
[x_demo,J_demo,P_demo,stat_demo(1)]=OCP_plot(EnvConstants,u_demo,'demo');
toc

%% 在线算法模拟：根据s_rule和s_greedy：
tic
Ts=0.01;
t=(EnvConstants.t(1):Ts:EnvConstants.t(end))';
t_24h=t-24*floor(t/24);
EnvConstants.t=t;
%strategy=@ s_greedy;
strategy=@ s_rule2;
[x_plot,J_plot,P_plot,u_act]=dyn_strategy(EnvConstants,strategy);
J_plot(end)
u_plot=u_act;

Design=EnvConstants.Design;

if t(end)==72
    xlim_max=100;
else
    xlim_max=36;
end
figure();
subplot(3,1,1);xlim([0 xlim_max])
yyaxis left;plot(t,u_plot(:,1));ylabel('\itv \rm(m·s^-^1)');ylim([0 35])
yyaxis right;hold on;plot(t,interp1(linspace(t(1),t(end),size(u_plot,1)),u_plot(:,2),t,"previous"),'-');plot(t,u_act(:,2),":",LineWidth=2);ylabel('\itr_b_a_t');ylim([0 1])
legend("\itv","\itr_b_a_t","actural \itr_b_a_t","Location","east")
set(gca,'Fontname', 'Times New Roman','FontWeight','norm','FontSize',10)
hold on;

subplot(3,1,2);xlim([0 xlim_max])
yyaxis left;hold on;plot(t,x_plot(:,1)./Design.E_bat_max,"-");plot(t,x_plot(:,2)./Design.E_RFC_max,':',LineWidth=2);ylabel('\itSOC');ylim([0 1])
yyaxis right;plot(t,x_plot(:,3),"-");ylabel('\itR_b_a_c_k \rm(km)');ylim([-100 600])
legend("\itSOC_b_a_t","\itSOC_R_F_C","\itR_b_a_c_k","Location","east")
set(gca,'Fontname', 'Times New Roman','FontWeight','norm','FontSize',10)
hold on;

subplot(3,1,3);xlabel('\itt \rm(h)');xlim([0 xlim_max])
yyaxis left;hold on;plot(t,P_plot(:,1),"-");plot(t,P_plot(:,2),':',LineWidth=2);ylabel('\itP \rm(kW)');ylim([0 200]);
yyaxis right;plot(t,P_plot(:,3),"-");plot(t,P_plot(:,4),':',LineWidth=2);ylabel('\itP \rm(kW)');ylim([-150 50]);
legend('\itP_s_o_l','\itP_r_e_q',"\itP_b_a_t","\itP_R_F_C","Location","east")
set(gca,'Fontname', 'Times New Roman','FontWeight','norm','FontSize',10)
hold on;

stat.final_battery_energy=x_plot(end,1);
stat.final_RFC_energy=x_plot(end,2);
stat.final_withdraw=x_plot(end,3);
stat.gas_loss=x_plot(end,4)-x_plot(1,4);
stat.average_volecity=mean(u_act(:,1));
stat.accumulated_distance_beyond_residency=Ts*sum((x_plot(:,3)>EnvConstants.r_max).*(x_plot(:,3)-EnvConstants.r_max));%km·h
stat.total_battery_consumption=sum(max(diff(x_plot(:,1)),0));
stat.total_RFC_consumption=sum(max(diff(x_plot(:,2)),0));
stat.curtailment_rate=100*sum(EnvConstants.Interp.tV2Psol(t_24h,u_act(:,1))-P_plot(:,1))./sum(EnvConstants.Interp.tV2Psol(t_24h,u_act(:,1)));%
stat.total_cost=J_plot(end);


function Action1=s_greedy(EnvConstants,LoggedSignals0,Action0,LoggedSignals1)% 规则：从之前的状态变化到控制



EnvConstants.J_temp_cal=@(t,x_cur,x_next) max(x_cur(:,1)-x_next(:,1),0)/200+max(x_cur(:,2)-x_next(:,2),0)/200+(max(0,x_cur(:,3)-EnvConstants.r_max)+max(0,x_next(:,3)-EnvConstants.r_max))/108+max(0,x_cur(:,4)-x_next(:,4))*24;

v_dis=linspace(0,EnvConstants.Design.v_max,100);
r_bat_dis=linspace(0,1,100);
[v_dis,r_bat_dis]=meshgrid(v_dis,r_bat_dis);
v_dis=v_dis(:);
r_bat_dis=r_bat_dis(:);
Action=[v_dis,r_bat_dis];
LoggedSignals1.t=repmat(LoggedSignals1.t,length(v_dis),1);
LoggedSignals1.State=repmat(LoggedSignals1.State,length(v_dis),1);
LoggedSignals1.Action_actual=repmat(LoggedSignals1.Action_actual,length(v_dis),1);
LoggedSignals1.P=repmat(LoggedSignals1.P,length(v_dis),1);

[~,Reward,~,~] = EM_step_pre(Action,LoggedSignals1,EnvConstants);%使用这个，因为预先不知道风场，基于平均风15m/s计算



[~,i_max]=max(Reward);
Action1(1)=v_dis(i_max);
Action1(2)=r_bat_dis(i_max);
end


function Action1=s_rule(EnvConstants,LoggedSignals0,Action0,LoggedSignals1)% 规则：从之前的状态变化到控制

% Unpack signals:
v0=Action0(1);
r_bat0=Action0(2);
t0=LoggedSignals0.t;
E_bat0=LoggedSignals0.State(1);
E_RFC0=LoggedSignals0.State(2);
R_back0=LoggedSignals0.State(3);

t1=LoggedSignals1.t;
t_24h1=(t1-24*floor(t1/24));
P_sol1=LoggedSignals1.P(1);
E_bat1=LoggedSignals1.State(1);
E_RFC1=LoggedSignals1.State(2);
R_back1=LoggedSignals1.State(3);



% 变化步长：
dv=15/80;dr_bat=0.01;

% v：白天P_sol>0,充电，增加v直至上限，否则降低v直至15；夜间控制速度到明天九点退飞
% r_bat:控制dE_RFC/dt在49或-24的高效率区间，假如E_bat正常
if P_sol1>0 %白天
    if E_bat1+E_RFC1>=E_bat0+E_RFC0 %&& (15-t_24h1)*(E_bat1+E_RFC1-E_bat0-E_RFC0)/(t1-t0)>EnvConstants.Design.E_RFC_max+EnvConstants.Design.E_bat_max-E_RFC1-E_bat1%能源增加,应该能充满
        v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
    else %能源减小
        v1=max(15,v0-dv); %减小v直至下限
    end      
else %夜间
    if t_24h1<9
        t_night_remain=9-t_24h1;
    else
        t_night_remain=24+9-t_24h1;
    end
    if t_night_remain*(R_back1-R_back0)/(t1-t0)>EnvConstants.r_max-R_back1 && t_night_remain*(E_bat1+E_RFC1-E_bat0-E_RFC0)/(t1-t0)<E_RFC1+E_bat1%要飞出，有余电
        v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
    else %还有余量
        v1=max(0,v0-dv); %减小v直至下限
    end
end
 
if E_RFC1+E_bat1>E_RFC0+E_bat0 %充电
    if (E_RFC1-E_RFC0)/(t1-t0)<49 %功率小
        r_bat1=r_bat0-dr_bat; 
    else
        r_bat1=r_bat0+dr_bat;
    end
else %放电
   if (E_RFC1-E_RFC0)/(t1-t0)>-24 %功率小
        r_bat1=r_bat0-dr_bat;
    else
        r_bat1=r_bat0+dr_bat;
    end
end
Action1(1)=v1;
Action1(2)=r_bat1;
end

function Action1=s_rule2(EnvConstants,LoggedSignals0,Action0,LoggedSignals1)% 规则：从之前的状态变化到控制


% Unpack signals:上一步的控制和状态+这步状态
v0=Action0(1);
r_bat0=Action0(2);
t0=LoggedSignals0.t;
t_24h0=(t0-24*floor(t0/24));
E_bat0=LoggedSignals0.State(1);
E_RFC0=LoggedSignals0.State(2);
R_back0=LoggedSignals0.State(3);


t1=LoggedSignals1.t;
t_24h1=(t1-24*floor(t1/24));
P_sol1=LoggedSignals1.P(1);
E_bat1=LoggedSignals1.State(1);
E_RFC1=LoggedSignals1.State(2);
R_back1=LoggedSignals1.State(3);
P_req1=LoggedSignals1.P(2);
P_bat1=LoggedSignals1.P(3);
P_RFC1=LoggedSignals1.P(4);

%% estimation and forecasting


T_He0=EnvConstants.Interp.tV2T_He(t_24h1,v0);

t_sunrise=fsolve(@(t) EnvConstants.Interp.tV2Psol(t,v0)-0.5,6.5,optimset('Display','off'));%接近6.5
t_chrg=fsolve(@(t) EnvConstants.Interp.tV2Psol(t,v0)-EnvConstants.Interp.V2Preq(v0),9,optimset('Display','off'));
t_full=15;
t_dischrg=fsolve(@(t) EnvConstants.Interp.tV2Psol(t,v0)-EnvConstants.Interp.V2Preq(v0),15,optimset('Display','off'));
t_sunset=fsolve(@(t) EnvConstants.Interp.tV2Psol(t,v0)-0.5,18.5,optimset('Display','off'));%接近18.5
E_sol_remain=integral(@(t) EnvConstants.Interp.tV2Psol(t,v0*ones(size(t))),t_24h1,24);
if E_bat1>E_bat0 %充电
    E_bat_remain=E_bat1*EnvConstants.Design.eff_dischrg_bat*0.95;% 模型计算折算到供电系统的放电量
    E_bat_tofull=P_bat1*(EnvConstants.Design.E_bat_max-E_bat1)/(E_bat0-E_bat1)*(t1-t0);%测量满充值
elseif E_bat1<E_bat0 %放电
    E_bat_remain=P_bat1*E_bat1/(E_bat0-E_bat1)*(t1-t0)*0.95;% 测量余电
    E_bat_tofull=(EnvConstants.Design.E_bat_max-E_bat1)/EnvConstants.Design.eff_chrg_bat;% 模型计算，折算到供电系统的到满充电量
else % 不动
    E_bat_remain=E_bat1*EnvConstants.Design.eff_dischrg_bat*0.95;% 模型计算折算到供电系统的放电量
    E_bat_tofull=(EnvConstants.Design.E_bat_max-E_bat1)/EnvConstants.Design.eff_chrg_bat;% 模型计算，折算到供电系统的到满充电量
end
if E_RFC1>E_RFC0 %充电
    E_RFC_remain=-0.1*EnvConstants.Design.E_RFC_max*E_RFC1/EnvConstants.Interp.P2dEdt_RFC(0.1*EnvConstants.Design.E_RFC_max)*0.95;
    E_RFC_tofull=P_RFC1*(EnvConstants.Design.E_RFC_max-E_RFC1)/(E_RFC0-E_RFC1)*(t1-t0);
elseif E_RFC1<E_RFC0 %放电
    E_RFC_remain=P_RFC1*E_RFC1/(E_RFC0-E_RFC1)*(t1-t0)*0.95;
    E_RFC_tofull=0.1*EnvConstants.Design.E_RFC_max*(EnvConstants.Design.E_RFC_max-E_RFC1)/EnvConstants.Interp.P2dEdt_RFC(-0.1*EnvConstants.Design.E_RFC_max);
else %不动
    E_RFC_remain=-0.1*EnvConstants.Design.E_RFC_max*E_RFC1/EnvConstants.Interp.P2dEdt_RFC(0.1*EnvConstants.Design.E_RFC_max)*0.95;
    E_RFC_tofull=0.1*EnvConstants.Design.E_RFC_max*(EnvConstants.Design.E_RFC_max-E_RFC1)/EnvConstants.Interp.P2dEdt_RFC(-0.1*EnvConstants.Design.E_RFC_max);
end
E_ESS_remain=E_bat_remain+E_RFC_remain;
E_ESS_tofull=E_bat_tofull+E_RFC_tofull;
if isnan(E_ESS_remain) || isnan(E_ESS_tofull)
    disp('E_ESS nan occured');
    return
end



T_thold_h=230;
T_thold_l=218;
dEdt_RFC_dischrg=-0.0267*900;
dEdt_RFC_chrg=0.0544*900;


% 变化步长：
dv=15/80;dr_bat=0.01;

% v：白天P_sol>0,充电，增加v直至上限，否则降低v直至15；夜间控制速度到明天九点退飞
% r_bat:控制dE_RFC/dt在49或-24的高效率区间，假如E_bat正常
if P_sol1>0.5 %白天
    if E_ESS_remain>(t_chrg-t_24h1)*P_req1 %白天不放空
        if R_back1<EnvConstants.r_max
            if (24+t_chrg-t_24h1)*(R_back1-R_back0)/(t1-t0)<EnvConstants.r_max-R_back1 %到第二天充电不飞出
                if t_24h1<t_full
                    if E_sol_remain-(t_dischrg-t_24h1)*P_req1>E_ESS_tofull %能充满
                        if EnvConstants.Interp.tV2T_He(t_24h1,v0+dv)>T_He0 && T_He0>T_thold_l
                            v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
                        elseif T_He0>T_thold_h
                            v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
                        else
                            v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
                        end
                    else
                        v1=max(0,v0-dv); %减小v直至0
                    end
                elseif t_24h1<t_dischrg% t_full到t_dischrg
                    if E_sol_remain-(t_dischrg-t_24h1)*P_req1>E_ESS_tofull %能充满
                        v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
                    else
                        v1=max(0,v0-dv); %减小v直至0
                    end
                else %t_dischrg 到夜间
                    if E_ESS_remain>EnvConstants.Interp.V2Preq(15)*(t_chrg+24-t_dischrg) %夜间余量不足
                        if EnvConstants.Interp.tV2T_He(t_24h1,v0-dv)<T_He0 && T_He0<T_thold_h %降温
                            v1=max(0,v0-dv); %减小v直至0
                        elseif T_He0<T_thold_l
                            v1=max(0,v0-dv); %减小v直至0
                        else
                            v1=v0;
                        end
                    else
                        v1=max(0,v0-dv); %减小v直至0
                    end
                end
            else %要飞出
                if E_sol_remain-(t_dischrg-t_24h1)*P_req1>E_ESS_tofull %能充满
                    v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
                elseif E_sol_remain+E_bat1+E_RFC1> (24+t_chrg-t_24h1)*P_req1 && E_bat1+E_RFC1>(24+t_chrg-t_24h1)*P_req1 %够第二天
                    v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
                else
                    v1=max(0,v0-dv); %减小v直至0
                end
            end
        else %已经飞出
            if (24+t_chrg-t_24h1)*(R_back1-R_back0)/(t1-t0)<EnvConstants.r_max-R_back1 %能回
                if E_sol_remain-(t_dischrg-t_24h1)*P_req1>E_ESS_tofull  %能充满
                    v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
                elseif E_sol_remain+E_bat1+E_RFC1> (24+t_chrg-t_24h1)*P_req1 && E_bat1+E_RFC1>(24+t_chrg-t_24h1)*P_req1 %够第二天
                    v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
                else
                    v1=max(0,v0-dv); %减小v直至0
                end
            else %会不来
                if E_sol_remain-(t_dischrg-t_24h1)*P_req1>E_ESS_tofull %能充满
                    v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
                elseif E_sol_remain+E_bat1+E_RFC1> (24+t_chrg-t_24h1)*P_req1 && E_bat1+E_RFC1>(24+t_chrg-t_24h1)*P_req1 %够第二天
                    v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
                else
                    v1=max(0,v0-dv); %减小v直至0
                end
            end
        end
    else %白天要放空
        v1=max(0,v0-dv); %减小v直至0
    end


else %夜间
    if E_ESS_remain>((t_chrg-t_24h1<0)*(24+t_chrg-t_24h1)+(t_chrg-t_24h1>=0)*(t_chrg-t_24h1))*P_req1 %有余电
        if ((t_chrg-t_24h1<0)*(24+t_chrg-t_24h1)+(t_chrg-t_24h1>=0)*(t_chrg-t_24h1))*(R_back1-R_back0)/(t1-t0)<EnvConstants.r_max-R_back1 %有余电不飞出
            v1=max(0,v0-dv); %减小v直至下限
        else %有余电飞出
            v1=min(EnvConstants.Design.v_max,v0+dv); %增加v直至上限
        end
    else %没余电
        v1=max(0,v0-dv); %减小v直至0
    end
end


if E_RFC1+E_bat1>E_RFC0+E_bat0 %充电
    if (E_RFC1-E_RFC0)/(t1-t0)<dEdt_RFC_chrg %功率小
        r_bat1=r_bat0-dr_bat; 
    else
        r_bat1=r_bat0+dr_bat;
    end
else %放电
   if (E_RFC1-E_RFC0)/(t1-t0)>dEdt_RFC_dischrg %功率小
        r_bat1=r_bat0-dr_bat;
    else
        r_bat1=r_bat0+dr_bat;
    end
end
Action1(1)=v1;
Action1(2)=r_bat1;
end


function [x,J,P,u_act]=dyn_strategy(EnvConstants,strategy)  %根据输入u开展仿真
% 输入：EnvConstants.t，EnvConstants.init, EnvConstants其它，控制变量序列u（速度、r_bat）
% 输出：x（锂电池电量E_bat、RFC电量E_RFC、退飞量R_back、漏气量m_He）和J


% unpack,统一EnvConstants.t和控制变量序列u:

t=EnvConstants.t(:);%时间
t_24h=t-24*floor(t/24);
t_N=length(t);
EnvConstants.Ts=t(2)-t(1);



% initialize:
x=zeros(t_N,4);% x四列分别为E_bat E_RFC R_back m_He
J=zeros(t_N,1);% cost=-Reward
P=zeros(t_N,4);% 分别为P_sol,P_req,P_bat,P_RFC
u_act=zeros(t_N,2);% 分别为v,r_bat_actual
v=zeros(t_N,1);
r_bat=zeros(t_N,1);

[InitialObservation, LoggedSignals_prev] = EM_reset(EnvConstants);
x(1,:)=InitialObservation;
Action_prev=[15 0.1429];% v r_bat
[~,~,~,LoggedSignals_cur] = EM_step(Action_prev,LoggedSignals_prev,EnvConstants);

for i =2:t_N
    
    %v(i-1)=predict(s_v,[t_24h(i-1),EnvConstants.Interp.tV2Psol(t_24h(i-1),v(i-1)),x(i-1,1:3)]);
    %r_bat(i-1)=predict(s_r_bat,[t_24h(i-1),EnvConstants.Interp.tV2Psol(t_24h(i-1),v(i-1)),x(i-1,1:3)]);

    Action_cur=strategy(EnvConstants,LoggedSignals_prev,Action_prev,LoggedSignals_cur);
    v(i-1)=Action_cur(1);
    r_bat(i-1)=Action_cur(2);

    v(i-1)=max(0,min(30,v(i-1))); %限幅0~30
    r_bat(i-1)=max(0,min(1,r_bat(i-1))); %限幅0~1

    Action(1)=v(i-1);
    Action(2)=r_bat(i-1);
    LoggedSignals.t=t(i-1);
    LoggedSignals.State=x(i-1,:);
    LoggedSignals_prev=LoggedSignals;% log
    Action_prev=Action;% log
    [NextObs,Reward,IsDone,LoggedSignals] = EM_step(Action,LoggedSignals,EnvConstants);
    LoggedSignals_cur=LoggedSignals;% log

    %代价：
    if IsDone
        J(i:end)=-EnvConstants.PenaltyForFalling*(t_N-i+1); %不允许的代价
        break;
    else
        J(i)=J(i-1)-Reward;
    end

    x(i,:)=NextObs;
    P(i-1,:)=LoggedSignals.P;%[P_sol,P_req,P_bat,P_RFC];
    u_act(i-1,:)=LoggedSignals.Action_actual;%[v,r_bat_actual];
    toc
end

end






function u=set2u(t,V_night,V_day,t_inc_V,t_dec_V,r_day,r_night,t_inc_r,t_dec_r) % 根据设定、夜间速度、昼间速度、加速时间、减速时间计算总代价

t_24h=t-24*floor(t/24);
u(:,1)=ones(size(t_24h)).*V_night+(t_24h>t_inc_V & t_24h<t_dec_V).*(V_day-V_night);%V
u(:,2)=ones(size(t_24h)).*r_night+(t_24h>t_inc_r & t_24h<t_dec_r).*(r_day-r_night);%r_bat

end

function J_f=xJ2J(EnvConstants,u) % 根据设定输入计算总代价
[~,J]=dyn(EnvConstants,u);
J_f=J(end);
end

function [x,J,P,u_act]=dyn(EnvConstants,u)  %根据输入u开展仿真
% 输入：EnvConstants.t，EnvConstants.init, EnvConstants其它，控制变量序列u（速度、r_bat）
% 输出：x（锂电池电量E_bat、RFC电量E_RFC、退飞量R_back、漏气量m_He）和J


% unpack,统一EnvConstants.t和控制变量序列u:

t=EnvConstants.t(:);%时间
t_N=length(t);
EnvConstants.Ts=t(2)-t(1);

v=interp1(linspace(t(1),t(end),size(u,1)),u(:,1),t,"previous");
r_bat=interp1(linspace(t(1),t(end),size(u,1)),u(:,2),t,"previous");

% initialize:
x=zeros(t_N,4);% x四列分别为E_bat E_RFC R_back m_He
J=zeros(t_N,1);% cost=-Reward
P=zeros(t_N,4);% 分别为P_sol,P_req,P_bat,P_RFC
u_act=zeros(t_N,2);% 分别为v,r_bat_actual
x(1,:)=[EnvConstants.init.E_bat EnvConstants.init.E_RFC 0 EnvConstants.init.m_He];
isFall=false;


for i =2:length(v)
    Action(1)=v(i-1);
    Action(2)=r_bat(i-1);
    LoggedSignals.t=t(i-1);
    LoggedSignals.State=x(i-1,:);
    [NextObs,Reward,IsDone,LoggedSignals] = EM_step(Action,LoggedSignals,EnvConstants);

    %代价：
    if IsDone
        J(i:end)=-EnvConstants.PenaltyForFalling*(t_N-i+1); %不允许的代价
        break;
    else
        J(i)=J(i-1)-Reward;
    end

    x(i,:)=NextObs;
    P(i-1,:)=LoggedSignals.P;%[P_sol,P_req,P_bat,P_RFC];
    u_act(i-1,:)=LoggedSignals.Action_actual;%[v,r_bat_actual];
end

% tic
% for i=1:100
% ode_fun=@(t,x) dynamics(t,x,V_sim_gre(:),EnvConstants);
% sol = ode45(ode_fun, EnvConstants.t_sim, [630 0]);
% end
% toc
% for i=1:100
% out_gre=sim_by_V(EnvConstants,V_sim_gre);%给定
% end
% toc
% plot(sol,x,sol.y(1,:));hold on;plot(t',out_gre.E);plot(t',out_gre.R_back);plot(sol,x,sol.y(2,:))
end

function [u_temp,dx1dt,dx2dt,isAllow]=state2ctrl(t_24h,x1_cur,x2_cur,x1_next,x2_next,EnvConstants) %根据状态变化计算控制变量和单步代价,向量格式
% 输入：时间t（h，值)和每步变量dt、前后状态x_cur和x_next（kWh，阵到值、值到阵、阵到阵）
% 输出：J_temp(m_cur,n_cur,m_next,n_next)
Interp=EnvConstants.Interp;
E_bat_cur=x1_cur;%第1列
E_RFC_cur=x2_cur;%第2列
E_bat_next=x1_next;%第1列
E_RFC_next=x2_next;%第2列


P_bat=-(E_bat_next-E_bat_cur>=0).*(E_bat_next-E_bat_cur)/EnvConstants.Ts/EnvConstants.Design.eff_chrg_bat+...
    -(E_bat_next-E_bat_cur<0).*(E_bat_next-E_bat_cur)/EnvConstants.Ts*EnvConstants.Design.eff_dischrg_bat;% kW,放电+
P_RFC=Interp.dEdt2P_RFC((E_RFC_next-E_RFC_cur)/EnvConstants.Ts);% kW,放电+
r_bat=P_bat./(P_bat+P_RFC);%这里0/0=nan
r_bat(isnan(r_bat))=0;%把0/0的nan变为0
% isAllow为允许的序号，通过先计算isAllow再插值减小计算量
isAllow=P_bat>-EnvConstants.Design.P_bat_max & P_bat<EnvConstants.Design.P_bat_max & P_RFC>-EnvConstants.Design.P_RFC_max & P_RFC<EnvConstants.Design.P_RFC_max;%不允许功率超过1C
isAllow=isAllow & (P_bat.*P_RFC>=0);%不允许电源互充
isAllow=isAllow & ((P_bat+P_RFC)<EnvConstants.Interp.t2Pmax(t_24h)) & ((P_bat+P_RFC)>EnvConstants.Interp.t2Pmin(t_24h));%不外插,满充
t_allow=t_24h*isAllow;
% v_allow为允许的v，比完整v小很多：约1.2%
v_allow=EnvConstants.Interp.tP2V(t_allow(isAllow),(P_bat(isAllow)+P_RFC(isAllow)));%飞行速度，m/s
v_allow(v_allow>EnvConstants.Design.v_max)=EnvConstants.Design.v_max;
v=nan(size(isAllow));
v(isAllow)=v_allow;
v(E_bat_next>=0.999*EnvConstants.Design.E_bat_max & E_RFC_next>=0.999*EnvConstants.Design.E_RFC_max & P_bat+P_RFC>EnvConstants.Interp.t2Pmax(t_24h))=EnvConstants.Design.v_max;% 满充是充电功率可能已被削减，需要允许该情况
isAllow=~isnan(v);

u_temp(:,1)=v;
r_bat(~isAllow)=nan;%nan一致
u_temp(:,2)=r_bat;

dx1dt=x1_next-x1_cur;
dx2dt=x2_next-x2_cur;
dx1dt(~isAllow)=nan;%nan一致
dx2dt(~isAllow)=nan;%nan一致
end




function [x_plot,J_plot,P_plot,stat]=OCP_plot(EnvConstants,u_plot,name)
Ts=0.01;
t=(EnvConstants.t(1):Ts:EnvConstants.t(end))';
t_24h=t-24*floor(t/24);
EnvConstants.t=t;
[x_plot,J_plot,P_plot,u_act]=dyn(EnvConstants,u_plot);
v_plot=interp1(linspace(t(1),t(end),size(u_plot,1)),u_plot(:,1),t,"previous");
u_act(end,:)=u_act(end-1,:);

Design=EnvConstants.Design;
disp(['cost of' name 'is' string(J_plot(end))]);

if t(end)==72
    xlim_max=100;
else
    xlim_max=36;
end
figure('Name',name);
subplot(3,1,1);xlim([0 xlim_max])
yyaxis left;plot(t,v_plot,"-");ylabel('\itv \rm(m·s^-^1)');ylim([0 35])
yyaxis right;hold on;plot(t,interp1(linspace(t(1),t(end),size(u_plot,1)),u_plot(:,2),t,"previous"),'-');plot(t,u_act(:,2),":",LineWidth=2);ylabel('\itr_b_a_t');ylim([0 1])
legend("\itv","\itr_b_a_t","actural \itr_b_a_t","Location","east")
set(gca,'Fontname', 'Times New Roman','FontWeight','norm','FontSize',10)
hold on;

subplot(3,1,2);xlim([0 xlim_max])
yyaxis left;hold on;plot(t,x_plot(:,1)./Design.E_bat_max,"-");plot(t,x_plot(:,2)./Design.E_RFC_max,':',LineWidth=2);ylabel('\itSOC');ylim([0 1])
yyaxis right;plot(t,x_plot(:,3),"-");ylabel('\itR_b_a_c_k \rm(km)');ylim([-100 600])
legend("\itSOC_b_a_t","\itSOC_R_F_C","\itR_b_a_c_k","Location","east")
set(gca,'Fontname', 'Times New Roman','FontWeight','norm','FontSize',10)
hold on;

subplot(3,1,3);xlabel('\itt \rm(h)');xlim([0 xlim_max])
yyaxis left;hold on;plot(t,P_plot(:,1),"-");plot(t,P_plot(:,2),':',LineWidth=2);ylabel('\itP \rm(kW)');ylim([0 200]);
yyaxis right;plot(t,P_plot(:,3),"-");plot(t,P_plot(:,4),':',LineWidth=2);ylabel('\itP \rm(kW)');ylim([-150 50]);
legend('\itP_s_o_l','\itP_r_e_q',"\itP_b_a_t","\itP_R_F_C","Location","east")
set(gca,'Fontname', 'Times New Roman','FontWeight','norm','FontSize',10)


figure('Name',[name '2'])
subplot(2,1,1);xlim([24 48])
yyaxis left;plot(t,v_plot,"-");ylabel('\itv \rm(m·s^-^1)');ylim([0 35])
yyaxis right;hold on;plot(t,P_plot(:,1),"-");ylabel('\itP_s_o_l \rm(kW)');ylim([0 200]);
%legend("\itv",'\itP_s_o_l',"Location","east")
set(gca,'Fontname', 'Times New Roman','FontWeight','norm','FontSize',10)
hold on;

subplot(2,1,2);xlabel('\itt \rm(h)');xlim([24 48])
yyaxis left;plot(t,v_plot,"-");ylabel('\itv \rm(m·s^-^1)');ylim([0 35])
yyaxis right;plot(t,EnvConstants.Interp.tV2T_He(t_24h,v_plot),"-");ylabel('\itT_H_e \rm(K)');
%legend("\itv",'$\dot m_H_e',"Location","east")
set(gca,'Fontname', 'Times New Roman','FontWeight','norm','FontSize',10)

stat.final_battery_energy=x_plot(end,1);
stat.final_RFC_energy=x_plot(end,2);
stat.final_withdraw=x_plot(end,3);
stat.gas_loss=x_plot(end,4)-x_plot(1,4);
stat.average_volecity=mean(u_act(:,1));
stat.accumulated_distance_beyond_residency=Ts*sum((x_plot(:,3)>EnvConstants.r_max).*(x_plot(:,3)-EnvConstants.r_max));%km·h
stat.total_battery_consumption=sum(max(diff(x_plot(:,1)),0));
stat.total_RFC_consumption=sum(max(diff(x_plot(:,2)),0));
stat.curtailment_rate=100*sum(EnvConstants.Interp.tV2Psol(t_24h,u_act(:,1))-P_plot(:,1))./sum(EnvConstants.Interp.tV2Psol(t_24h,u_act(:,1)));%
stat.total_cost=J_plot(end);
   

end



function [NextObs,Reward,IsDone,LoggedSignals] = EM_step(Action,LoggedSignals,EnvConstants)
% step function to construct energy management environment for the function handle case.
% This function applies the given action to the environment and evaluates the system dynamics for one simulation step.

% Unpack the state vector from the logged signals.
Design=EnvConstants.Design;
Interp=EnvConstants.Interp;
Ts=EnvConstants.Ts;
[atm.T, ~, atm.p, atm.rho] = atmoscoesa(EnvConstants.H);%1976 COESA model
v=Action(:,1);
r_bat=Action(:,2);
t_24h=(LoggedSignals.t-24*floor(LoggedSignals.t/24));

E_bat=LoggedSignals.State(:,1);
E_RFC=LoggedSignals.State(:,2);
R_back=LoggedSignals.State(:,3);
m_He=LoggedSignals.State(:,4);

% Check if the given action is valid.
if any(v<0) || any(v>EnvConstants.Design.v_max) || any(r_bat<0) || any(r_bat>1)
    error('Action value error.');
end

% Apply motion equations.
% 能量平衡：P_bat+P_RFC+P_sol=P_req,P_bat=r_bat*(P_bat+P_RFC),未满充时P_sol=P_sol_max，都满充后P_sol匹配负载
P_sol=Interp.tV2Psol(t_24h,v);% kW
dm_He=Interp.tV2Hedm(t_24h,v);% kg/s
P_req=Interp.V2Preq(v);%kW总消耗功率
P_bat=r_bat.*(P_req-P_sol);
P_RFC=(1-r_bat).*(P_req-P_sol);


IsFall=false(size(Action,1),1);%是否因某原因任务失败

P_bat_lb=max(-Design.P_bat_max,(E_bat-Design.E_bat_max)./Ts/Design.eff_chrg_bat);%锂电池功率低边界（最大充电功率），倍率<1C
P_bat_ub=min(Design.P_bat_max,(E_bat)./Ts*Design.eff_dischrg_bat);%锂电池功率高边界（最大放电功率），倍率<1C
P_RFC_lb=max(-Design.P_RFC_max,Interp.dEdt2P_RFC((Design.E_RFC_max-E_RFC)./Ts));%RFC功率低边界（最大充电功率），倍率<1C
P_RFC_ub=min(Design.P_RFC_max,Interp.dEdt2P_RFC(-E_RFC./Ts));%RFC功率高边界（最大放电功率），倍率<1C

% 以下情况r_bat不适用:
% if P_bat<P_bat_lb %bat充电功率太大或满充
%     if P_RFC<P_RFC_lb %都满充
%         P_bat=P_bat_lb;
%         P_RFC=P_RFC_lb;
%         P_sol=P_req-P_bat-P_RFC;
%     elseif P_RFC>P_RFC_ub %RFC放电功率太大或空
%         IsFall = true; %不允许的代价
% 
%     else %RFC正常
%         P_bat=P_bat_lb;
%         P_RFC=max(P_req-P_sol-P_bat,P_RFC_lb);%算上bat多余功率
%         P_sol=P_req-P_bat-P_RFC;
%     end
% elseif P_bat>P_bat_ub %bat放电功率太大或放空
%     if P_RFC<P_RFC_lb %RFC充电功率太大或满充
%         IsFall = true; %不允许的代价
%     elseif P_RFC>P_RFC_ub %RFC放电功率太大或空
%         IsFall = true; %不允许的代价
%     else %RFC正常
%         P_bat=P_bat_ub;
%         P_RFC=P_req-P_sol-P_bat;
%         if P_RFC>P_RFC_ub
%             IsFall = true; %不允许的代价
%         end
%     end
% 
% else %bat正常
%     if P_RFC<P_RFC_lb %RFC充电功率太大或满充
%         P_RFC=P_RFC_lb;
%         P_bat=max(P_req-P_sol-P_RFC,P_bat_lb);%算上RFC多余功率
%     elseif P_RFC>P_RFC_ub %RFC放电功率太大或空
%         P_RFC=P_RFC_ub;
%         P_bat=P_req-P_sol-P_RFC;
%         if P_bat>P_bat_ub
%             IsFall = true; %不允许的代价
%         end
%     end
% end
% 限幅：
P_bat(P_bat<P_bat_lb)=P_bat_lb(P_bat<P_bat_lb);
P_bat(P_bat>P_bat_ub)=P_bat_ub(P_bat>P_bat_ub);
P_RFC(P_RFC<P_RFC_lb)=P_RFC_lb(P_RFC<P_RFC_lb);
P_RFC(P_RFC>P_RFC_ub)=P_RFC_ub(P_RFC>P_RFC_ub);
% 补充：
P_bat(P_bat+P_RFC<P_req-P_sol & P_bat<P_bat_ub)=P_req(P_bat+P_RFC<P_req-P_sol & P_bat<P_bat_ub)-P_sol(P_bat+P_RFC<P_req-P_sol & P_bat<P_bat_ub)-P_RFC(P_bat+P_RFC<P_req-P_sol & P_bat<P_bat_ub);
P_RFC(P_bat+P_RFC<P_req-P_sol & P_RFC<P_RFC_ub)=P_req(P_bat+P_RFC<P_req-P_sol & P_RFC<P_RFC_ub)-P_sol(P_bat+P_RFC<P_req-P_sol & P_RFC<P_RFC_ub)-P_bat(P_bat+P_RFC<P_req-P_sol & P_RFC<P_RFC_ub);
P_bat(P_bat+P_RFC>P_req-P_sol & P_bat>P_bat_lb)=P_req(P_bat+P_RFC>P_req-P_sol & P_bat>P_bat_lb)-P_sol(P_bat+P_RFC>P_req-P_sol & P_bat>P_bat_lb)-P_RFC(P_bat+P_RFC>P_req-P_sol & P_bat>P_bat_lb);
P_RFC(P_bat+P_RFC>P_req-P_sol & P_RFC>P_RFC_lb)=P_req(P_bat+P_RFC>P_req-P_sol & P_RFC>P_RFC_lb)-P_sol(P_bat+P_RFC>P_req-P_sol & P_RFC>P_RFC_lb)-P_bat(P_bat+P_RFC>P_req-P_sol & P_RFC>P_RFC_lb);
% 补充后限幅：
P_bat(P_bat<P_bat_lb)=P_bat_lb(P_bat<P_bat_lb);
P_bat(P_bat>P_bat_ub)=P_bat_ub(P_bat>P_bat_ub);
P_RFC(P_RFC<P_RFC_lb)=P_RFC_lb(P_RFC<P_RFC_lb);
P_RFC(P_RFC>P_RFC_ub)=P_RFC_ub(P_RFC>P_RFC_ub);
% 满充太阳能减少：
P_sol(P_bat+P_RFC>P_req-P_sol)=P_req(P_bat+P_RFC>P_req-P_sol)-P_bat(P_bat+P_RFC>P_req-P_sol)-P_RFC(P_bat+P_RFC>P_req-P_sol);
%IsFall(abs(P_req-P_sol-P_bat-P_RFC)>1e-6)=true;
IsFall(P_req-P_sol-P_bat-P_RFC>20)=true;
% 
% %bat充电功率太大或满充:
% P_bat(P_bat<P_bat_lb & P_RFC<P_RFC_lb)=P_bat_lb(P_bat<P_bat_lb & P_RFC<P_RFC_lb);%都满充
% P_RFC(P_bat<P_bat_lb & P_RFC<P_RFC_lb)=P_RFC_lb(P_bat<P_bat_lb & P_RFC<P_RFC_lb);%都满充
% IsFall(P_bat<P_bat_lb & P_RFC>P_RFC_ub)=true;%RFC放电功率太大或空
% P_bat(P_bat<P_bat_lb & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub)=P_bat_lb(P_bat<P_bat_lb & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub);%bat充电功率太大或满充,%RFC正常算上bat多余功率
% P_RFC(P_bat<P_bat_lb & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub)=max(P_req(P_bat<P_bat_lb & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub)-P_sol(P_bat<P_bat_lb & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub)-P_bat(P_bat<P_bat_lb & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub),P_RFC_lb(P_bat<P_bat_lb & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub));%算上bat多余功率
% %bat放电功率太大或放空:
% IsFall(P_bat>P_bat_ub & P_RFC<P_RFC_lb) = true; %bat放电功率太大或放空%RFC充电功率太大或满充
% IsFall(P_bat>P_bat_ub & P_RFC>P_RFC_ub) = true; %都空
% P_bat(P_bat>P_bat_ub & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub)=P_bat_ub(P_bat>P_bat_ub & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub);
% P_RFC(P_bat>P_bat_ub & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub)=P_req(P_bat>P_bat_ub & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub)-P_sol(P_bat>P_bat_ub & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub)-P_bat(P_bat>P_bat_ub & P_RFC>=P_RFC_lb & P_RFC<=P_RFC_ub);
% IsFall(P_bat>P_bat_ub & P_RFC>P_RFC_ub)=true;
% %bat正常:
% P_RFC(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC<P_RFC_lb)=P_RFC_lb(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC<P_RFC_lb);
% P_bat(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC<P_RFC_lb)=max(P_req(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC<P_RFC_lb)-P_sol(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC<P_RFC_lb)-P_RFC(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC<P_RFC_lb),P_bat_lb(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC<P_RFC_lb));%算上RFC多余功率
% 
% P_RFC(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC>P_RFC_ub)=P_RFC_ub(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC>P_RFC_ub);
% P_bat(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC>P_RFC_ub)=P_req(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC>P_RFC_ub)-P_sol(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC>P_RFC_ub)-P_RFC(P_bat>=P_bat_lb & P_bat<=P_bat_ub & P_RFC>P_RFC_ub);
% IsFall(P_RFC>P_RFC_ub & P_bat>P_bat_ub)=true;



dEdt_bat=-(P_bat>=0).*P_bat/Design.eff_dischrg_bat-(P_bat<0).*P_bat*Design.eff_chrg_bat;
dEdt_RFC=Interp.P2dEdt_RFC(P_RFC);
dRdt_back=max((EnvConstants.Interp.wind(LoggedSignals.t)-v)*3.6,(-EnvConstants.r_max-R_back)/Ts);%因退飞的驻留点累积偏移量km/h，风速wind，若正飞到-50则拐弯;
dmdt_He=dm_He;
r_bat_actual=P_bat./(P_bat+P_RFC);%这里0/0=nan
r_bat_actual(isnan(r_bat_actual))=0;%把0/0的nan变为0


% Perform Euler integration.
t_previous=LoggedSignals.t;
State_previous=LoggedSignals.State;
LoggedSignals.t=LoggedSignals.t+Ts;
LoggedSignals.State = LoggedSignals.State + Ts.*[dEdt_bat,dEdt_RFC,dRdt_back,dmdt_He];
LoggedSignals.Action_actual=[v,r_bat_actual];%u=[v,r_bat];
LoggedSignals.P=[P_sol,P_req,P_bat,P_RFC];


% Transform state to observation.
NextObs = LoggedSignals.State;

% Get reward.
EnvConstants.J_temp_cal=@(t,x_cur,x_next) max(x_cur(:,1)-x_next(:,1),0)/200+max(x_cur(:,2)-x_next(:,2),0)/200+(max(0,x_cur(:,3)-EnvConstants.r_max)+max(0,x_next(:,3)-EnvConstants.r_max))/108+max(0,x_cur(:,4)-x_next(:,4))*24;

% if ~IsFall && LoggedSignals.t<=EnvConstants.t(end)
%     IsDone=false;
%     Reward = -EnvConstants.J_temp_cal(t_previous,State_previous,LoggedSignals.State);
% else    
%     IsDone=true;
%     Reward = EnvConstants.PenaltyForFalling;
% end
IsDone=IsFall | LoggedSignals.t>EnvConstants.t(end);
Reward=IsDone.*EnvConstants.PenaltyForFalling+(~IsDone).*-EnvConstants.J_temp_cal(t_previous,State_previous,LoggedSignals.State);

end

function [InitialObservation, LoggedSignals] = EM_reset(EnvConstants)
% Reset function to place EM environment.

init=EnvConstants.init;

% Return initial environment state variables as logged signals.
LoggedSignals.t=EnvConstants.t(1);
LoggedSignals.State=[init.E_bat,init.E_RFC,init.R_back,init.m_He];
LoggedSignals.Action_actual=[15,0.1429];%u=[v,r_bat];
LoggedSignals.P=[0,0,0,0];%P=[P_sol,P_req,P_bat,P_RFC];

InitialObservation = LoggedSignals.State;
end


function [NextObs,Reward,IsDone,LoggedSignals] = EM_step_pre(Action,LoggedSignals,EnvConstants)
% step function to construct energy management environment for the function handle case.
% This function applies the given action to the environment and evaluates the system dynamics for one simulation step.

% Unpack the state vector from the logged signals.
Design=EnvConstants.Design;
Interp=EnvConstants.Interp;
Ts=EnvConstants.Ts;
[atm.T, ~, atm.p, atm.rho] = atmoscoesa(EnvConstants.H);%1976 COESA model
v=Action(:,1);
r_bat=Action(:,2);
t_24h=(LoggedSignals.t-24*floor(LoggedSignals.t/24));

E_bat=LoggedSignals.State(:,1);
E_RFC=LoggedSignals.State(:,2);
R_back=LoggedSignals.State(:,3);
m_He=LoggedSignals.State(:,4);

% Check if the given action is valid.
if any(v<0) || any(v>EnvConstants.Design.v_max) || any(r_bat<0) || any(r_bat>1)
    error('Action value error.');
end

% Apply motion equations.
% 能量平衡：P_bat+P_RFC+P_sol=P_req,P_bat=r_bat*(P_bat+P_RFC),未满充时P_sol=P_sol_max，都满充后P_sol匹配负载
P_sol=Interp.tV2Psol(t_24h,v);% kW
dm_He=Interp.tV2Hedm(t_24h,v);% kg/s
P_req=Interp.V2Preq(v);%kW总消耗功率
P_bat=r_bat.*(P_req-P_sol);
P_RFC=(1-r_bat).*(P_req-P_sol);


IsFall=false(size(Action,1),1);%是否因某原因任务失败

P_bat_lb=max(-Design.P_bat_max,(E_bat-Design.E_bat_max)./Ts/Design.eff_chrg_bat);%锂电池功率低边界（最大充电功率），倍率<1C
P_bat_ub=min(Design.P_bat_max,(E_bat)./Ts*Design.eff_dischrg_bat);%锂电池功率高边界（最大放电功率），倍率<1C
P_RFC_lb=max(-Design.P_RFC_max,Interp.dEdt2P_RFC((Design.E_RFC_max-E_RFC)./Ts));%RFC功率低边界（最大充电功率），倍率<1C
P_RFC_ub=min(Design.P_RFC_max,Interp.dEdt2P_RFC(-E_RFC./Ts));%RFC功率高边界（最大放电功率），倍率<1C

% 以下情况r_bat不适用:

% 限幅：
P_bat(P_bat<P_bat_lb)=P_bat_lb(P_bat<P_bat_lb);
P_bat(P_bat>P_bat_ub)=P_bat_ub(P_bat>P_bat_ub);
P_RFC(P_RFC<P_RFC_lb)=P_RFC_lb(P_RFC<P_RFC_lb);
P_RFC(P_RFC>P_RFC_ub)=P_RFC_ub(P_RFC>P_RFC_ub);
% 补充：
P_bat(P_bat+P_RFC<P_req-P_sol & P_bat<P_bat_ub)=P_req(P_bat+P_RFC<P_req-P_sol & P_bat<P_bat_ub)-P_sol(P_bat+P_RFC<P_req-P_sol & P_bat<P_bat_ub)-P_RFC(P_bat+P_RFC<P_req-P_sol & P_bat<P_bat_ub);
P_RFC(P_bat+P_RFC<P_req-P_sol & P_RFC<P_RFC_ub)=P_req(P_bat+P_RFC<P_req-P_sol & P_RFC<P_RFC_ub)-P_sol(P_bat+P_RFC<P_req-P_sol & P_RFC<P_RFC_ub)-P_bat(P_bat+P_RFC<P_req-P_sol & P_RFC<P_RFC_ub);
P_bat(P_bat+P_RFC>P_req-P_sol & P_bat>P_bat_lb)=P_req(P_bat+P_RFC>P_req-P_sol & P_bat>P_bat_lb)-P_sol(P_bat+P_RFC>P_req-P_sol & P_bat>P_bat_lb)-P_RFC(P_bat+P_RFC>P_req-P_sol & P_bat>P_bat_lb);
P_RFC(P_bat+P_RFC>P_req-P_sol & P_RFC>P_RFC_lb)=P_req(P_bat+P_RFC>P_req-P_sol & P_RFC>P_RFC_lb)-P_sol(P_bat+P_RFC>P_req-P_sol & P_RFC>P_RFC_lb)-P_bat(P_bat+P_RFC>P_req-P_sol & P_RFC>P_RFC_lb);
% 补充后限幅：
P_bat(P_bat<P_bat_lb)=P_bat_lb(P_bat<P_bat_lb);
P_bat(P_bat>P_bat_ub)=P_bat_ub(P_bat>P_bat_ub);
P_RFC(P_RFC<P_RFC_lb)=P_RFC_lb(P_RFC<P_RFC_lb);
P_RFC(P_RFC>P_RFC_ub)=P_RFC_ub(P_RFC>P_RFC_ub);
% 满充太阳能减少：
P_sol(P_bat+P_RFC>P_req-P_sol)=P_req(P_bat+P_RFC>P_req-P_sol)-P_bat(P_bat+P_RFC>P_req-P_sol)-P_RFC(P_bat+P_RFC>P_req-P_sol);
IsFall(abs(P_req-P_sol-P_bat-P_RFC)>1e-6)=true;
% 


dEdt_bat=-(P_bat>=0).*P_bat/Design.eff_dischrg_bat-(P_bat<0).*P_bat*Design.eff_chrg_bat;
dEdt_RFC=Interp.P2dEdt_RFC(P_RFC);
dRdt_back=max((15-v)*3.6,(-EnvConstants.r_max-R_back)/Ts);%因退飞的驻留点累积偏移量km/h，风速wind，若正飞到-50则拐弯;
dmdt_He=dm_He;
r_bat_actual=P_bat./(P_bat+P_RFC);%这里0/0=nan
r_bat_actual(isnan(r_bat_actual))=0;%把0/0的nan变为0


% Perform Euler integration.
t_previous=LoggedSignals.t;
State_previous=LoggedSignals.State;
LoggedSignals.t=LoggedSignals.t+Ts;
LoggedSignals.State = LoggedSignals.State + Ts.*[dEdt_bat,dEdt_RFC,dRdt_back,dmdt_He];
LoggedSignals.Action_actual=[v,r_bat_actual];%u=[v,r_bat];
LoggedSignals.P=[P_sol,P_req,P_bat,P_RFC];


% Transform state to observation.
NextObs = LoggedSignals.State;

% Get reward.
EnvConstants.J_temp_cal=@(t,x_cur,x_next) max(x_cur(:,1)-x_next(:,1),0)/200+max(x_cur(:,2)-x_next(:,2),0)/200+(max(0,x_cur(:,3)-EnvConstants.r_max)+max(0,x_next(:,3)-EnvConstants.r_max))/108+max(0,x_cur(:,4)-x_next(:,4))*24;

% if ~IsFall && LoggedSignals.t<=EnvConstants.t(end)
%     IsDone=false;
%     Reward = -EnvConstants.J_temp_cal(t_previous,State_previous,LoggedSignals.State);
% else    
%     IsDone=true;
%     Reward = EnvConstants.PenaltyForFalling;
% end
IsDone=IsFall | LoggedSignals.t>EnvConstants.t(end);
Reward=IsDone.*EnvConstants.PenaltyForFalling+(~IsDone).*-EnvConstants.J_temp_cal(t_previous,State_previous,LoggedSignals.State);

end



