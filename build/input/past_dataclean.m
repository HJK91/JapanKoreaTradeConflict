clc
clear

%% Load the data

data = readtable('KLoWF_ind.csv');
idx=isnan(table2array(data(:,'P18CN02')));

data(idx,:)=[];

%% Label the data
% ID 
OPID = data.OPID; % Individual id numbers
Pwave = data.Pwave; % wave of surveys
age = data.P99AG;

% Education before schooling
ischild = data.P18CN02;
rel_child = data.P18RCA;
secondary = data.P18CA23A;
facility = data.P18CC01A;
alonetime= data.P18AT01A;

% Relationship with parents/husband's parents
alive = [data.P22AP data.P23AP];
Page = [data.P22AGG data.P22AGH data.P23AGI data.P23AGJ];
colive = [data.P22LT03 data.P23LT04];
dist = [data.P22HF06 data.P23HF07];
Pcare = [data.P22CA08 data.P23CA08];
livtgt= data.P10LA01;

% Labor supply, wage, consumption and income
Mworking = data.W01MD01;
sal = data.JOBSAL;
cons = data.HH04LE17;
faminc =data.HH03TI03;
asset = data.HH05FA14 + data.HH05TA11;

%% Clean the data from strings to numbers
% Pwave
for i = 1:7
    %work(work(:)=="실업자")={'0'}
    char = [num2str(i) '차'];
    idx = strcmp(Pwave,char);
    Pwave(idx)={i};
end
Pwave = [Pwave{:,:}]';

%% Parents alive
char_vec = {('두 분 모두 돌아가셨다');('어머님만 계신다');('아버님만 계신다');...
    ('두 분 모두 계신다')};
for i=1:size(char_vec,1)
    char = char_vec{i,:};
    idx = strcmp(alive,char);
    alive(idx) = {num2str(i)};
end
char_vec = {('두 분 모두 돌아가셨다');('시어머님만 계신다');('시아버님만 계신다');...
    ('두 분 모두 계신다')};
for i= 2:3
    char = char_vec{i,:};
    idx = strcmp(alive,char);
    alive(idx) = {num2str(i)};    
end
alive = str2double(alive);

%% Distance
char_vec = {('도보로 이동 가능한 거리');('차로 30분 이상 1시간 이내 거리');...
    ('차로 1시간 이상 2시간 이내 거리'); ('차로 2시간 이상 거리')};
for i=1:size(char_vec,1)
    char = char_vec{i,:};
    idx = strcmp(dist,char);
    dist(idx) = {num2str(i-1)};
end
dist = str2double(dist);

%% Colive 
char_vec = {('나와 함께 살고 계신다'); ('남편의 다른 형제자매와 함께 살고 계신다');...
    ('둘 다 아니다');('다른 형제자매와 함께 살고 계신다')};
char = char_vec{1,:};
idx = strcmp(colive,char);
colive(idx) = {num2str(1)};
colive(~idx) = {num2str(0)};
colive = str2double(colive);

%% Live together with husband 
char_vec = {('그렇다'); ('아니다')};
char = char_vec{1,:};
idx = strcmp(livtgt,char);
livtgt(idx) = {num2str(1)};
livtgt(~idx) = {num2str(0)};
livtgt = str2double(livtgt);

%% Binary
char_vec = unique(secondary);
char_vec = {char_vec{end-1,:};char_vec{end,:};char_vec{2,:}};
for i=1:size(char_vec,1)
    char = char_vec{i,:};
    idx_Pcare = strcmp(Pcare,char);
    idx_secondary = strcmp(secondary,char);
    Pcare(idx_Pcare) = {num2str(i-1)};
    secondary(idx_secondary) = {num2str(i-1)};
end
Pcare = str2double(Pcare);
secondary = str2double(secondary);

%% Relationship with child
char_vec = unique(rel_child);
char = char_vec{3,:};
idx = strcmp(rel_child,char);
rel_child(idx) = {num2str(1)};
rel_child(~idx) = {num2str(0)};
rel_child = str2double(rel_child);

%% Woman's labor supply
char_vec = unique(Mworking);
char1 = char_vec{12,:}; % employed
char0 = {char_vec{5,:};char_vec{13,:}}; % unemployed
idx1 = strcmp(Mworking,char1);
idx00 = strcmp(Mworking,char0{1,:});
idx01 = strcmp(Mworking,char0{2,:});
Mworking(idx1) = {num2str(1)};
Mworking(idx00) = {num2str(0)};
Mworking(idx01) = {num2str(0)};
Mworking(~(idx1+idx00+idx01)) = {num2str(2)};
Mworking = str2double(Mworking);

%% Facility 
char_vec = unique(facility);
char = {char_vec{6,:};char_vec{8,:};char_vec{9,:}}; 
% no response/don't know/not using any facility
idx1 = strcmp(facility,char{1,:});
idx2 = strcmp(facility,char{2,:});
idx3 = strcmp(facility,char{3,:});
facility(idx1) = {num2str(0)};
facility(idx2) = {num2str(0)};
facility(idx3) = {num2str(0)};
facility(~(idx1+idx2+idx3)) = {num2str(1)};
facility = str2double(facility);



%% clear redundant variables


 clear idx char ans char0 char1 char_vec idx00 idx01 idx1 idx_Pcare ...
     idx_secondary i idx2 idx3
save('KLoWF_cleaned.mat')
%% obtain the summary stats
load('KLoWF_cleaned.mat')

%%
% Numbers cases
MworkPcare      = sum(Mworking(Pcare(:,1)==1)==1);
MworkPcareN     = sum(Mworking(Pcare(:,1)==0)==1);
MworkHPcare     = sum(Mworking(Pcare(:,2)==1)==1);
MworkHPcareN    = sum(Mworking(Pcare(:,2)==0)==1);
MolfPcare       = sum(Mworking(Pcare(:,1)==1)==2);
MolfHPcare      = sum(Mworking(Pcare(:,2)==1)==2);
MolfPcareN      = sum(Mworking(Pcare(:,1)==0)==2);
MolfHPcareN     = sum(Mworking(Pcare(:,2)==0)==2);

MworkNocare     = sum(Mworking(Pcare(:,1)==0 | Pcare(:,2)==0)==1);
MolfNocare      = sum(Mworking(Pcare(:,1)==0 | Pcare(:,2)==0)==2);


% Labor market participation rates across waves

lpr = zeros(size(unique(Pwave),1),3);

for i = 1:size(lpr,1)
    
    lpr(i,1)=sum(Mworking(Pcare(:,1)==1 & Pwave==i)==1)/...
        (sum(Mworking(Pcare(:,1)==1 & Pwave==i)==1)+...
        sum(Mworking(Pcare(:,1)==1 & Pwave==i)==2));
    
    lpr(i,2)=sum(Mworking(Pcare(:,2)==1 & Pwave==i)==1)/...
        (sum(Mworking(Pcare(:,2)==1 & Pwave==i)==1)+...
        sum(Mworking(Pcare(:,2)==1 & Pwave==i)==2));
    
    lpr(i,3)=sum(Mworking(Pcare(:,1)==0 & Pwave==i)==1)/...
        (sum(Mworking(Pcare(:,1)==0 & Pwave==i)==1)+...
        sum(Mworking(Pcare(:,1)==0 & Pwave==i)==2));
    
    lpr(i,4)=sum(Mworking(Pcare(:,2)==0 & Pwave==i)==1)/...
        (sum(Mworking(Pcare(:,2)==0 & Pwave==i)==1)+...
        sum(Mworking(Pcare(:,2)==0 & Pwave==i)==2));    
    
    lpr(i,5)=sum(Mworking(Pwave==i)==1)/...
        (sum(Mworking( Pwave==i)==1)+...
        sum(Mworking( Pwave==i)==2));      
end

%%

f1 = figure(1);
hold on 
for i = [1, 3, 5]
    plot(unique(Pwave),lpr(:,i))
end
legend('Pcare','No Pcare','Whole sample','location','best')
xlabel('Wave')
ylabel('labor market participation rate')
saveas(f1,'fig1.png')
ylim([0 0.8])
hold off

f2 = figure(2);
hold on 
for i = [2, 4, 5]
    plot(unique(Pwave),lpr(:,i))
end
xlabel('Wave')
ylabel('labor market participation rate')
legend('HPcare','No HPcare','Whole sample','location','best')
saveas(f2,'fig2.png')
ylim([0 0.8])
hold off