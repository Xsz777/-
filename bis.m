%%  NN

%% ======================= Part 1: Plotting =======================
fprintf('プロット Data1 ...\n')
data1_A = csvread('data1-A.txt');
data1_B = csvread('data1-B.txt');
datal_X = csvread('data1-X.txt');

X_A = data1_A(:, 1);
Y_A = data1_A(:, 2);
X_B = data1_B(:, 1);
Y_B = data1_B(:, 2);
X_X = datal_X(:, 1);
Y_X = datal_X(:, 2);

plotData(X_A,Y_A,X_B,Y_B,X_X,Y_X);

legend('A社','B社','サンプル');
%% ======================= Part 2: ベイズ data1 =======================
Mean_A = mean(data1_A);
Mean_B = mean(data1_B);
Cov_A = cov(data1_A);
Cov_B = cov(data1_B);
fprintf('A社dataの平均:\n');
disp(Mean_A);
fprintf('B社dataの平均:\n');
disp(Mean_B);
fprintf('A社dataのCov:\n');
disp(Cov_A);
fprintf('B社dataのCov:\n');
disp(Cov_B);

fprintf('予测y:\n');
a = zeros(length(X_X),1);
b = zeros(length(X_X),1);

for i = 1 :length(X_X)
    a(i,1) = ((-1/2)*(datal_X(i,:) - Mean_A)) * inv(Cov_A) * (datal_X(i,:) - Mean_A)'-1/2 * log(abs(det(Cov_A)))-length(X_X)/2*log(2 * pi)+log(1/2);
    b(i,1) = ((-1/2)*(datal_X(i,:) - Mean_B)) * inv(Cov_B) * (datal_X(i,:) - Mean_B)'-1/2 * log(abs(det(Cov_B)))-length(X_X)/2*log(2 * pi)+log(1/2);
    if(a(i,1) > b(i,1))
        fprintf('P(Wa|A): %f P(Wb|B): %f',exp(a(i,1)),exp(b(i,1)));
        fprintf('サンプルはA社です\n');
    else
        fprintf('P(Wa|A): %f P(Wb|B): %f',exp(a(i,1)),exp(b(i,1)));
        fprintf('サンプルはB社です\n');
    end
end
fprintf('\nプログラム停止。Enterを押します。data2をインポ一トします。\n');
pause;

%% ======================= Part 3: ベイズ data2 =======================
fprintf('プロット Data2 ...\n')
data1_A = csvread('data2-A.txt');
data1_B = csvread('data2-B.txt');
datal_X = csvread('data2-X.txt');

X_A = data1_A(:, 1);
Y_A = data1_A(:, 2);
X_B = data1_B(:, 1);
Y_B = data1_B(:, 2);
X_X = datal_X(:, 1);
Y_X = datal_X(:, 2);

plotData(X_A,Y_A,X_B,Y_B,X_X,Y_X);

legend('A社','B社','サンプル');
Mean_A = mean(data1_A);
Mean_B = mean(data1_B);
Cov_A = cov(data1_A);
Cov_B = cov(data1_B);
fprintf('A社dataの平均:\n');
disp(Mean_A);
fprintf('B社dataの平均:\n');
disp(Mean_B);
fprintf('A社dataのCov:\n');
disp(Cov_A);
fprintf('B社dataのCov:\n');
disp(Cov_B);

fprintf('予测y:\n');
a = zeros(length(X_X),1);
b = zeros(length(X_X),1);

for i = 1 :length(X_X)
    a(i,1) = ((-1/2)*(datal_X(i,:) - Mean_A)) * inv(Cov_A) * (datal_X(i,:) - Mean_A)'
             -1/2 * log(abs(det(Cov_A)))-length(X_X)/2*log(2 * pi)+log(1/2);
    b(i,1) = ((-1/2)*(datal_X(i,:) - Mean_B)) * inv(Cov_B) * (datal_X(i,:) - Mean_B)'
             -1/2 * log(abs(det(Cov_B)))-length(X_X)/2*log(2 * pi)+log(1/2);
    if(a(i,1) > b(i,1))
        fprintf('P(Wa|A): %f P(Wb|B): %f',exp(a(i,1)),exp(b(i,1)));
        fprintf('サンプルはA社です\n');
    else
        fprintf('P(Wa|A): %f P(Wb|B): %f',exp(a(i,1)),exp(b(i,1)));
        fprintf('サンプルはB社です\n');
    end
end

%%
function plotData(x1,y1,x2,y2,x3,y3)
figure('name','Data'); % open a new figure window
hold on;
plot(x1,y1,'r.', 'MarkerSize', 10);
plot(x2,y2,'b.', 'MarkerSize', 10);
plot(x3,y3,'k.', 'MarkerSize', 10);
hold off;
end
