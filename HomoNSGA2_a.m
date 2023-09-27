function HomoNSGA2_a
clear all
% make it pseudo-random
seed = sum(100*clock);
rand('twister', seed);
randn('state', seed);
NCC=[];
for innn=1:1
A=imread('ukbench09772.jpg'); A=rgb2gray(A); [M,N]=size(A);
B=imread('ukbench09775.jpg'); B=rgb2gray(B);
ptsOriginal  = detectSURFFeatures(A);
ptsDistorted = detectSURFFeatures(B);
[featuresOriginal,  validPtsOriginal]  = extractFeatures(A,  ptsOriginal);
[featuresDistorted, validPtsDistorted] = extractFeatures(B, ptsDistorted);
indexPairs = matchFeatures(featuresOriginal, featuresDistorted);
matchedOriginal  = validPtsOriginal(indexPairs(:,1));
matchedDistorted = validPtsDistorted(indexPairs(:,2));
X1=matchedOriginal.Location'; X2=matchedDistorted.Location';
X1o=ptsOriginal.Location'; X2o=ptsDistorted.Location';
%figure; 
tam1=size(X1o,2); tam2=size(X2o,2); tam=size(X1,2);
tam1=round(min(tam1,tam2)); in1=randperm(tam1); temp=X2o(:,in1); X2o=temp;
%showMatchedFeatures(A,B,matchedOriginal,matchedDistorted);
%title('Putatively matched points (including outliers)');
X = [X1,X1o(:,1:tam1); X2,X2o(:,1:tam1)];
% DIFFERENTIAL EVOLUTION: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nStarting DE'); tic;
% get minimal subset size and model codimension:
k=4; %Para la homografía.
fprintf('\nMinimal sample set dimension = %d', k);
% total number of elements
N = size(X, 2);
% get the noise threshold via Chi squared distribution:
T_noise_squared=chi2inv(0.9999,k);
fprintf('\nInitial Squared noise threshold = %f', T_noise_squared);
% Np=90;          %Tamano de la poblacion.90
% F=0.25;         %Factor de escalamiento
% Cr=0.8;         %Probabilidad de cruza
% tmax=300;       %Contador de iteraciones
% J=0; Jbest=0; it=0;
max_range=[N,N,N,N,chi2inv(0.9999,k)];%[tam/4,tam/2,3*tam/4,tam];
min_range=[1,1,1,1,chi2inv(0.1,k)]; 
d=5;
%% General parameters of bioinspired algorithm: 
M=2; 
Np=200; gen=200;    % tamano de poblacion, numero de generaciones
%% Initialize the population:
Pt = initPopH_a(Np, d, min_range, max_range,X);
%% Sort the initialized Pt (MultiObjective):
Pt = nonDomSort(Pt, M, d);
figure;
plot(Pt(:,d + 1),Pt(:,d + 2),'g*'); xlabel('B','FontSize',20),ylabel('\epsilon','FontSize',20);
%% Start the evolution process:
for i = 1 : gen
    %% Particular bioinspired algorithm (GA): 
    pool = round(Np/2); tour = 2; mu = 20; mum = 20;
    parents = tournamentSelection(Pt, pool, tour);
    Qt=GA_a(parents,M,d,mu,mum,min_range,max_range,X);
    %%% End of particular bioinspired algorithm.     
    %% Multiobjective part: 
    [sizePt,~] = size(Pt);      % size of ofiginal Pt
    [sizeQt,~] = size(Qt);      % size of secondary Pt
    % Rt is a concatenation of current Pt and the Qt:
    Rt(1:sizePt,:) = Pt;
    Rt(sizePt+1 : sizePt+sizeQt,1 : M+d) = Qt;
    % Non-domination-sort and Crowding distance of intermediate Pt:
    if size(Qt,1)==0 || size(Rt,1)==0
        Qt=Qt;
    end
    Rt = nonDomSort(Rt, M, d);
    if size(Rt,1)<Np
       Pt=Pt;
    else
       Pt = replacePop(Rt, M, d, Np); 
    end
    %Pt=Rt(1:Np,:);
    % Perform Selection:
    %
    %% Report of advances:
    if ~mod(i,5)
      fprintf('%d generations completed, tamano Pt: %d\n',i,size(Pt,1));
%       [temp, inn]=sort(Pt(:,d + 1)); Pt(:,d + 1)=temp;
%       temp=Pt(inn(1,1:end),d + 2); Pt(inn(1,1:end),d + 2)=temp;
      plot(Pt(:,d + 1),Pt(:,d + 2),'or','LineWidth',3); xlabel('B','FontSize',20),ylabel('\epsilon','FontSize',20);
      %refreshdata(h) % Evaluate y in the function workspace
      drawnow
    end
end
%% Visualize results:
% hold on
% if M == 2
%     plot(Pt(:,d + 1),Pt(:,d + 2),'k*');
% elseif M ==3
%     plot3(Pt(:,d + 1),Pt(:,d + 2),Pt(:,d + 3),'*');
% end
[~,in1]=sort(Pt(:,d + 2));
extremo1=Pt(in1(1),1:4); extremo2=Pt(in1(end),1:4); medio=Pt(in1(40),1:4);
[Theta1, ~] = estimate_homography(X, round(extremo1));
[Theta2, ~] = estimate_homography(X, round(extremo2));
[Theta3, ~] = estimate_homography(X, round(medio));
img1 = repmat(uint8(0),size(B)); [M,N]=size(A);% Create an m by n array of 0's
img2 = repmat(uint8(0),size(B));
img3 = repmat(uint8(0),size(B));
H1=reshape(Theta1, 3, 3);H2=reshape(Theta2, 3, 3);H3=reshape(Theta3, 3, 3);
for in1=1:M
    for in2=1:N
    CS1=(H1(2,1)*in2+H1(2,2)*in1+H1(2,3))/(H1(3,1)*in2+H1(3,2)*in1+H1(3,3));
    CS2=(H1(1,1)*in2+H1(1,2)*in1+H1(1,3))/(H1(3,1)*in2+H1(3,2)*in1+H1(3,3));
    CS1=round(CS1); CS2=round(CS2);
    % Solo considera los puntos que caen dentro de la imagen: %%%%%%%%%%%%%    
    if CS1>0 && CS1<=M && CS2>0 && CS2<=N        
        img1(CS1,CS2)=A(in1,in2);      
    end
    CS1=(H2(2,1)*in2+H2(2,2)*in1+H2(2,3))/(H2(3,1)*in2+H2(3,2)*in1+H2(3,3));
    CS2=(H2(1,1)*in2+H2(1,2)*in1+H2(1,3))/(H2(3,1)*in2+H2(3,2)*in1+H2(3,3));
    CS1=round(CS1); CS2=round(CS2);
    % Solo considera los puntos que caen dentro de la imagen: %%%%%%%%%%%%%    
    if CS1>0 && CS1<=M && CS2>0 && CS2<=N        
        img2(CS1,CS2)=A(in1,in2);      
    end
    CS1=(H3(2,1)*in2+H3(2,2)*in1+H3(2,3))/(H3(3,1)*in2+H3(3,2)*in1+H3(3,3));
    CS2=(H3(1,1)*in2+H3(1,2)*in1+H3(1,3))/(H3(3,1)*in2+H3(3,2)*in1+H3(3,3));
    CS1=round(CS1); CS2=round(CS2);
    % Solo considera los puntos que caen dentro de la imagen: %%%%%%%%%%%%%    
    if CS1>0 && CS1<=M && CS2>0 && CS2<=N        
        img3(CS1,CS2)=A(in1,in2);      
    end
    end
end
K1 = imadd(img1,B,'uint16');K2 = imadd(img2,B,'uint16');K3 = imadd(img3,B,'uint16');
NK1 = NormalizedCrossCorrelation(img1, B);
NK2 = NormalizedCrossCorrelation(img2, B);
NK3 = NormalizedCrossCorrelation(img3, B);
NCC(innn,:)=[NK1,NK2,NK3];
end
figure,imshow(K1,[]),title(num2str(NK1)),
figure,imshow(K2,[]),title(num2str(NK2)),
figure,imshow(K3,[]),title(num2str(NK3)),
% % h1=figure;
% % imshow(K1,[]);
% % ti = get(gca,'TightInset');
% % set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
% % set(gca,'units','centimeters')
% % pos = get(gca,'Position');
% % ti = get(gca,'TightInset');
% % set(gcf, 'PaperUnits','centimeters');
% % set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% % set(gcf, 'PaperPositionMode', 'manual');
% % set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% % saveas(h1,'Fig5aGA.eps');
% % h1=figure;
% % imshow(K2,[]);
% % ti = get(gca,'TightInset');
% % set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
% % set(gca,'units','centimeters')
% % pos = get(gca,'Position');
% % ti = get(gca,'TightInset');
% % set(gcf, 'PaperUnits','centimeters');
% % set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% % set(gcf, 'PaperPositionMode', 'manual');
% % set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% % saveas(h1,'Fig5bGA.eps');
% % h1=figure;
% % imshow(K3,[]);
% % ti = get(gca,'TightInset');
% % set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
% % set(gca,'units','centimeters')
% % pos = get(gca,'Position');
% % ti = get(gca,'TightInset');
% % set(gcf, 'PaperUnits','centimeters');
% % set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% % set(gcf, 'PaperPositionMode', 'manual');
% % set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
% % saveas(h1,'Fig5cGA.eps');
% % imshow(A),figure,imshow(K,[]),figure,imshow(B)

