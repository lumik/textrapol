function h=textrapol_functions

h={@textrapol_load, @faktorka, @func, @residuum, @strfromvec, @lmmin};

function [status,filename,file_path,spectra]=textrapol_load(do_test,start_directory,...
FilterSpec,dialog_title,multiselect)
%--------------------------------------------------------------------------
% Nacteni spekter ze souboru (spektra ve sloupcich, 1. sloupec je x-ova skala) 
%--------------------------------------------------------------------------
% Syntaxe funkce:
% [status,filename,file_path,spectra]=readdata(do_test,start_directory,...
% FilterSpec,dialog_title,multiselect)
%--------------------------------------------------------------------------
% Vstupni parametry:
% do_test -> promenna do_test nabyva hodnoty 0,1. Pro hodnotu 1 se testuje,
% jestli je x-ova skala setridena (viz poznamka).
% start_directory -> pocatecni adresar, ktery se vypise v dialogovem okne
% FilterSpec -> 'Cell Array' pro specifikaci typu souboru pro otevreni (tak 
% jak je definovano ve funkci Matlabu 'uigetfile') 
% dialog_title -> retezec, ktery se vypisuje v dialogovem okne 
% multiselect -> 1: vybirame vice souboru, 0: lze vybrat pouze jeden soubor
%--------------------------------------------------------------------------
% Vystupni parametry:
% status -> v pripade stisku cancel pri nacitani dat ma status hodnotu 0.
% V opacnem pripade ma status hodnotu 1.
% filename -> jmeno nacteneho souboru (v pripade, ze se vybira vice souboru
% v rezimu multiselect=1, jsou jmena nactenych souboru v promenne typu Cell
% array).
% file_path -> absolutni cesta k souboru
% spectra -> sloupce v matici spectra odpovidaji jednotlivym nactenym spektrum 
% (prvni sloupec v matici je x-ova skala). v pripade stisku cancel pri
% nacitani dat je spectra=[]
%--------------------------------------------------------------------------
% Poznamka: 1) Pokud prvni sloupec dat (x-ove hodnoty) neni setriden (coz muze
% znacit chybu), pak je nabidnuto jeho vzestupne setrideni. Podle tohoho
% prvniho sloupce se pak setridi ostatni sloupce (spektralni intenzity). S
% takto setridenymi spektry se lepe pracuje pri vykreslovani grafu,
% orezavani, apod.
% 2)  Nacitana Spektra mohou byt ulozena v textovem souboru nebo v binarnim
% souboru (format matlabu: pripona mat). 
%--------------------------------------------------------------------------

% Vyber souboru:
[filename,file_path] = uigetfile(FilterSpec,dialog_title,...
fullfile(start_directory,'*.*'),'Multiselect',multiselect);
if (isequal(filename,0) || isequal(file_path,0)) % stisknuto cancel
 status=0;
 spectra=[];
else
 status=1;   
 spec=strcat(file_path,filename); % specifikace pro nasledne nacteni dat ze souboru
 try
  spectra=load(spec); % nacteni spekter ze souboru (textovy nebo binarni soubor "mat").
 catch
  status=0;
  filename=0;
  file_path=0;
  spectra=[];
  return
 end
 extension=filename((strfind(filename,'.')+1):1:end); % pripona souboru
 if strcmp(extension,'mat') % Spektra se nactou z "mat" souboru 
  polozky_spektra=fieldnames(spectra);   
  spectra=eval(['spectra.' polozky_spektra{1}]); % v souboru "mat" je ulozena
  % pouze jedna promenna (pozadovana spektra)
 end
 if do_test==1 &&  ~issorted(spectra(:,1)) && ~issorted(flipdim(spectra(:,1),1))
  % x-ove hodnoty nejsou setrideny ani vzestupne, ani sestupne -> muze 
  % signalizovat chybu (je nabidnuta moznost spektra vzestupne setridit)
  vypis=['Spectra have unsorted x-values (probably due to the error). '...
  'Do you want to sort them in increasing order ?'];    
  tridit = questdlg(vypis,'Loading data','Yes','No','Yes');
  if strcmp(tridit,'Yes')
   [spectra(:,1),indexy]=sort(spectra(:,1)); % tridi se x-ovy hodnoty
   for i=2:1:size(spectra,2) % y-ovy hodnoty pro kazdy spektrum se tridi podle x
    y=spectra(:,i);
    spectra(:,i)=y(indexy);
   end        
  end
 end
end

%% funkce pro vypocet faktorove analyzy
function [U,W,V,E]=faktorka(spektra)
%--------------------------------------------------------------------------
% Faktorova analyza serie spekter
%--------------------------------------------------------------------------
% Syntaxe funkce: [U,W,V,E]=faktorka(spektra)
%--------------------------------------------------------------------------
% Vstupni parametry:
% spektra -> sloupce teto matice odpovidaji exp. spektrum (bez x-oveho
% sloupce)
%--------------------------------------------------------------------------
% Vystupni parametry:
% U -> matice, jejiz sloupce jsou subspektra 
% V -> matice, jejiz sloupce jsou koeficienty
% W -> sloupec singularnich hodnot
% E -> sloupec rezidualnich chyb
%--------------------------------------------------------------------------

% Singular Value Decomposition (SVD). Urceni subspekter U, koeficientu V a
% singularnich hodnot W:
[U,W,V] = svd(spektra,'econ'); 
W=diag(W); % vektor singularnich hodnot

% vypocet rezidualni chyby E v zavislosti na faktorove dimenzi: 
spektra_pocet=size(spektra,2); % pocet spekter  
pocet=size(spektra,1); % pocet spektralnich bodu v jednom spektru
SK=W.^2;% Kvadrat residualni chyby
if spektra_pocet==1 
 E=0; % Pro jedno spektrum nelze definovat rezidualni chybu. E je v tomto
 % pripade definovano jako 0.
else 
 pocet_subspekter=size(U,2); % V pripade, ze pocet spektrálnich bodu >=pocet  
 % spekter, odpovida pocet spekter poctu subspekter. Potom funkce
 % svd(spektra,'econ') je ekvivalentni svd(spektra,0). V opacnem pripade 
 % svd(spektra,'econ') poskytuje mensi pocet subspekter. Pocet subspekter je
 % roven mensimu z parametru (pocet spekter, pocet spektralnich bodu)
 E=zeros(pocet_subspekter-1,1);
 for m=1:(pocet_subspekter-1)
  vyraz_1=sum(SK((m+1):(pocet_subspekter)));
  vyraz_2=pocet*((pocet_subspekter)-m);
  E(m)=sqrt(vyraz_1/vyraz_2); % hodnota rezidualni chyby
 end
end

function y=func(x, p)
y = p(1) * exp(p(2) * (x - p(3)) ) + ...
    p(4) * exp(p(5) * (x - p(3)) ) + p(6) * x + p(7);

function resid = residuum(func, x, V, W, factor_num, p)

M = size(V,1);

pind = 1:7;
resid = zeros(M*length(factor_num) ,1);
for ii=1:length(factor_num)
    if ii > 1
        pind(1) = 8 + (ii-2) * 4;
        pind(4) = 9 + (ii-2) * 4;
        pind(6) = 10 + (ii-2) * 4;
        pind(7) = 11 + (ii-2) * 4;
    end
    size(resid((ii-1)*M + 1:ii*M))
    size(V(:,factor_num(ii)))
    size(func(x, p(pind)))
    resid((ii-1)*M + 1:ii*M) = (V(:,factor_num(ii)) - func(x, p(pind))) / W(factor_num(ii));
end

function s = strfromvec(v)

if length(v) == 1
    s = num2str(v);
    return;
end

s = '';
jj = 1;
kk = 1;
brackets = false;
for ii = 2:length(v)
    if v(ii) ~= jj + 1
        if jj == kk
            if ~isempty(s)
                s = strcat(s, ', ');
                brackets = true;
            end
            s = strcat(s, num2str(jj));
        elseif jj - 1 == kk
            if ~isempty(s)
                s = strcat(s, ', ');
            end
            brackets = true;
            s = strcat(s, num2str(kk), ', ', num2str(jj));
        else
            if ~isempty(s)
                s = strcat(s, ', ');
                brackets = true;
            end
            s = strcat(s, num2str(kk), ':', num2str(jj));
        end
        kk = v(ii);
    end
    jj = v(ii);
end
if jj == kk
    if ~isempty(s)
        s = strcat(s, ', ');
        brackets = true;
    end
    s = strcat(s, num2str(jj));
elseif jj - 1 == kk
    if ~isempty(s)
        s = strcat(s, ', ');
    end
    brackets = true;
    s = strcat(s, num2str(kk), ', ', num2str(jj));
else
    if ~isempty(s)
        s = strcat(s, ', ');
        brackets = true;
    end
    s = strcat(s, num2str(kk), ':', num2str(jj));
end

if brackets
    s = strcat('[',s,']');
end

function [p,sigma,iter]...
    =lmmin(func,dfunc,sumfunc,x,p0,maxit,eps,lambda,mu)
% funkce minalizuje funkci func vuci parametrum p
% Levenberg-Marquardtovou metodou
% pouziti: [p,sigma,pravdepodobnost,iter]...
% =lmmin(func,dfunc,x,p0,sigma,maxit,eps,lambda,mu);
% kde func = fitovana funkce
%     dfunc = funkce pocitajici Jacobian, pokud je zadana prazdna matice []
%             vypocte se jakobian numericky pomoci adaptivniho algoritmu
%             zalozeneho na Rombergove extrapolaci
%     sumfunc = funkce, ktera se pouzije na sumaci rezidui, durazne se
%               doporucuje suma ctvercu, aby spravne fungovala vsechna
%               statistika a normalni rozdeleni chyb
%     x = x-ove souradnice fitovanych bodu
%     p0 = odhad parametru
%     maxit = maximalni pocet iteraci, pocatecni hodnota je 100 (nepovinny parametr)
%     eps = max. odchylka sum rezidui nasledujich iteraci, pocatecni hodnota je 1e-6 (nepovinne)
%     lambda = tlumici faktor pocatecni hodnota je 1e-3 (damping faktor)
%     mu = korekcni faktor pro tlumici faktor pocatecni hodnota je 10
%
% funkce vraci dva sloupcove vektory - p je vektor parametru,
% sigma vektor smerodatnych odchylek, a pocet probehlych iteraci
%

% nastaveni ticheho rezimu, kdy se nic krome chybovych hlasek nevypisuje
% na konzoli:
display_results=0; % 1 pro vypsani vysledku, 0 pro mlceni
extent=1e1; % parametr pro funkci findump, ktery urcuje nejvyssi (+extent)
% a nejnizsi mocnitel (-extent), pomoci kterych se hleda nove lambda
% (lambda*mu^k, kde k=-extent:extent)
 
N=length(x);   % pocet dat
M=length(p0);  % pocet vypresnovanych parametru
 
p0=p0(:);      % timto udelam z p0 sloupcovy vektor
dfunc_control=1;% preinicializace indikatoru, zdali byl na zacatku zadan
                % Jakobian (hodnota 1), nebo se Jakobian bude pocitat
                % numericky (hodnota 0).
 
if (nargin<6) maxit=100;end; % maximalni pocet iteraci
if (nargin<7) eps=1e-6;end; % max. odchylka sum rezidui nasledujich iteraci 
if (nargin<8) lambda=1e-3;end; % tlumici faktor (damping faktor)
if (nargin<9) mu=10;end; % korekce pro tlumici faktor 

if ischar(func)
  func = str2func(func);
end

% pokud je dfunc prazdna matice, vypocte se jakobian numericky
if isempty(dfunc)
 dfunc=@(x,p) jacobianest(@(p) func(x,p),p,eps(1.0))';
 dfunc_control=0;
else
    if ischar(dfunc)
        dfunc = str2func(dfunc);
    end
end
 
% prealokuji pamet na matice a vektory potrebne pro vypocty
J=zeros(M,N);     % Jacobian (M radku, N sloupcu)
 
% nyni spoctu sumu rezidui pro prvotni odhad p0
S0=sumfunc(func(x,p0));
%  S0=sum((y-feval(func,x,p0)).^2);
 deltaS=1;  % prednastavena hodnota rozdilu rezidui
 iter=0;    % pocitadlo iteraci

 % prvni iterace
  p=p0;
  Sn=S0;
 % cyklus vypresnovani
 if ~dfunc_control && display_results
  fprintf('Iterace cislo:\n');
 elseif dfunc_control && display_results
  fprintf('Iterace cislo: 1');
 end
 while (iter<maxit && deltaS>eps);
    if ~dfunc_control && display_results
     fprintf('%d\n',iter+1);
     fprintf('Chy2=%e\n',Sn/(N-M));
     fprintf('deltaS=%e\n',deltaS);
     fprintf('lambda=%e\n',lambda);
    elseif dfunc_control && display_results
     if rem((iter+1),10)==0
      fprintf(',%d',iter+1)
     end
     if rem((iter+1+10),210)==0 && iter+1~=20
      fprintf('\n%d',iter+1)
     end
    end
    % cyklus je ukoncen po dosazeni max.poctu iteraci nebo zadane presnosti
    
    % spoctu matici Jacobianu
    J=dfunc(x,p);
    % spoctu pomocne matice LeveStrany a PraveStrany
    PS=J*(-func(x,p));
    LS0=J*J';
    [p,Sn,lambda]=findamp(func,sumfunc,x,p,PS,LS0,lambda,mu,extent,maxit);
    % odchylka sum rezidui
    deltaS=abs(S0-Sn);
    iter=iter+1;
    S0=Sn;
 end;
 % pokud je dosazeno maximalniho poctu iteraci, zobrazi se
 % varovani
 if (iter>=maxit)
   if dfunc_control && display_results
    fprintf('\n');
   end
   disp('Vypocet zastaven po dosazeni maximalniho poctu iteraci.');
   disp('Zrejme nebylo dosazeno minima.');
   disp('Zadejte jiny pocatecni odhad p0.')
 end;
 % cyklus ukoncen, pocitaji se smerodatne odchylky
 
 % pokud je chikvadrat vetsi nez 1, tak byly chyby nedocenene a tak je
 % treba normovat chyby tak, aby chykvadrat vysel 1
 nu=N-M; % pocet stupnu volnosti
 LS=J*J';
 sigma=diag(sqrt(inv(LS)));
 
 % vypocet rozdilu chy2 od jednicky
%   chi_vb=abs(sum(func(x,p).^2)/nu-1);
 
 
% vypocet pravdepodobnosti chy2
%  pravdepodobnost=1-chi2cdf(nu*(1+chi_vb),nu)+chi2cdf(nu*(1-chi_vb),nu);
% no, a to je vse, zbyva jen vypsat nejak uhledne vysledky
 if ~dfunc_control && display_results
  fprintf('---konec iterovani---\n');
 else
  %fprintf('\n');   
 end
 if display_results
  fprintf('\n');
  fprintf('Po %d iteracich jsem dospel k vysledku:\n',iter);
  for k=1:M
   fprintf('p(%d)=%f +/- %f\n',k,p(k),sigma(k));
  end;
  fprintf('Finalni damping factor je: %e\n\n',lambda);
 end

%%%%%%%%%%%%%%%%%%
% subfunkce pro nelezeni nejlepsiho lambda
%%%%%%%%%%%%%%%%%%

function [p,Sn,lambda]=findamp(func,sumfunc,x,p0,PS,LS0,lambda,mu,extent,maxit)
% [p,Sn,lambda]=findamp(func,x,p0,PS,LS0,lambda,lambda0,mu,extent,
% extent0,maxit)
% nalezne nejlepsi dumping factor lambda pro danou iteraci lmfitu:
%   func = fitovana funkce
%   x = x-ove souradnice fitovanych bodu
%   p0 = parametry pred iteraci
%   PS = vektor prave strany reseni problemu normalnich rovnic
%   LS0 = vektor leve strany reseni problemu normalnich rovnic pouze
%         Newtonovou metodou
%   lambda = tlumici faktor (damping faktor) pred iteraci
%   mu = korekcni faktor pro tlumici faktor
%   extent = urcuje nejvyssi a nejnizsi mocnitel (-extent), pomoci kterych
%         se hleda nove lambda (lambda*mu^k, kde k=-extent:extent)
%   maxit = maximalni pocet iteraci
%
% funkce vraci dva sloupcove vektory - p je vektor parametru,

M=length(p0);
S0=sumfunc(func(x,p0));
lambdan=zeros(2*extent+1,1);
dpn=zeros(M,2*extent+1);
S=ones(2*extent+1,1);
M=eye(M);

% napleni vektoru lambda faktoru, ktere se budou testovat na nejlepsi
% vysledek
for ii=-extent:extent
 index=ii+extent+1;
 lambdan(index)=lambda*mu^ii;
 LSn=LS0+M*lambdan(index);
 dpn(:,index)=LSn\PS;  % tj. inv(LS)*PS
 S(index)=sumfunc(func(x,p0+dpn(:,index)));
end

if isempty(find(S0>S,1)) % neuspel jsem, je treba hledat dal (lambda->lambda*mu^k nebo lambda*mu^-k)
 Sn1=S(end);
 k=extent+1;
 
 while (Sn1>=S0 & k<maxit+extent) % dokud jsem horsi nez predtim
  % pocitam jen levou stranu, prava zustava stejna 
  LS=LS0+M*lambda*mu^k;
  % a nyni uz pocitam posunuti parametru smerem k minimu
  dp1=LS\PS;  % tj. inv(LS)*PS
  Sn1=sumfunc(func(x,p0+dp1));
  k=k+1;
 end
 % nove parametry
 if Sn1<S0
  p=p0+dp1;
  lambda=lambda*mu^(k-1); % nova hodnota lambda
  Sn=Sn1;
 else
  Sn1=S(1);
  k=extent+1;
  while (Sn1>=S0 & k<maxit+extent) % dokud jsem horsi nez predtim
   % pocitam jen levou stranu, prava zustava stejna 
   LS=LS0+M*lambda*mu^-k;
   % a nyni uz pocitam posunuti parametru smerem k minimu
   dp1=LS\PS;  % tj. inv(LS)*PS
   Sn1=sumfunc(func(x,p0+dp1));
   k=k+1;
  end
  % nove parametry
  if Sn1<S0
   p=p0+dp1;
   lambda=lambda*mu^(k-1); % nova hodnota lambda
   Sn=Sn1;
  else % zbyva dodelat pripad, co kdyz se nedari nalezt zadne reseni
   p=p0;
   lambda=lambda;
   Sn=S0;
  end
 end
else % alespon jedno priblizeni je lepsi
 [Sn,ind]=min(S);
 p=p0+dpn(:,ind);
 lambda=lambdan(ind);
end

function J=jacobianest(func,p0,epsilon)
p0=p0(:);   % udela z p0 sloupcovy vektor
M = length(p0);
%epsilon = eps(1.0)
dp = sqrt(epsilon);
delta = zeros(M,1);
delta(1) = dp;
dy1 = (func(p0+delta)-func(p0-delta))/(2*dp);
delta(1) = 0;
N = length(dy1);
J = zeros(N,M);
J(:,1) = dy1;
for ii = 2:M
    delta(ii) = dp;
    J(:,ii) = (func(p0+delta)-func(p0-delta))/(2*dp);
    delta(ii)= 0;
end