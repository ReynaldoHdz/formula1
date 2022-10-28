clear
clc
% Evidencia final - Simulación
% Reynaldo Hernández González A00829814


% Valores ingresados por usuario
mu = input("Ingresar valor de coeficiente de fricción (0.65 a 1.00): ");
masa = input("Ingresa la masa del automóvil (mínimo 703 kg): ");
velKh1 = input("Velocidad en km/h con la que el auto ingresa la primera zona de derrape: ");
velKh2 = input("Velocidad en km/h con la que el auto ingresa la segunda zona de derrape: ");

velMs1 = velKh1/3.6; %convertir en metros por segundo
velMs2 = velKh2/3.6; %convertir en metros por segundo

g = -9.81; %aceleración de gravedad

% PI(30,230) punto inicial
x1 = 30;
y1 = 230;

% P2(80,180) punto variado
x2 = 80;
y2 = 180;

% P3(150,10) punto variado
x3 = 150;
y3 = 10;

% PF(260,80) punto final
x4 = 260;
y4 = 80;


%% CASO CUBICO: a3*x^3 + a2*x^2 + a1*x + a0
% Definir vectores y matrices
A = [x1^3 x1^2 x1 1 ;
    x2^3 x2^2 x2 1;
    x3^3 x3^2 x3 1;
    x4^3 x4^2 x4 1];

y = [y1 ; y2 ; y3 ; y4];


% Calcular coeficientes
x = inv(A)*y; %x=A/y , Vector de coeficientes [a3 ; a2 ; a1 ; a0]


% Generar polinomio
% Podemos definir una funcion anonima para f(x)
fx = @(z) x(1).*(z.^3) + x(2).*(z.^2) + x(3).*z + x(4); %funcion
der = @(z) 3.*x(1).*(z.^2) + 2.*x(2).*z + x(3); %primera derivada
der1 = @(z) 6.*x(1).*z + 2.*x(2); %segunda derivada

%puntos_x = linspace(x1-50,x4+10,1000);
puntos_x = x1-50:0.1:x4+10;
puntos_y = fx(puntos_x);


%% Punto de inflexion
punto_max_y = max(puntos_y); %encontrar punto maximo en y
punto_min_y = min(puntos_y); %encontrar punto minimo en y
inflec_y = linspace(punto_min_y,punto_max_y,3); %vector de 3 puntos (max, medio, min)
punto_inflec_y = inflec_y(2); %punto medio en y

punto_inflec_x = 0; %inicializar variable
for i = 1:length(puntos_y)   
    if puntos_y(i) >= punto_inflec_y-1 && puntos_y(i) <= punto_inflec_y+1
        punto_inflec_x = puntos_x(i);
        break
    end
end


%% Encontrar longitud de la curva
derL = @(z) sqrt((der(z)).^2 + 1); %derivada de la funcion adaptada para longitud de curva
L = integral(derL,x1,x4);


%% Encontrar radio de curvatura
ind_pmaxy = find(puntos_y == punto_max_y); %encontramos el indice del punto max en y para encontrar x del punto max
ind_pminy = find(puntos_y == punto_min_y); %encontramos el indice del punto min en y para encontrar x del punto min
punto_max_x = puntos_x(ind_pmaxy); % x del punto max
punto_min_x = puntos_x(ind_pminy); % x del punto min
nummax = der(punto_max_x); %evaluamos la primera derivada usando punto maximo
denmax = der1(punto_max_x); %evaluamos la segunda derivada usando punto maximo
nummin = der(punto_min_x); %evaluamos la primera derivada usando punto minimo
denmin = der1(punto_min_x); %evaluamos la segunda derivada usando punto minimo

Rmax = ((1+(nummax)^2)^(3/2))/abs(denmax); %radio usando punto maximo
Rmin = ((1+(nummin)^2)^(3/2))/abs(denmin); %radio usando punto minimo

R = Rmax;

if Rmin < Rmax %se elige el radio mas pequeño de los dos
    R = Rmin;
end

% Iteración para verificación de valor de radio crítico en cada punto con
deltaT = 0.1;
radCritx = []; %inicializar vector
radCrity = []; %inicializar vector
for i = x1-50:deltaT:x4+10
    radio = ((1+(der(i))^2)^(3/2))/abs(der1(i)); %aplicando fórmula para encontrar radio
    if radio < 50
        radCritx = [radCritx, i]; %puntos x donde hay radio crítico
        radCrity = [radCrity, fx(i)]; %puntos y donde hay radio crítico
    end
end


%% Calcular velocidad máxima
v_max = []; %inicializar vector
for j = radCritx
    radio = ((1+(der(j))^2)^(3/2))/abs(der1(j)); %aplicando fórmula para encontrar radio
    vOptima = sqrt(mu*abs(g)*radio); %aplicando fórmula para encontrar velocidad optima
    v_max = [v_max, vOptima]; %las velocidades optimas en cada punto de la pista
end


%% Rectas tangentes
ind_rectasmax = ind_pmaxy-200:10:ind_pmaxy+150; %indices para puntos antes del maximo
ind_rectasmin = ind_pminy-200:10:ind_pminy+150; %indices para puntos antes del minimo
prxmax = puntos_x(ind_rectasmax); %vector de puntos en x antes del maximo
prxmin = puntos_x(ind_rectasmin); %vector de puntos en x antes del minimo
prymax = puntos_y(ind_rectasmax); %vector de puntos en y antes del maximo
prymin = puntos_y(ind_rectasmin); %vector de puntos en y antes del minimo
m1 = der(puntos_x); %valores de derivada de y con puntos antes del maximo

% rectas tangentes
for j = 1:length(prxmax)+1
    if j==length(prxmax)
        break
    end                    %el 20 determina la longitud de la recta tangente
    xmax = prxmax(j):prxmax(j)+20; %modificar longitud de rectas en punto max   
    atxmax = prxmax(j); %definir un punto x en max (cuando x=número)
    indmax = find(puntos_x==atxmax); %encontrar indice para cuando x=número en punto max
    bmax = puntos_y(indmax)-m1(indmax)*puntos_x(indmax); %encontrar b de y=mx+b en punto max
    yrectmax = m1(indmax).*xmax+bmax; %evaluar puntos cuando x=número usando y=mx+b en punto max

    xmin = prxmin(j):prxmin(j)+20; %modificar longitud de rectas en punto min
    atxmin = prxmin(j); %definir un punto x en min (cuando x=número)
    indmin = find(puntos_x==atxmin); %encontrar indice para cuando x=número en punto min
    bmin = puntos_y(indmin)-m1(indmin)*puntos_x(indmin); %encontrar b de y=mx+b en punto min
    yrectmin = m1(indmin).*xmin+bmin; %evaluar puntos cuando x=número usando y=mx+b en punto min
    
    % para las gradas:
    if atxmax==median(prxmax)-3.5
        pendmax = m1(indmax);
        pendmin = m1(indmin);
        equismax = xmax;
        equismin = xmin;    %hacemos copias de las rectas en un punto especifico para despues
        b_max = bmax;       %simplemente dezplazar la linea a un lado
        b_min = bmin;
        gradasmax = pendmax.*equismax+b_max;
        gradasmin = pendmin.*equismin+b_min;
    end
end


%% Evaluación de condicion de derrape para la zona de curvas
derrape = false;
k = 0;
% Curva 1
for h = 200:720 %posición 200 hasta posición 720 (para hacer una animación más rápida)
    if h == length(puntos_x)/2 %evitamos sobrepasar indice
        break
    end
    if ismembertol(puntos_x(h),radCritx) == true %checar si el punto actual es uno de radio crítico
        k = k + 1;
        if velMs1 > v_max(k) %si la velocidad actual es mayor a la velocidad máxima recomendada
            derrape = true;
            fprintf("El vehiculo se derrapa en el Punto de Derrape %d = [%0.1f, %0.3f]\n",k,radCritx(k),radCrity(k));
            tanVel = @(x) der(radCritx(k)).*x - (radCritx(k)*der(radCritx(k))) + fx(radCritx(k)); %encontrar vel tangencial
            T = -velMs1/(mu*g); %tiempo de recorrido
            s = 0.5*velMs1*T;
            dist = 0;
            incr = 0.1;
            while dist < s
                dist = sqrt(((radCritx(k)+incr) - radCritx(k))^2 + (tanVel(radCritx(k)+incr) - tanVel(radCritx(k)))^2);
                lim1 = radCritx(k) + incr;
                incr = incr + 0.1;
            end
            ptsDerrapeX = radCritx(k):0.1:lim1;
            ptsDerrapeY = tanVel(ptsDerrapeX);
            plot(puntos_x,puntos_y);
            hold on
            plot(equismax-25,gradasmax-3,'b','LineWidth',3)
            plot(equismin-25,gradasmin+3,'b','LineWidth',3)
            plot(puntos_x(h),puntos_y(h),'or')
            title("Simulación")
            xlabel("Eje x")
            ylabel("Eje y")
            ylim([punto_min_y-10 punto_max_y+10])
            plot(ptsDerrapeX,ptsDerrapeY,"color","r");
            calor = 0.5*masa*(velMs1^2);
            fprintf("La distancia recorrida por el auto al derraparse antes de detenerse es de %0.3f m\n",s);
            fprintf("La energía total perdida debido al derrape es de %0.3f J\n",calor);

            break
        end
    end
    % crea un punto que simula movimiento
    plot(puntos_x,puntos_y)
    hold on
    plot(equismax-25,gradasmax-3,'b','LineWidth',3)
    plot(equismin-25,gradasmin+3,'b','LineWidth',3)
    plot(puntos_x(h),puntos_y(h),'or')
    title("Simulación")
    xlabel("Eje x")
    ylabel("Eje y")
    ylim([punto_min_y-10 punto_max_y+10])
    pause(1/10000)
    if k ~= 718
        clf
    end
end

% Trayecto de en medio   
for aj = 721:20:2043 %esto es para hacer la animación más rápida
    if derrape == true %evitamos correr este for loop si ya se detectó un derrape
        break
    end
    plot(puntos_x,puntos_y)
    hold on
    plot(equismax-25,gradasmax-3,'b','LineWidth',3)
    plot(equismin-25,gradasmin+3,'b','LineWidth',3)
    plot(puntos_x(aj),puntos_y(aj),'or')
    ylim([punto_min_y-10 punto_max_y+10])
    pause(1/10000)
    clf
end

% Curva 2
h = (length(radCritx)/2)-0.5;
for k = 2044:length(puntos_x)-400 %posición 2044 hasta posición 2501 (para hacer una animación más rápida)
    if derrape == true %evitamos correr este for loop si ya se detectó un derrape
        break
    end
    if k == length(puntos_x)/2 %evitamos sobrepasar indice
        break
    end
    if ismembertol(puntos_x(k),radCritx) == true %checar si el punto actual es uno de radio crítico
        h = h + 1;
        if velMs2 > v_max(h) %si la velocidad actual es mayor a la velocidad máxima recomendada
            fprintf("El vehiculo se derrapa en el Punto de Derrape %d = [%0.1f, %0.3f]\n",h,radCritx(h),radCrity(h));
            tanVel = @(x) der(radCritx(h)).*x - (radCritx(h)*der(radCritx(h))) + fx(radCritx(h));
            T = -velMs2/(mu*g);
            s = 0.5*velMs2*T;
            dist = 0;
            incr = 0.1;
            while dist < s
                dist = sqrt(((radCritx(h)+incr) - radCritx(h))^2 + (tanVel(radCritx(h)+incr) - tanVel(radCritx(h)))^2);
                lim1 = radCritx(h) + incr;
                incr = incr + 0.1;
            end
            ptsDerrapeX = radCritx(h):0.1:lim1;
            ptsDerrapeY = tanVel(ptsDerrapeX);
            plot(puntos_x,puntos_y);
            hold on
            plot(equismax-25,gradasmax-3,'b','LineWidth',3)
            plot(equismin-25,gradasmin+3,'b','LineWidth',3)
            plot(puntos_x(k),puntos_y(k),'or')
            title("Simulación")
            xlabel("Eje x")
            ylabel("Eje y")
            ylim([punto_min_y-10 punto_max_y+10])
            plot(ptsDerrapeX,ptsDerrapeY,"color","r");
            calor = 0.5*masa*(velMs2^2);
            fprintf("La distancia recorrida por el auto al derraparse antes de detenerse es de %0.3f m\n",s);
            fprintf("La energía total perdida debido al derrape es de %0.3f J\n",calor);
            break
        end
    end
    % crea un punto que simula movimiento
    plot(puntos_x,puntos_y)
    hold on
    plot(equismax-25,gradasmax-3,'b','LineWidth',3)
    plot(equismin-25,gradasmin+3,'b','LineWidth',3)
    plot(puntos_x(k),puntos_y(k),'or')
    title("Simulación")
    xlabel("Eje x")
    ylabel("Eje y")
    ylim([punto_min_y-10 punto_max_y+10])
    pause(1/10000)    
    if k ~= length(puntos_x)-400
        clf
    end
end


%% Graficar
% figure(1)
% subplot(1,2,1)
% plot(puntos_x,puntos_y)
% hold on
% plot(x1,y1,'-or');
% plot(x2,y2,'-or');
% plot(x3,y3,'-or');
% plot(x4,y4,'-or');
% plot(punto_inflec_x,punto_inflec_y,'-_k')
% title("Figura 1 - Zona de curvas")
% xlabel("Eje x")
% ylabel("Eje y")
% ylim([punto_min_y-10 punto_max_y+10])
% 
% subplot(1,2,2)
% plot(puntos_x,puntos_y)
% hold on
% plot(puntos_x(ind_pmaxy),punto_max_y,'-or')
% plot(puntos_x(ind_pminy),punto_min_y,'-or')
% title("Figura 2 - Puntos máximo y mínimo")
% xlabel("Eje x")
% ylabel("Eje y")
% ylim([punto_min_y-10 punto_max_y+10])
% 
% 
% figure(2)
% subplot(1,2,1)
% plot(puntos_x,puntos_y)
% hold on
% plot(prxmax,prymax,'xr')
% plot(prxmax,prymax,'xr')
% plot(prxmin,prymin,'xr')
% plot(prxmin,prymin,'xr')
% title("Figura 3 - Puntos de derrape")
% xlabel("Eje x")
% ylabel("Eje y")
% ylim([punto_min_y-10 punto_max_y+10])
% 
% subplot(1,2,2)
% plot(puntos_x,puntos_y)
% hold on
% title("Figura 4 - Rectas tangentes")
% xlabel("Eje x")
% ylabel("Eje y")
% ylim([punto_min_y-20 punto_max_y+20])
% 
% for j = 1:length(prxmax)+1 %utilizamos ciclo for para obtener varias rectas tangentes
%     if j==length(prxmax)
%         break
%     end                    %el 20 determina la longitud de la recta tangente
%     xmax = prxmax(j):prxmax(j)+20; %modificar longitud de rectas en punto max   
%     atxmax = prxmax(j); %definir un punto x en max (cuando x=número)
%     indmax = find(puntos_x==atxmax); %encontrar indice para cuando x=número en punto max
%     bmax = puntos_y(indmax)-m1(indmax)*puntos_x(indmax); %encontrar b de y=mx+b en punto max
%     yrectmax = m1(indmax).*xmax+bmax; %evaluar puntos cuando x=número usando y=mx+b en punto max
%     plot(xmax,yrectmax,'r') %graficar rectas tangentes en punto max
% 
%     xmin = prxmin(j):prxmin(j)+20; %modificar longitud de rectas en punto min
%     atxmin = prxmin(j); %definir un punto x en min (cuando x=número)
%     indmin = find(puntos_x==atxmin); %encontrar indice para cuando x=número en punto min
%     bmin = puntos_y(indmin)-m1(indmin)*puntos_x(indmin); %encontrar b de y=mx+b en punto min
%     yrectmin = m1(indmin).*xmin+bmin; %evaluar puntos cuando x=número usando y=mx+b en punto min
%     plot(xmin,yrectmin,'r') %graficar rectas tangentes en punto min
% end
% 
% 
% figure(3)
% función
% plot(puntos_x,puntos_y)
% hold on
% 
% % puntos
% plot(x1,y1,'-ok');
% plot(x2,y2,'-ok');
% plot(x3,y3,'-ok');
% plot(x4,y4,'-ok');
% plot(punto_inflec_x,punto_inflec_y,'-_k')
% 
% % punto máximo y mínimo
% plot(puntos_x(ind_pmaxy),punto_max_y,'*m')
% plot(puntos_x(ind_pminy),punto_min_y,'*m')
% 
% % puntos de derrape
% plot(prxmax,prymax,'xr')
% plot(prxmax,prymax,'xr')
% plot(prxmin,prymin,'xr')
% plot(prxmin,prymin,'xr')
% 
% % rectas tangentes
% for j = 1:length(prxmax)+1
%     if j==length(prxmax)
%         break
%     end                    %el 20 determina la longitud de la recta tangente
%     xmax = prxmax(j):prxmax(j)+20; %modificar longitud de rectas en punto max   
%     atxmax = prxmax(j); %definir un punto x en max (cuando x=número)
%     indmax = find(puntos_x==atxmax); %encontrar indice para cuando x=número en punto max
%     bmax = puntos_y(indmax)-m1(indmax)*puntos_x(indmax); %encontrar b de y=mx+b en punto max
%     yrectmax = m1(indmax).*xmax+bmax; %evaluar puntos cuando x=número usando y=mx+b en punto max
%     plot(xmax,yrectmax,':r') %graficar rectas tangentes en punto max
% 
%     xmin = prxmin(j):prxmin(j)+20; %modificar longitud de rectas en punto min
%     atxmin = prxmin(j); %definir un punto x en min (cuando x=número)
%     indmin = find(puntos_x==atxmin); %encontrar indice para cuando x=número en punto min
%     bmin = puntos_y(indmin)-m1(indmin)*puntos_x(indmin); %encontrar b de y=mx+b en punto min
%     yrectmin = m1(indmin).*xmin+bmin; %evaluar puntos cuando x=número usando y=mx+b en punto min
%     plot(xmin,yrectmin,':r') %graficar rectas tangentes en punto min
%     
%     % para las gradas:
%     if atxmax==median(prxmax)-3.5
%         pendmax = m1(indmax);
%         pendmin = m1(indmin);
%         equismax = xmax;
%         equismin = xmin;    %hacemos copias de las rectas en un punto especifico para despues
%         b_max = bmax;       %simplemente dezplazar la linea a un lado
%         b_min = bmin;
%         gradasmax = pendmax.*equismax+b_max;
%         gradasmin = pendmin.*equismin+b_min;
%     end
% end
% 
% % gradas (el -25 dezplaza las rectas tangentes a la izquierda, -3 hacia abajo, +3 hacia arriba)
% plot(equismax-25,gradasmax-3,'b','LineWidth',3)
% plot(equismin-25,gradasmin+3,'b','LineWidth',3)
% 
% title("Figura 6 - Plano")
% xlabel("Eje x")
% ylabel("Eje y")
% ylim([punto_min_y-20 punto_max_y+20])
% grid on



%% Mostrar coeficientes, longitud de curva, radio de curvatura
indi = 3;
for j = 1:length(x)
    fprintf("Coeficiente a%d: %3.4f\n",indi,x(j))
    indi = indi - 1;
end

fprintf("Longitud de curva: 300 u < %4.2f u < 500 u\nRadio de curvatura: %2.4f\n",L,R)

%fprintf("Radio punto máximo: %2.4f\nRadio punto mínimo: %2.4f\n", Rmax, Rmin)

%fprintf("Indice punto max: %d\nIndice punto min: %d\nPunto max: (%3.4f,%3.4f)\nPunto min: (%3.4f,%3.4f)\n",ind_pmaxy,ind_pminy,punto_max_x,punto_max_y,punto_min_x,punto_min_y)






