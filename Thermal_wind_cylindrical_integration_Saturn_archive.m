clearvars; close all; format long e; 

%%%%%%%% parameters of Saturn %%%%%%%%
R_saturn = 3760.051;  Cp_saturn = 10570.18; %%%%%% gas constants; 
Req_saturn = 60268*1000; Rpo_saturn = 54364*1000; %%%%%% equatorial and polar radii; 
g_1bar_saturn = 11.19; %%%%%%%  gravity at 1-bar surface (global average). 
day_saturn = 10+34/60.0+13/3600.0;%%%%%%%% rotational period; 
Omega_saturn = 2*pi/(day_saturn*3600.0); %%%%%%%%% angular velocity;  

%%%%%%%%% read temperature field and other related parameters %%%%%%%
%%%%%%%%% such as dz_lev_lat, g_lev_lat, and r_lev_lat; %%%%%%%%                   
read_dir =  '/Users/Liming/Work/Laptop_2021/Work/VIMS_winds/Study_2023/Fletcher_T_data_2023/'; 
save_name = strcat(read_dir,'Saturn_Fletcher_2023_latc181_lev98_up_parameters_New_A_500mbar.mat');
load(save_name);       
lat = latc_181;                                          

%%%%%%%% Using the new equations in "Grid_conversion_2021_Nov_10_B.pdf" to
%%%%%%%% get the new latitudes at different levels when integrating downward 
%%%%%%%% (see eq.(2) in my note); 
lat_old_lev_lat = ones(length(P_lev),length(lat))*-999.9; 
for i = 1:length(P_lev) 
lat_old_lev_lat(i,:) = lat; 
end 

coslat_lev_lat = cos(lat_old_lev_lat*pi/180.0); 
lat_new_lev_lat = lat_old_lev_lat; 

%%%%% for NH; 
lat_new_lev_lat(:,92:181) = acos((r_lev_lat(:,92:181)./(r_lev_lat(:,92:181)+dz_lev_lat(:,92:181)))...
                                 .*coslat_lev_lat(:,92:181))*180.0/pi; 
%%%%% for SH; 
lat_new_lev_lat(:,1:90) = acos((r_lev_lat(:,1:90)./(r_lev_lat(:,1:90)+dz_lev_lat(:,1:90)))...
                                 .*coslat_lev_lat(:,1:90))*180.0/pi*(-1); %%% for SH; 

%%%%%%%% set 2004_09 time-mean temperature as the 15th time for thermal winds. 
tmp_T_lev_lat = NaN(length(P_lev),length(latc_181));
for i = 1:length(P_lev) 
    for j = 1:length(latc_181)
        tmp_N = 0; tmp_data = 0.0; 
        for k = 1:6
            if T_yr_lev_lat(k,i,j) > 0 
                tmp_N = tmp_N + 1; 
                tmp_data = tmp_data+T_yr_lev_lat(k,i,j); 
            end 
        end 
        if tmp_N > 0
            tmp_T_lev_lat(i,j) = tmp_data/tmp_N; 
        end 
    end 
end 
T_yr_lev_lat(15,:,:) = tmp_T_lev_lat; 



%%%%%%%% Compute the meridional temperature gradient dTdy_lev_lat for each year

[n1,n2,n3] = size(T_yr_lev_lat); 
dTdy_yr_lev_lat = ones(n1,n2,n3)*-999.9; 
 
%%%%%% Equatorial and polar radii functioning with height. 
Req_saturn_r = Req_saturn + cumsum(dz_lev_lat(:,90)); %%%% for equator;   
Rpo_saturn_r = Rpo_saturn + cumsum(dz_lev_lat(:,1)); %%%% use 70N for pole; 

%%%%%% radius functioning with latitude and height. 
r_lev_lat = zeros(length(P_lev),length(lat));  
for i = 1:length(P_lev);
    for j = 1:length(lat); 
        lat_rad = lat(j,1)*pi/180.0;  
        tmp_1 = (Rpo_saturn_r(i,1)*cos(lat_rad)).^2; 
        tmp_2 = (Req_saturn_r(i,1)*sin(lat_rad)).^2; 
        r_lev_lat(i,j) = (Req_saturn_r(i,1)*Rpo_saturn_r(i,1))./((tmp_1+tmp_2).^0.5); 
    end 
end 

delta_2grid_lat = (lat(3,1)-lat(1,1))*pi/180*r_lev_lat;
delta_grid_lat = (lat(2,1) - lat(1,1))*pi/180*r_lev_lat;   

for i = 1:n1
    for j = 1:n2
        for k = 2:(n3-1)
            dTdy_yr_lev_lat(i,j,k) = (T_yr_lev_lat(i,j,k+1)-T_yr_lev_lat(i,j,k-1))/delta_2grid_lat(j,k); 
        end 
        % note: latitude is from -89.5 to 89.5. 
    end
end 
       
for i = 1:n1
    tmp_a = squeeze(T_yr_lev_lat(i,:,2)); tmp_a = tmp_a'; 
    tmp_b = squeeze(T_yr_lev_lat(i,:,1)); tmp_b = tmp_b'; 
    dTdy_yr_lev_lat(i,:,1) = (tmp_a-tmp_b)./delta_grid_lat(:,1);
    tmp_a = squeeze(T_yr_lev_lat(i,:,n3)); tmp_a = tmp_a'; 
    tmp_b = squeeze(T_yr_lev_lat(i,:,n3-1)); tmp_b = tmp_b'; 
    dTdy_yr_lev_lat(i,:,n3) = (tmp_a-tmp_b)./delta_grid_lat(:,end);
end

dTdy_yr_lev_lat_mid = ones(n1,n2,n3)*-999.9; 
T_yr_lev_lat_mid = ones(n1,n2,n3)*-999.9; 
for ii = 1:(n2-1)
    dTdy_yr_lev_lat_mid(:,ii,:) = (dTdy_yr_lev_lat(:,ii,:)+dTdy_yr_lev_lat(:,ii+1,:))/2.0;
    T_yr_lev_lat_mid(:,ii,:) = (T_yr_lev_lat(:,ii,:)+T_yr_lev_lat(:,ii+1,:))/2.0; 
end 
dTdy_yr_lev_lat_mid(:,n2,:) = dTdy_yr_lev_lat_mid(:,n2-1,:); 
T_yr_lev_lat_mid(:,n2,:)= T_yr_lev_lat_mid(:,n2-1,:);

dTdy_yr_lev_lat_raw = dTdy_yr_lev_lat; 

%%%%% Fill missing data in temperature; 
T_yr_lev_lat = fillmissing(T_yr_lev_lat,'linear',1,'EndValues','nearest'); 

%%%%%% Fill missing data in temepture gradient; 
dTdy_yr_lev_lat = fillmissing(dTdy_yr_lev_lat,'linear',1,'EndValues','nearest'); 

%%%%%%%% using eq.(4) in my note to get dH from dR (dz = dr = dR);  
sinlat_old_lev_lat = sind(lat_old_lev_lat); 
sinlat_new_lev_lat = sind(lat_new_lev_lat); 

dH_lev_lat = (r_lev_lat + dz_lev_lat).*sinlat_new_lev_lat ...
              -r_lev_lat.*sinlat_old_lev_lat; 
   
dH_lev_lat(:,1:90) = dH_lev_lat(:,1:90)*(-1); %%%% set both NH and SH positive;          
dH_lev_lat(1,:) = 0.0; %%%%% there is a very small value for level 1 
%%%%%%%% comes from the very small difference between lat_new and lat_old
%%%%%%%% (further associated to the acos in line 188. 

%%%%%%%%% Input the ISS CB winds as the boudnary condition at 500 mbar; 
dir_a = '/Users/Liming/Work/Laptop_2021/Work/Grand_Finale/'; 
tmp_a =strcat(dir_a, 'ISS_CB23_2004_08_global_wind_A.mat'); 
load(tmp_a);

CB23_U_saturn = interp1(latc_2004_08_glo,wind_2004_08_glo,lat,'linear','extrap');
U_lev1 = CB23_U_saturn; %%%%%% set Cassini winds at level 1 around 0.5 bar; 
% 
U_lev_lat = ones(length(P_lev),length(lat))*-999.9; 
U_lev_lat(1,:) = U_lev1'; 

[n1,n2,n3] = size(dTdy_yr_lev_lat);
U_yr_lev_lat = ones(n1,n2,n3)*-999.9;
for i10 = 1:n1
    dTdy_lev_lat = squeeze(dTdy_yr_lev_lat(i10,:,:));
    T_lev_lat = squeeze(T_yr_lev_lat(i10,:,:));

    %%%%% NH;
    for i = 1:length(P_lev)-1
        tmp_x = lat_old_lev_lat(i+1,92:181);
        tmp_y = dTdy_lev_lat(i+1,92:181);
        tmp_x_N = lat_new_lev_lat(i+1,92:181);
        tmp_y_N = interp1(tmp_x,tmp_y,tmp_x_N,'linear','extrap');
        dTdy_mid = dTdy_lev_lat(i,92:181)*1/2 + tmp_y_N*1/2;
        T_mid = T_lev_lat(i,92:181)*1/2 + T_lev_lat(i+1,92:181)*1/2;
        g_mid = g_lev_lat(i,92:181)*1/2 + g_lev_lat(i+1,92:181)*1/2;
        %%%%% for above T_mid and g_mid, we may follow dT_dy with linear
        %%%%% interpolation (but no big difference);
        
        %%%%% for NH, we have dH as positive for upward (dH from above is positive globally);
        tmp_fac = g_mid./T_mid;
         
        for ii = 92:181
            tmp_fac_ii = tmp_fac(ii-91);
            tmp_fac_ii_A = (-1)*dH_lev_lat(i+1,ii)*tmp_fac_ii*dTdy_mid(ii-91);
            tmp_U_square = U_lev_lat(i,ii)^2/(r_lev_lat(i,ii)*coslat_lev_lat(i,ii))+2*U_lev_lat(i,ii)*Omega_saturn;
            tmp_U_square_fac = tmp_U_square + tmp_fac_ii_A;
            
            tmp_key = ones(1,10000)*0.01;
            tmp_key = cumsum(tmp_key);
            if tmp_fac_ii_A < 0
                tmp_key = -1*tmp_key;
            end
            
            tmp_U_square_N = (U_lev_lat(i,ii)+tmp_key).^2/(r_lev_lat(i+1,ii)*coslat_lev_lat(i+1,ii))+2*(U_lev_lat(i,ii)+tmp_key)*Omega_saturn;
            tmp_dif = abs(tmp_U_square_N-tmp_U_square_fac);
            tmp_min = min(tmp_dif);
            tmp_delta_U = tmp_key(1,tmp_dif == tmp_min);
            
            %%%%% Note: due to precision, sometimes the above sentences
            %%%%% generates multiple values (they are close to each other);
            if length(tmp_delta_U)> 1
                tmp_delta_U = mean(tmp_delta_U);
            end
            U_lev_lat(i+1,ii) = U_lev_lat(i,ii) + tmp_delta_U;
        end
        
        %%%%% interpolate computed winds at new latitudes to original/regular latitude grid for next level;
        tmp_U_N = U_lev_lat(i+1,92:181);
        tmp_U = interp1(tmp_x_N,tmp_U_N,tmp_x,'linear','extrap');
        U_lev_lat(i+1,92:181) = tmp_U;
        %end
    end
    %%%%% SH;
    for i = 1:length(P_lev)-1
        %for j = 1:900
        tmp_x = lat_old_lev_lat(i+1,1:90);
        tmp_y = dTdy_lev_lat(i+1,1:90);
        tmp_x_N = lat_new_lev_lat(i+1,1:90);
        tmp_y_N = interp1(tmp_x,tmp_y,tmp_x_N,'linear','extrap');
        dTdy_mid = dTdy_lev_lat(i,1:90)*1/2 + tmp_y_N*1/2;
        T_mid = T_lev_lat(i,1:90)*1/2 + T_lev_lat(i+1,1:90)*1/2;
        g_mid = g_lev_lat(i,1:90)*1/2 + g_lev_lat(i+1,1:90)*1/2;
        %%%%% for above T_mid and g_mid, we may follow dT_dy with linear
        %%%%% interpolation (but no big difference);
        
        %%%%% for SH, we have dH as negative for upward. So we have the factor "-1".
        %%%%% because dH from above is always positive in two hemispheres.
        %%%%% Note: in the classical TWE, dz is postive for both hemispheres
        %%%%% when upward, but the the sign of coriolis parameter (f) changes
        %%%%% from NH to SH. Here, dH changes sign from NH to SH because there
        %%%%% is no f parameter in the new equation.
        tmp_fac = (-1)*g_mid./T_mid;
        %U_lev_lat(i+1,1:900) = U_lev_lat(i,1:900)+ dH_lev_lat(i+1,1:900).*tmp_fac.*dTdy_mid;
        
        for ii = 1:90
            tmp_fac_ii = tmp_fac(ii);
            tmp_fac_ii_A = (-1)*dH_lev_lat(i+1,ii)*tmp_fac_ii*dTdy_mid(ii);
            tmp_U_square = U_lev_lat(i,ii)^2/(r_lev_lat(i,ii)*coslat_lev_lat(i,ii))+2*U_lev_lat(i,ii)*Omega_saturn;
            tmp_U_square_fac = tmp_U_square + tmp_fac_ii_A;
            
            tmp_key = ones(1,10000)*0.01;
            tmp_key = cumsum(tmp_key);
            if tmp_fac_ii_A < 0
                tmp_key = -1*tmp_key;
            end
            
            tmp_U_square_N = (U_lev_lat(i,ii)+tmp_key).^2/(r_lev_lat(i+1,ii)*coslat_lev_lat(i+1,ii))+2*(U_lev_lat(i,ii)+tmp_key)*Omega_saturn;
            tmp_dif = abs(tmp_U_square_N-tmp_U_square_fac);
            tmp_min = min(tmp_dif);
            tmp_delta_U = tmp_key(1,tmp_dif == tmp_min);
            %%%%% Note: due to precision, sometimes the above sentences
            %%%%% generates multiple values (they are close to each other);
            if length(tmp_delta_U)> 1
                tmp_delta_U = mean(tmp_delta_U);
            end
            
            U_lev_lat(i+1,ii) = U_lev_lat(i,ii) + tmp_delta_U;
        end
       
        %%%%% interpolate computed winds at new latitudes to original/regular latitude grid for next level;
        tmp_U_N = U_lev_lat(i+1,1:90);
        tmp_U = interp1(tmp_x_N,tmp_U_N,tmp_x,'linear','extrap');
        U_lev_lat(i+1,1:90) = tmp_U;

    end
    
    U_lev_lat(2:end,91) = U_lev_lat(2:end,90)*1/2+U_lev_lat(2:end,92)*1/2;
    
    U_yr_lev_lat(i10,:,:) = U_lev_lat;
end


%%%%% Testing high spatial resolution in latiude; 
lat_high = -90:0.1:90; 
lat_old_lev_lat_high = ones(length(P_lev),length(lat_high))*-999.9; 
lat_new_lev_lat_high = ones(length(P_lev),length(lat_high))*-999.9; 
dTdy_yr_lev_lat_high = ones(n1,length(P_lev),length(lat_high))*-999.9; 
T_yr_lev_lat_high = ones(n1,length(P_lev),length(lat_high))*-999.9; 
g_lev_lat_high = ones(length(P_lev),length(lat_high))*-999.9; 
U_lev_lat_high = ones(length(P_lev),length(lat_high))*-999.9; 
dH_lev_lat_high = ones(length(P_lev),length(lat_high))*-999.9; 
r_lev_lat_high = ones(length(P_lev),length(lat_high))*-999.9; 
coslat_lev_lat_high = ones(length(P_lev),length(lat_high))*-999.9; 

for i = 1:length(P_lev) 
    tmp_a = lat; tmp_b = U_lev_lat(i,:); 
    tmp_bb = interp1(tmp_a,tmp_b,lat_high,'linear','extrap'); 
    U_lev_lat_high(i,:) = tmp_bb; 
    
    tmp_a = lat; tmp_b = g_lev_lat(i,:); 
    tmp_bb = interp1(tmp_a,tmp_b,lat_high,'linear','extrap'); 
    g_lev_lat_high(i,:) = tmp_bb;  
    
    tmp_a = lat; tmp_b = dH_lev_lat(i,:); 
    tmp_bb = interp1(tmp_a,tmp_b,lat_high,'linear','extrap'); 
    dH_lev_lat_high(i,:) = tmp_bb;
    
    tmp_a = lat; tmp_b = lat_old_lev_lat(i,:); 
    tmp_bb = interp1(tmp_a,tmp_b,lat_high,'linear','extrap'); 
    lat_old_lev_lat_high(i,:) = tmp_bb; 
    
    tmp_a = lat; tmp_b = lat_new_lev_lat(i,:); 
    tmp_bb = interp1(tmp_a,tmp_b,lat_high,'linear','extrap'); 
    lat_new_lev_lat_high(i,:) = tmp_bb; 
    
    tmp_a = lat; tmp_b = r_lev_lat(i,:); 
    tmp_bb = interp1(tmp_a,tmp_b,lat_high,'linear','extrap'); 
    r_lev_lat_high(i,:) = tmp_bb; 
    
    tmp_a = lat; tmp_b = coslat_lev_lat(i,:); 
    tmp_bb = interp1(tmp_a,tmp_b,lat_high,'linear','extrap'); 
    coslat_lev_lat_high(i,:) = tmp_bb; 
    
end 

for i = 1:n1
for j = 1:length(P_lev)
    tmp_a = lat; tmp_b = squeeze(T_yr_lev_lat(i,j,:)); 
    tmp_bb = interp1(tmp_a,tmp_b,lat_high,'linear','extrap'); 
    T_yr_lev_lat_high(i,j,:) = tmp_bb; 
    
    tmp_a = lat; tmp_b = squeeze(dTdy_yr_lev_lat(i,j,:)); 
    tmp_bb = interp1(tmp_a,tmp_b,lat_high,'linear','extrap'); 
    dTdy_yr_lev_lat_high(i,j,:) = tmp_bb;    
end 
end 

lat_old_lev_lat = lat_old_lev_lat_high; 
lat_new_lev_lat = lat_new_lev_lat_high; 
dTdy_yr_lev_lat = dTdy_yr_lev_lat_high; 
T_yr_lev_lat = T_yr_lev_lat_high; 
g_lev_lat = g_lev_lat_high; 
U_lev_lat = U_lev_lat_high; 
dH_lev_lat = dH_lev_lat_high; 
r_lev_lat = r_lev_lat_high; 
coslat_lev_lat = coslat_lev_lat_high; 


lat = lat_high; 

[n1,n2,n3] = size(dTdy_yr_lev_lat);
U_yr_lev_lat_N = ones(n1,n2,n3)*-999.9;
for i11 = 1:n1
    dTdy_lev_lat = squeeze(dTdy_yr_lev_lat(i11,:,:));
    T_lev_lat = squeeze(T_yr_lev_lat(i11,:,:));    
    
    %%%%% NH;
    for i = 1:length(P_lev)-1
        tmp_x = lat_old_lev_lat(i+1,902:1801);
        tmp_y = dTdy_lev_lat(i+1,902:1801);
        tmp_x_N = lat_new_lev_lat(i+1,902:1801);
        tmp_y_N = interp1(tmp_x,tmp_y,tmp_x_N,'linear','extrap');
        dTdy_mid = dTdy_lev_lat(i,902:1801)*1/2 + tmp_y_N*1/2;
        T_mid = T_lev_lat(i,902:1801)*1/2 + T_lev_lat(i+1,902:1801)*1/2;
        g_mid = g_lev_lat(i,902:1801)*1/2 + g_lev_lat(i+1,902:1801)*1/2;
        %%%%% for above T_mid and g_mid, we may follow dT_dy with linear
        %%%%% interpolation (but no big difference);
        
        %%%%% for NH, we have dH as positive for upward (dH from above is positive globally);
        tmp_fac = g_mid./T_mid;
        %U_lev_lat(i+1,902:1801) = U_lev_lat(i,902:1801)+ (-1)*dH_lev_lat(i+1,902:1801).*tmp_fac.*dTdy_mid;
        
        for ii = 902:1801
            tmp_fac_ii = tmp_fac(ii-901);
            tmp_fac_ii_A = (-1)*dH_lev_lat(i+1,ii)*tmp_fac_ii*dTdy_mid(ii-901);
            tmp_U_square = U_lev_lat(i,ii)^2/(r_lev_lat(i,ii)*coslat_lev_lat(i,ii))+2*U_lev_lat(i,ii)*Omega_saturn;
            tmp_U_square_fac = tmp_U_square + tmp_fac_ii_A;
            
            tmp_key = ones(1,10000)*0.01;
            tmp_key = cumsum(tmp_key);
            if tmp_fac_ii_A < 0
                tmp_key = -1*tmp_key;
            end
            
            tmp_U_square_N = (U_lev_lat(i,ii)+tmp_key).^2/(r_lev_lat(i+1,ii)*coslat_lev_lat(i+1,ii))+2*(U_lev_lat(i,ii)+tmp_key)*Omega_saturn;
            tmp_dif = abs(tmp_U_square_N-tmp_U_square_fac);
            tmp_min = min(tmp_dif);
            tmp_delta_U = tmp_key(1,tmp_dif == tmp_min);
            
            %%%%% Note: due to precision, sometimes the above sentences
            %%%%% generates multiple values (they are close to each other);
            if length(tmp_delta_U)> 1
                tmp_delta_U = mean(tmp_delta_U);
            end
            U_lev_lat(i+1,ii) = U_lev_lat(i,ii) + tmp_delta_U;
        end
        
        %%%%% interpolate computed winds at new latitudes to original/regular latitude grid for next level;
        tmp_U_N = U_lev_lat(i+1,902:1801);
        tmp_U = interp1(tmp_x_N,tmp_U_N,tmp_x,'linear','extrap');
        U_lev_lat(i+1,902:1801) = tmp_U;
        
        %end
    end
    
    %%%%% SH;
    for i = 1:length(P_lev)-1
        tmp_x = lat_old_lev_lat(i+1,1:900);
        tmp_y = dTdy_lev_lat(i+1,1:900);
        tmp_x_N = lat_new_lev_lat(i+1,1:900);
        tmp_y_N = interp1(tmp_x,tmp_y,tmp_x_N,'linear','extrap');
        dTdy_mid = dTdy_lev_lat(i,1:900)*1/2 + tmp_y_N*1/2;
        T_mid = T_lev_lat(i,1:900)*1/2 + T_lev_lat(i+1,1:900)*1/2;
        g_mid = g_lev_lat(i,1:900)*1/2 + g_lev_lat(i+1,1:900)*1/2;
        %%%%% for above T_mid and g_mid, we may follow dT_dy with linear
        %%%%% interpolation (but no big difference);
        
        %%%%% for SH, we have dH as negative for upward. So we have the factor "-1".
        %%%%% because dH from above is always positive in two hemispheres.
        %%%%% Note: in the classical TWE, dz is postive for both hemispheres
        %%%%% when upward, but the the sign of coriolis parameter (f) changes
        %%%%% from NH to SH. Here, dH changes sign from NH to SH because there
        %%%%% is no f parameter in the new equation.
        tmp_fac = (-1)*g_mid./T_mid;
        %U_lev_lat(i+1,1:900) = U_lev_lat(i,1:900)+ dH_lev_lat(i+1,1:900).*tmp_fac.*dTdy_mid;
        
        for ii = 1:900
            tmp_fac_ii = tmp_fac(ii);
            tmp_fac_ii_A = (-1)*dH_lev_lat(i+1,ii)*tmp_fac_ii*dTdy_mid(ii);
            tmp_U_square = U_lev_lat(i,ii)^2/(r_lev_lat(i,ii)*coslat_lev_lat(i,ii))+2*U_lev_lat(i,ii)*Omega_saturn;
            tmp_U_square_fac = tmp_U_square + tmp_fac_ii_A;
            
            tmp_key = ones(1,10000)*0.01;
            tmp_key = cumsum(tmp_key);
            if tmp_fac_ii_A < 0
                tmp_key = -1*tmp_key;
            end
            
            tmp_U_square_N = (U_lev_lat(i,ii)+tmp_key).^2/(r_lev_lat(i+1,ii)*coslat_lev_lat(i+1,ii))+2*(U_lev_lat(i,ii)+tmp_key)*Omega_saturn;
            tmp_dif = abs(tmp_U_square_N-tmp_U_square_fac);
            tmp_min = min(tmp_dif);
            tmp_delta_U = tmp_key(1,tmp_dif == tmp_min);
            %%%%% Note: due to precision, sometimes the above sentences
            %%%%% generates multiple values (they are close to each other);
            if length(tmp_delta_U)> 1
                tmp_delta_U = mean(tmp_delta_U);
            end
            
            U_lev_lat(i+1,ii) = U_lev_lat(i,ii) + tmp_delta_U;
        end
        
        %%%%% interpolate computed winds at new latitudes to original/regular latitude grid for next level;
        tmp_U_N = U_lev_lat(i+1,1:900);
        tmp_U = interp1(tmp_x_N,tmp_U_N,tmp_x,'linear','extrap');
        U_lev_lat(i+1,1:900) = tmp_U;
        
        %end
    end
    
    
    U_lev_lat(2:end,901) = U_lev_lat(2:end,900)*1/2+U_lev_lat(2:end,902)*1/2;
    
    U_yr_lev_lat_N(i11,:,:) = U_lev_lat;
end

U_yr_lev_lat_lat1 = U_yr_lev_lat; 
U_yr_lev_lat_lat01 = U_yr_lev_lat_N; 

T_yr_lev_lat_fill = T_yr_lev_lat; 
dTdy_yr_lev_lat_fill = dTdy_yr_lev_lat; 

latc_1801 = lat; 

read_dir =  '/Users/Liming/Work/Laptop_2021/Work/VIMS_winds/Study_2023/Fletcher_T_data_2023/'; 
save_name = strcat(read_dir,'Saturn_Fletcher_2023_latc181_lev98_up_winds_A_500mbar_archive.mat');
save(save_name,'latc_181','latc_1801','P_lev','T_yr_lev_lat_raw','T_yr_lev_lat_fill',...
               'g_lev_lat','r_lev_lat','dz_lev_lat','dH_lev_lat','lat_old_lev_lat',...
               'lat_new_lev_lat','dTdy_yr_lev_lat_raw','dTdy_yr_lev_lat_fill',...
               'U_yr_lev_lat_lat1','U_yr_lev_lat_lat01'); 
