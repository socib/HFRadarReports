clear all
clc

link = 'http://thredds.socib.es/thredds/dodsC/hf_radar/hf_radar_ibiza-scb_codarssproc001/L1/2016/dep0001_hf-radar-ibiza_scb-codarssproc001_L1_2016-04.nc'

time = ncread(link, 'time');
u = ncread(link, 'U');
v = ncread(link, 'V');
lat = ncread(link, 'LAT');
lon = ncread(link, 'LON');



desired_constituents = ['K1  '; 'M2  '; 'S2  ']
m_max = length(lat)
n_max = length(lon)

figure(1)
figure(2)
figure(3)
for m=1:m_max
    cur_lat = lat(m);
    for n=1:n_max
        cur_u = squeeze(u(n,m,:));
        cur_v = squeeze(v(n,m,:));
        check1 = length(find(isnan(cur_u)))/length(cur_u)*100;
        check2 = length(find(isnan(cur_v)))/length(cur_v)*100;
        if check1 >99 || check2 >99
            continue;
        end
        nan_idx_u = isnan(cur_u);
        nan_idx_v = isnan(cur_v);
        cur_u_interpolated = interp1(time(~nan_idx_u), cur_u(~nan_idx_u), time);
        cur_v_interpolated = interp1(time(~nan_idx_v), cur_v(~nan_idx_v), time);
        complex_u_v = cur_u_interpolated + 1i * cur_v_interpolated; 
        [const_names, const_freqs, tide_const, prediction] = t_tide(complex_u_v, 'interval', 1,  'start_time', datenum(posixtime2utc(time(1))), 'latitude', cur_lat, 'output', 'none');
        idx1 = strmatch('K1  ', const_names);
        idx2 = strmatch('M2  ', const_names);
        idx3 = strmatch('S2  ', const_names);
        maj_e1 = tide_const(idx1, 1);
        maj_e2 = tide_const(idx2, 1);
        maj_e3 = tide_const(idx3, 1);
        min_e1 = tide_const(idx1, 3);
        min_e2 = tide_const(idx2, 3);
        min_e3 = tide_const(idx3, 3);
        angle_e1 = tide_const(idx1, 5);
        angle_e2 = tide_const(idx2, 5);
        angle_e3 = tide_const(idx3, 5);
        figure(1)
        hold on
        ellipse(maj_e1,min_e1,deg2rad(angle_e1),lon(n),lat(m),'k');
        figure(2)
        hold on
        ellipse(maj_e2,min_e2,deg2rad(angle_e2),lon(n),lat(m),'k');
        figure(3)
        hold on
        ellipse(maj_e3,min_e3,deg2rad(angle_e3),lon(n),lat(m),'k');
    end
end
figure(1)
axis equal
figure(2)
axis equal
figure(3)
axis equal