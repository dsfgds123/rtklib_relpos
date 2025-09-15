function obs = readrinexobs3(file)
% readrinexobs3: reads RINEX 3.02 observation file for GPS & BDS
% Version 3: Outputs a structure that includes both GPST (sow) and the full time vector.

obs = [];
obs_count = 0;

fid = fopen(file);
if fid == -1, error('Cannot open file: %s', file); end

% --- Header Reading ---
obs_types = struct();
header_ended = false;
while ~header_ended
    line = fgetl(fid);
    if isempty(line), break; end
    if contains(line, 'SYS / # / OBS TYPES')
        sys_char = line(1);
        if ~isfield(obs_types, sys_char)
            num_types = str2double(line(4:6));
            types = strsplit(strtrim(line(7:60)));
            for i = 1:ceil(num_types/13)-1
                line_cont = fgetl(fid);
                types = [types, strsplit(strtrim(line_cont(7:60)))];
            end
            obs_types.(sys_char) = types(1:num_types);
        end
    elseif contains(line, 'END OF HEADER')
        header_ended = true;
    end
end

% --- Data Body Reading ---
current_epoch_data = [];
while ~feof(fid)
    line = fgetl(fid);
    if isempty(line) || ~ischar(line), continue; end
    if startsWith(line, '>')
        if ~isempty(current_epoch_data)
            obs_count = obs_count + 1;
            obs(obs_count).data = current_epoch_data;
            obs(obs_count).time = current_sow;
            obs(obs_count).time_vec = current_time_vec; % <-- 新增：保存完整时间向量
        end
        
        current_epoch_data = [];
        time_vec_in = sscanf(line(3:21), '%f');
        current_time_vec = [time_vec_in(1), time_vec_in(2), time_vec_in(3), time_vec_in(4), time_vec_in(5), time_vec_in(6)];
        
        jd = juliandate(datetime(current_time_vec));
        [~, sow] = jd2gpst(jd);
        current_sow = sow;
    else
        sys_char = line(1);
        prn_str = line(2:3);
        [sys, prn] = satsys(prn_str, sys_char);
        if sys == 0 || ~isfield(obs_types, sys_char), continue; end
        
        sat_id = prn;
        current_sys_types = obs_types.(sys_char);
        obs_vals = nan(1, 4);
        data_str = line(4:end);
        
        for i = 1:length(current_sys_types)
            start_pos = (i-1)*16 + 1;
            end_pos = i*16;
            if end_pos > length(data_str), break; end
            field = data_str(start_pos:end_pos);
            val = str2double(field(1:14));
            if i <= 4
                obs_vals(i) = val;
            end
        end
        
        new_row = [sat_id, obs_vals(1), obs_vals(2), 0, 0, 0, 0, sys];
        current_epoch_data = [current_epoch_data; new_row];
    end
end
if ~isempty(current_epoch_data)
    obs_count = obs_count + 1;
    obs(obs_count).data = current_epoch_data;
    obs(obs_count).time = current_sow;
    obs(obs_count).time_vec = current_time_vec; % <-- 新增：保存最后一个历元的时间向量
end
fclose(fid);
end

function [week, sow] = jd2gpst(jd)
    sow = mod(jd - 2444244.5, 7) * 86400;
    week = floor((jd - 2444244.5) / 7);
end