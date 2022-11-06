function datout = convertIntanToDat(data, name, time, samplerate)

ichanlabel = name;
ichanval = [];
itstart = time(1);
itend = time(end);

if contains(name,'DIN')

        iunits = 's';
        
        iData_rise = time(diff(data)>0);
        iData_fall = time(diff(data)<0);
        
        iData_rise = iData_rise(:);
        iData_fall = iData_fall(:);
        
        % First column is rise times, second column is fall times
        try
            
            if isempty(iData_rise)
                iData_rise = time(1);
            end
            if isempty(iData_fall)
                iData_fall = time(end);
            end
            if iData_rise(1)>iData_fall(1)
                iData_rise = [NaN; iData_rise];
            end
            if iData_fall(end)<iData_rise(end)
                iData_fall = [iData_fall; NaN];
            end
        catch
            warning('Error in digital process')
            
        end        
        datout = dat([iData_rise iData_fall], ichanlabel,ichanval,'dig',itstart,itend,iunits);     
        
else        
    if ~isempty(time) && ~exist('samplerate','var')
        samplerate = 1/mean(diff(time));
    end
        iunits = 'uV';
%         data = single(data) *0.195; % Convert to uV;% Changed 6/1/2019 -- convert outside of this function
        datout = dat(data,ichanlabel,ichanval, samplerate,itstart,itend,iunits);                
end
