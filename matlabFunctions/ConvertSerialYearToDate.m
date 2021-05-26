function [num] = ConvertSerialYearToDate( y )
  year = floor(y);
  partialYear = mod(y,1);
  date0 = datenum(num2str(year),'yyyy');
  date1 = datenum(num2str(year+1),'yyyy');
  daysInYear = date1 - date0;
  if any(size(y)>[1 1])
    num = date0 + partialYear.*daysInYear;
  else
    num = date0 + partialYear*daysInYear;
  end      
end