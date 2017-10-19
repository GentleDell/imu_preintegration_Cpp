importdata('IMUcam_RawData.txt');

% gyro
ans.data(:,16) = ans.data(:,16);
ans.data(:,17) = -ans.data(:,17);
ans.data(:,18) = -ans.data(:,18);

% acc
ans.data(:,28) = ans.data(:,28);
ans.data(:,29) = -ans.data(:,29);
ans.data(:,30) = -ans.data(:,30);

fid = fopen('testdata1.txt', 'wt');
fprintf(fid,[repmat('%d\t', 1, size(ans.data,2)), '\n'], ans.data');
fclose(fid);