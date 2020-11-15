function ldvData = importLDVdata_SinglePointAVG(filename)
filenameFull = sprintf("%s.STATISTICS.MSEBP.FILTERED.txt",filename);

data = importdata(filenameFull);

ldvData.meanVel = data(1);
ldvData.rms = data(2);
ldvData.snr = data(3);
ldvData.snrrms = data(4);
ldvData.datarate = data(5);
ldvData.nsamples = data(6);
ldvData.duration = data(7);
ldvData.internalLimit = data(8);




end