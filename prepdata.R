metafname = "~/data/weather/ground_measurements/nrel_hawaii/FieldDefinitions.txt"
meta = read.csv(metafname, sep = "\t")

datapath = "~/data/weather/ground_measurements/nrel_hawaii/"
fname = "20100318.txt"
singledata = read.csv(file.path(datapath, fname), col.names = meta$Instrument, na = '-99999')

files = list.files(path = "~/data/weather/ground_measurements/nrel_hawaii/", pattern="20*.txt")
readdata = function(fn){
  read.csv(file.path(datapath, fn), col.names = meta$Instrument, na = '-99999')
}
data = do.call(rbind, lapply(files[1:2], readdata))
