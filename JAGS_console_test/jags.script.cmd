#jags.test file for simulating metacommunities

cd "/Users/JRMihaljevic/Documents/Thesis Research/Occ.Model.EMS/JAGS_console_test"
model in "OccMod_SingleYear.txt"
data in "data.R"
compile, nchains(3)
parameters in "inits.R"
initialize
adapt 1000
update 5000
monitor lpsiMean, thin(10)
monitor lpMean, thin(10)
monitor z, thin(10)
update 10000
coda *
exit