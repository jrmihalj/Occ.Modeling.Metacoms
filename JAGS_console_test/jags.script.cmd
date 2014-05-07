#jags.test file for simulating metacommunities

model in "OccMod_SingleYear.txt"
data in "data.R"
compile, nchains(3)
parameters in "inits.R"
initialize
adapt 1000
update 5000
monitor p, thin(10)
monitor z, thin(10)
update 10000
coda *
exit