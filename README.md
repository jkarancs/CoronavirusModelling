# CoronavirusModelling

Data-driven coronavirus epidemiology modelling based on simple projections.

### &#x1F539; Data:

The code uses the up to date daily comfirmed case and death data from ECDC:
https://www.ecdc.europa.eu/en/geographical-distribution-2019-ncov-cases

It downlaods this csv:
https://opendata.ecdc.europa.eu/covid19/casedistribution/csv

It also uses the latest daily testing data:
https://ourworldindata.org/covid-testing

It downloads this csv:
https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/testing/covid-testing-all-observations.csv

### &#x1F539; setup ROOT if you don't yet have it:

cd ~
git clone http://root.cern.ch/git/root.git
cd root
git tag -l
git checkout -b v5-34-08 v5-34-08
./configure
make -j 8
. bin/thisroot.sh

### &#x1F539; Recipe to run it:

First, set up ROOT environments above and then run the python script:

git clone https://github.com/jkarancs/CoronavirusModelling.git
cd CoronavirusModelling
python corona.py

The code does all the rest (download/update csv files, creat plots etc.)

See the result png files that are generated.
