import sys, math, time, operator, urllib2, os, stat, csv, ROOT
import numpy as np

trends = 0

extra = 4
ndayspread = 1

# Auto select worst 8 countries
# Criteria: worst fraction in 14 days, 10 mill+ pop, 1k+ current cases
worst  = 0
min_pop = 1e7
case_threshold = 1e3
days_pred = 14
# Auto select best 8 countries
# Criteria: best current rate, 1 mill+ pop, 1k+ current cases
best   = 0
best_min_pop = 1e6
best_case_threshold = 1e3
# set both best/worst to off to manually select countries

save = 1

# Plotting and linear projection settings
if best:
    # linear project dialy ratio back to n days
    ndaysback  = 7
    ndaysback2 = 4
    # plot options
    zoom = 0
    ndayspread = 0 if zoom else 2
    xrange  = [80,110] if zoom else [0,120]
    yrange  = [1.01e0,1e6]
    logy2  = 0
    yrange2 = [1.01e-2,1e2] if logy2 else ([0,10] if zoom else [0,80])
    yrange3  = [0,100]
    yrange4  = [1e2,1e8]
    y2leg = 0.5 if zoom else 0.95
elif worst:
    # linear project dialy ratio back to n days
    ndaysback  = 14
    ndaysback2 = 6
    # plot options
    xrange  = [30,120]
    yrange  = [1.01e0,1e7]
    logy2   = 0
    yrange2 = [1.01e-2,1e2] if logy2 else [0,100]
    yrange3  = [0,100]
    yrange4  = [1e2,1e9]
    y2leg = 0.95
else:
    # Manually select countries
    # linear project dialy ratio back to n days
    ndaysback  = 5
    ndaysback2 = 4
    # plot options
    xrange   = [60,120]
    yrange   = [1.01e0,1e7]
    yrange2  = [0,100]
    yrange3  = [0,60]
    yrange4  = [1e2,1e8]
    logy2  = 0
    y2leg = 0.95
    countries = [
        #"WRL",
        #"CHN",
        #"KOR",
        #"USA",
        #"FRA",
        ##"ESP",
        #"ITA",
        #"AUT",
        #"SWE",
        #"FIN",
        #"NOR",
        #"DEU",
        #"GBR",
        #"UKR",
        #"SRB",
        #"AUT",
        #"JPN",
        #"TUR",
        #"BEL",
        # Hungary and its neighbours
        "HRV",
        #"SVK",
        "AUT",
        "HUN",
        "SVN",
        "UKR",
        "ROU",
        "SRB",
        # RazorBoost
        #"CHE",
        #"TUR",
        #"USA",
        #"KOR",
        #"HUN",
        # CFR Spike
        #"FRA",
        #"HUN",
        #"ROU",
        #"BEL",
        #"IRL",
    ]

# plot settings: marker, hollow marker, color, size
settings = [
    #[20, 24,   1, 1.0],
    #[21, 25,   3, 1.0],
    #[22, 26,   2, 1.0],
    #[23, 32,   4, 1.0],
    #[33, 27,   6, 1.2],
    #[34, 28, 800, 1.2],
    #[29, 30, 804, 1.4],
    [22, 26, 601, 1.0, 3002],
    [21, 25, 632, 1.0, 3007],
    [20, 24,   1, 1.0, 3020],
    [23, 32, 417, 1.0, 3021],
    [33, 27, 617, 1.2, 3005],
    [34, 28, 433, 1.2, 3004],
    [29, 30, 401, 1.4, 3006],
    [20, 24, 801, 1.0, 3244],
]

# default plot style in ROOT
def set_default_style_():
    ROOT.gStyle.SetPaperSize(20.,20.)
    ROOT.gStyle.SetTitleFont(42,"xyz")
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(0)
    ROOT.gStyle.SetErrorX(0)
    ROOT.gStyle.SetFrameBorderMode(0)
    #ROOT.gStyle.SetFrameFillColor(0)
    #ROOT.gStyle.SetFrameFillStyle(0)
    ROOT.gStyle.SetFrameLineWidth(2)
    #ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetTitleFillColor(0)
    ROOT.gStyle.SetTitleStyle(0)
    ROOT.gStyle.SetTitleX(1)
    ROOT.gStyle.SetTitleY(1)
    ROOT.gStyle.SetTitleAlign(33)

# Methods
def daily_growth(h, dh):
    count = 0
    vdiff = []
    growth = 0
    for binx in range(1,h.GetNbinsX()+1):
        if h.GetBinContent(binx)>0:
            count += 1
            if count>extra:
                if h.GetMaximum()>20:
                    fit = ROOT.TF1("fit","expo",(binx-extra)*86400,(binx+extra)*86400)
                    h.Fit("fit","RMWQ0" if "WRL" in h.GetName() else "RMQ0")
                    vdiff.append((fit.Eval(h.GetBinCenter(binx))-h.GetBinContent(binx))/fit.Eval(h.GetBinCenter(binx)))
                    growth = (math.exp(fit.GetParameter(1)*86400)-1)*100
                    #growth = (h.GetBinContent(binx)-h.GetBinContent(binx-1))/max(h.GetBinContent(binx-1),1)*100
                dh. SetBinContent(binx,growth)
                dh. SetBinError  (binx,0)
        else:
            dh. SetBinContent(binx,0)
            dh. SetBinError  (binx,0)
    return np.mean(vdiff[-extra:]) if len(vdiff) else 0
    

def daily_growth_and_diff(h, dh, ddh, dddh):
    mean = daily_growth(h, dh)
    growth = 0
    for binx in range(1,h.GetNbinsX()+1):
        if dh.GetBinContent(binx)>0:
            prev_growth = dh.GetBinContent(binx-1)
            growth      = dh.GetBinContent(binx)
            #prev_growth = (dh.GetBinContent(binx-2)+dh.GetBinContent(binx+1)+dh.GetBinContent(binx))/3
            #growth      = (dh.GetBinContent(binx-1)+dh.GetBinContent(binx)+dh.GetBinContent(binx+1))/3
            growth_change = growth - prev_growth
            #growth_change = (growth - prev_growth)/prev_growth*100 if prev_growth!=0 else growth - prev_growth
            ddh.SetBinContent(binx,growth_change)
    for binx in range(1,h.GetNbinsX()+1):
        if ddh.GetBinContent(binx)>0:
            prev_growth = ddh.GetBinContent(binx-1)
            growth      = ddh.GetBinContent(binx)
            #prev_growth = (dh.GetBinContent(binx-2)+dh.GetBinContent(binx+1)+dh.GetBinContent(binx))/3
            #growth      = (dh.GetBinContent(binx-1)+dh.GetBinContent(binx)+dh.GetBinContent(binx+1))/3
            growth_change = growth - prev_growth
            #growth_change = (growth - prev_growth)/prev_growth*100 if prev_growth!=0 else growth - prev_growth
            dddh.SetBinContent(binx,growth_change)
            #if "ITA" in h.GetName():
            #    print ("bin=%d  cont=%d  fit=%.1f  growth=%.1f  change=%.1f " % (binx, h.GetBinContent(binx), fit.Eval(h.GetBinCenter(binx)),
            #          
    return mean

def get_pred(dh, h, linear_back, c1=-9999, draw=False):
    maxbin = 0
    maxcont = 0
    for binx in range(1,dh.GetNbinsX()+1):
        if dh.GetBinContent(binx)>0:
            maxbin = binx
            maxcont = h.GetBinContent(binx)
    if c1==-9999:
        fit = ROOT.TF1(dh.GetName()+"_fit_"+str(linear_back),"pol1",(maxbin-linear_back)*86400,(maxbin+28)*86400)
        name = h.GetName()+"_pred_"+str(linear_back)
        dh.Fit(dh.GetName()+"_fit_"+str(linear_back),"RMWQ0")
        if draw:
            fit.SetLineColor(dh.GetMarkerColor())
            fit.SetLineWidth(2)
            fit.Draw("SAME")
    else:
        fit = ROOT.TF1(dh.GetName()+"_curfew_fit_"+str(linear_back),"pol1",(maxbin-2)*86400,(maxbin+28)*86400)
        fit.FixParameter(1,c1)
        name = h.GetName()+"_pred_curfew_"+str(linear_back)
        dh.Fit(dh.GetName()+"_curfew_fit_"+str(linear_back),"RMQ0")
        if draw:
            fit.SetLineColor(dh.GetMarkerColor())
            fit.SetLineStyle(7)
            fit.SetLineWidth(2)
            fit.Draw("SAME")
    h_pred = h.Clone(name)
    h_pred.SetMarkerColor(h.GetMarkerColor())
    for binx in range(1,h_pred.GetNbinsX()+1):
        if h_pred.GetBinContent(binx)>0:
            h_pred.SetBinContent(binx,0)
        elif binx>maxbin:
            pred = maxcont*(1+(max(0, fit.Eval(h.GetBinCenter(binx)))/100))
            maxcont = pred
            h_pred.SetBinContent(binx, pred)
    if c1==-9999:
        return [h_pred, fit]
    else:
        return [h_pred, fit]

def get_pred_and_error(h, dh, diff):
    h_pred_nom = h.Clone(h.GetName()+"_pred_nom")
    h_pred_err = h.Clone(h.GetName()+"_err")
    vhpred = []
    vfpred = []
    if h.GetMaximum()>20:
        #for ndays in range(max(ndaysback-ndayspread*2,2),ndaysback+ndayspread*2+1):
        for ndays in range(-ndayspread,ndayspread+1):
            hpred,fpred = get_pred(dh,h,ndaysback+ndays*2)
            vhpred.append(hpred)
            vfpred.append(fpred)
        for binx in range(1,hpred.GetNbinsX()+1):
            if vhpred[0].GetBinContent(binx)==0:
                h_pred_nom.SetBinContent(binx,0)
                h_pred_nom.SetBinError  (binx,0)
                h_pred_err.SetBinContent(binx,0)
                h_pred_err.SetBinError  (binx,0)
            else:
                # calculate mean
                logs = []
                vals = []
                for ifit in range(len(vhpred)):
                    logs.append(math.log(vhpred[ifit].GetBinContent(binx)))
                    vals.append(vhpred[ifit].GetBinContent(binx))
                #print binx
                #print logs
                #print vals
                #mean = math.exp(np.mean(logs))
                mean = np.mean(vals)
                #print "mean1="+str(mean)
                #print "mean2="+str(mean2)
                #print "stdev="+str(np.std(logs))
                #print "stdv2="+str(np.std(vals))
                h_pred_nom.SetBinContent(binx, mean)
                #h_pred_nom.SetBinError  (binx, np.std(vals))
                h_pred_nom.SetBinError  (binx,0)
                # calculate error up/down separately
                vup = []
                vdn = []
                for ifit in range(len(vals)):
                    if vals[ifit]>mean:
                        vup.append(vals[ifit])
                    else:
                        vdn.append(vals[ifit])
                err_up = ((np.mean(vup) if len(vup) else 0) - mean)/mean
                err_up = ((err_up**2 + diff**2)**0.5)*mean
                err_dn = (mean - (np.mean(vdn) if len(vdn) else 0))/mean
                err_dn = ((err_dn**2 + diff**2)**0.5)*mean
                #print binx
                #print mean
                #print (np.mean(vup)+np.mean(vdn))/2
                #print (np.mean(vup)-np.mean(vdn))/2
                #h_pred_err.SetBinContent(binx, (np.mean(vup)+np.mean(vdn))/2)
                #h_pred_err.SetBinError  (binx, (np.mean(vup)-np.mean(vdn))/2)
                #print err_up
                #print err_dn
                #print mean + (err_up-err_dn)/2
                #print (err_up-err_dn)/2
                h_pred_err.SetBinContent(binx, mean + (err_up-err_dn)/2)
                h_pred_err.SetBinError  (binx, (err_up+err_dn)/2)
    return [h_pred_nom, h_pred_err, vfpred]

# Coronavirus data
data = {}
data["WRL"]  = [ ROOT.TH1D("WRL_dchg", ";Date;2nd derivative of growth rate (%)", 366,0,366*86400),
                 ROOT.TH1D("WRL_chg",  ";Date;1st derivative of growth rate (%)", 366,0,366*86400),
                 ROOT.TH1D("WRL_gro",  ";Date;Daily incrase (%)",                 366,0,366*86400),
                 ROOT.TH1D("WRL",      "World;Date;Total comfirmed cases",        366,0,366*86400)]
death = {}
death["WRL"] = [ ROOT.TH1D("WRL_death_dchg", ";Date;2nd der. of death growth rate (%)", 366,0,366*86400),
                 ROOT.TH1D("WRL_death_chg",  ";Date;1st der. of death growth rate (%)", 366,0,366*86400),
                 ROOT.TH1D("WRL_death_gro",  ";Date;Daily death incrase (%)",           366,0,366*86400),
                 ROOT.TH1D("WRL_death",      "World;Date;Total comfirmed death",        366,0,366*86400) ]
monthtoday = { 1:0, 2:31, 3:60, 4:91, 5:121, 6:152, 7:182, 8:213, 9:244, 10:274, 11:305, 12:335 }
population = {"WRL":0}
total = {}
dtotal = {}
lastday = {}
lastcount = {}
lastdcount = {}
lastcountry = ""
allcountries = []
valid_day = 0
valid_count = 0
countrynames = {}
# ---------------- ECDC ----------------
# Update case data if older than an hour
update = 0
if not os.path.exists('corona_data.csv'): update = 1
elif (time.time() - os.path.getmtime('corona_data.csv')) > 3600: update = 1
if update:
    print "Updating data from ECDC servers"
    filedata = urllib2.urlopen('https://opendata.ecdc.europa.eu/covid19/casedistribution/csv')
    datatowrite = filedata.read()
    with open('corona_data.csv', 'wb') as f:
        f.write(datatowrite)
print "Start reading world covid-19 data (ECDC)"
for line in reversed(open("corona_data.csv").readlines()[1:]):
    info = line.rstrip().split(",")
    country = info[8]
    population[country] = int(info[9]) if info[9] != "" else 1
    # correcting
    # missing info
    if country == "":
        geoId = info[7]
        if geoId == "FK":
            country = "FLK"
            population[country] = 2840
        elif geoId == "CZ":
            country = "CZE"
            population[country] = 10650000 # 2019
        elif geoId == "AI":
            country = "AIA"
            population[country] = 15094
        else:
            country = geoId
            population[country] = 1
    # wrong daily stat for Hungary
    if country == "HUN":
        if info[0] == "31/03/2020": info[4] = 45
        if info[0] == "01/04/2020": info[4] = 33
        if info[0] == "02/04/2020": info[4] = 60
        if info[0] == "03/04/2020": info[4] = 38
        if info[0] == "04/04/2020": info[4] = 55
    if country != lastcountry:
        if country not in allcountries:
            allcountries.append(country)
            population["WRL"] += (int(info[9]) if info[9] != "" else 0)
            dddh = ROOT.TH1D(country+"_dchg",                 ";Date;2nd drivative of growth rate (%)", 366,0,366*86400)
            ddh  = ROOT.TH1D(country+"_chg",                  ";Date;1st drivative of growth rate (%)", 366,0,366*86400)
            dh   = ROOT.TH1D(country+"_gro",                  ";Date;Daily growth (%)",                 366,0,366*86400)
            h    = ROOT.TH1D(country,info[6].replace("_"," ")+";Date;Total comfirmed cases",            366,0,366*86400)
            data[country] = [dddh, ddh, dh, h]
            dddd = ROOT.TH1D(country+"_death_dchg",                    ";Date;2nd der. of death growth rate (%)", 366,0,366*86400)
            ddd  = ROOT.TH1D(country+"_death_chg",                     ";Date;1st der. of death growth rate (%)", 366,0,366*86400)
            dd   = ROOT.TH1D(country+"_death_gro",                     ";Date;Daily death growth (%)",            366,0,366*86400)
            d    = ROOT.TH1D(country+"_death",info[6].replace("_"," ")+";Date;Total comfirmed death",             366,0,366*86400)
            death[country] = [dddd, ddd, dd, d]
            countrynames[country] = info[6].replace("_"," ")
        if country not in total:
            total[country] = 0
            dtotal[country] = 0
            lastday[country] = 0
            lastcount[country] = 0
            lastdcount[country] = 0
    total[country] += int(info[4])
    dtotal[country] += int(info[5])
    count = total[country]
    dcount = dtotal[country]
    day = monthtoday[int(info[2])]+int(info[1])
    for currday in range(lastday[country]+1,day):
        h.SetBinContent(currday, lastcount[country])
        d.SetBinContent(currday, lastdcount[country])
        #h.SetBinError  (currday, 0)
    h.SetBinContent(day, count)
    d.SetBinContent(day, dcount)
    #h.SetBinError  (day, 0)
    h2 = data["WRL"][3]
    h2.Fill((day-1)*86400, count)
    h2.SetBinError(day,0)
    d2 = death["WRL"][3]
    d2.Fill((day-1)*86400, dcount)
    d2.SetBinError(day,0)
    lastcountry = country
    lastday[country] = day
    lastcount[country] = count
    lastdcount[country] = dcount
    if country == "HUN":
        valid_day = day
        valid_count = count
# --------------- Testing --------------
# Update testing data if older than 4 hours
update = 0
if not os.path.exists('covid-testing-all-observations.csv'): update = 1
elif (time.time() - os.path.getmtime('covid-testing-all-observations.csv')) > 3600: update = 1
if update:
    print "Updating testing data from github"
    filedata = urllib2.urlopen('https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/testing/covid-testing-all-observations.csv')
    datatowrite = filedata.read()
    with open('covid-testing-all-observations.csv', 'wb') as f:
        f.write(datatowrite)
        # I need to manually add the data for Hungary
        f.write("Hungary - tests performed,2020-03-04,http://abouthungary.hu/,Ministry of Interior,,230,,,\n")
        f.write("Hungary - tests performed,2020-03-05,http://abouthungary.hu/,Ministry of Interior,,230,,,\n")
        f.write("Hungary - tests performed,2020-03-06,http://abouthungary.hu/,Ministry of Interior,,269,,,\n")
        f.write("Hungary - tests performed,2020-03-07,http://abouthungary.hu/,Ministry of Interior,,319,,,\n")
        f.write("Hungary - tests performed,2020-03-08,http://abouthungary.hu/,Ministry of Interior,,362,,,\n")
        f.write("Hungary - tests performed,2020-03-09,http://abouthungary.hu/,Ministry of Interior,,362,,,\n")
        f.write("Hungary - tests performed,2020-03-10,http://abouthungary.hu/,Ministry of Interior,,531,,,\n")
        f.write("Hungary - tests performed,2020-03-11,http://abouthungary.hu/,Ministry of Interior,,609,,,\n")
        f.write("Hungary - tests performed,2020-03-12,http://abouthungary.hu/,Ministry of Interior,,730,,,\n")
        f.write("Hungary - tests performed,2020-03-13,http://abouthungary.hu/,Ministry of Interior,,858,,,\n")
        f.write("Hungary - tests performed,2020-03-14,http://abouthungary.hu/,Ministry of Interior,,1014,,,\n")
        f.write("Hungary - tests performed,2020-03-15,http://abouthungary.hu/,Ministry of Interior,,1236,,,\n")
        f.write("Hungary - tests performed,2020-03-16,http://abouthungary.hu/,Ministry of Interior,,1470,,,\n")
        f.write("Hungary - tests performed,2020-03-17,http://abouthungary.hu/,Ministry of Interior,,1587,,,\n")
        f.write("Hungary - tests performed,2020-03-18,http://abouthungary.hu/,Ministry of Interior,,1803,,,\n")
        f.write("Hungary - tests performed,2020-03-19,http://abouthungary.hu/,Ministry of Interior,,2322,,,\n")
        f.write("Hungary - tests performed,2020-03-20,http://abouthungary.hu/,Ministry of Interior,,3007,,,\n")
        f.write("Hungary - tests performed,2020-03-21,http://abouthungary.hu/,Ministry of Interior,,3477,,,\n")
        f.write("Hungary - tests performed,2020-03-22,http://abouthungary.hu/,Ministry of Interior,,4443,,,\n")
        f.write("Hungary - tests performed,2020-03-23,http://abouthungary.hu/,Ministry of Interior,,5515,,,\n")
        f.write("Hungary - tests performed,2020-03-24,http://abouthungary.hu/,Ministry of Interior,,6113,,,\n")
        f.write("Hungary - tests performed,2020-03-25,http://abouthungary.hu/,Ministry of Interior,,6817,,,\n")
        f.write("Hungary - tests performed,2020-03-26,http://abouthungary.hu/,Ministry of Interior,,8005,,,\n")
        f.write("Hungary - tests performed,2020-03-27,http://abouthungary.hu/,Ministry of Interior,,9275,,,\n")
        f.write("Hungary - tests performed,2020-03-28,http://abouthungary.hu/,Ministry of Interior,,10303,,,\n")
        f.write("Hungary - tests performed,2020-03-29,http://abouthungary.hu/,Ministry of Interior,,12148,,,\n")
        f.write("Hungary - tests performed,2020-03-30,http://abouthungary.hu/,Ministry of Interior,,13301,,,\n")
        f.write("Hungary - tests performed,2020-03-31,http://abouthungary.hu/,Ministry of Interior,,14146,,,\n")
        f.write("Hungary - tests performed,2020-04-01,http://abouthungary.hu/,Ministry of Interior,,15208,,,\n")
        f.write("Hungary - tests performed,2020-04-02,http://abouthungary.hu/,Ministry of Interior,,16401,,,\n")
        f.write("Hungary - tests performed,2020-04-03,http://abouthungary.hu/,Ministry of Interior,,17769,,,\n")
        f.write("Hungary - tests performed,2020-04-04,http://abouthungary.hu/,Ministry of Interior,,19424,,,\n")
        f.write("Hungary - tests performed,2020-04-05,http://abouthungary.hu/,Ministry of Interior,,21250,,,\n")
        f.write("Hungary - tests performed,2020-04-06,http://abouthungary.hu/,Ministry of Interior,,22282,,,\n")
        f.write("Hungary - tests performed,2020-04-07,http://abouthungary.hu/,Ministry of Interior,,23746,,,\n")
        f.write("Hungary - tests performed,2020-04-08,http://abouthungary.hu/,Ministry of Interior,,25748,,,\n")
        f.write("Hungary - tests performed,2020-04-09,http://abouthungary.hu/,Ministry of Interior,,27826,,,\n")
        f.write("Hungary - tests performed,2020-04-10,http://abouthungary.hu/,Ministry of Interior,,29948,,,\n")
        f.write("Hungary - tests performed,2020-04-11,http://abouthungary.hu/,Ministry of Interior,,31961,,,\n")
        f.write("Hungary - tests performed,2020-04-12,http://abouthungary.hu/,Ministry of Interior,,33532,,,\n")
        f.write("Hungary - tests performed,2020-04-13,http://abouthungary.hu/,Ministry of Interior,,34819,,,\n")
        f.write("Hungary - tests performed,2020-04-14,http://abouthungary.hu/,Ministry of Interior,,35825,,,\n")
        f.write("Hungary - tests performed,2020-04-15,http://abouthungary.hu/,Ministry of Interior,,37326,,,\n")
        f.write("Hungary - tests performed,2020-04-16,http://abouthungary.hu/,Ministry of Interior,,38489,,,\n")
        f.write("Hungary - tests performed,2020-04-16,http://abouthungary.hu/,Ministry of Interior,,41590,,,\n")
print "Start reading world testing data (github)"
testcountries = []
testdata = {}
testdata["WRL"]  = [ ROOT.TH1D("WRL_test_gro",  ";Date;Daily growth in tests (%)", 366,0,366*86400),
                     ROOT.TH1D("WRL_test", "World;Date;Total tests done",          366,0,366*86400),
                     "sample"]
lasttestday = {}
lasttestcount = {}
lasttestcountry = ""
testtotal = {}
nline = 0
with open('covid-testing-all-observations.csv') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    for info in csvreader:
        if nline != 0 and info[0] != '' and info[1] != '' and info[5] != '':
            countryname = info[0].split(" - ")[0]
            if countryname == "Czech Republic": countryname = "Czechia"
            testinginfo = info[0].split(" - ")[1]
            datatype = "sample"
            norm = 1
            if "people" in testinginfo or "cases" in testinginfo:
                datatype = "case"
                norm = 1.5
            country = ""
            for short, name in countrynames.items():
                if name == countryname:
                    country = short
                elif countryname in name:
                    country = short
            if countryname == "Hong Kong": country = "HKG"
            # corrections
            if country != lasttestcountry:
                #if lasttestcountry != '': print countrynames[lasttestcountry]+" "+str(lasttestcount[lasttestcountry])
                if country not in testcountries:
                    #print "TESTINFO: "+country+" "+testinginfo
                    testcountries.append(country)
                    dtest   = ROOT.TH1D(country+"_test_gro",        ";Date;Daily growth in tests (%)", 366,0,366*86400)
                    test    = ROOT.TH1D(country+"_test",countryname+";Date;Total tests done",          366,0,366*86400)
                    testdata[country] = [dtest, test, datatype]
                if country not in testtotal:
                    testtotal[country] = 0
                    lasttestday[country] = 0
                    lasttestcount[country] = 0
            testtotal[country] = int(info[5]) * norm
            testcount = testtotal[country]
            day = monthtoday[int(info[1].split("-")[1])]+int(info[1].split("-")[2])
            for currday in range(lasttestday[country]+1,day):
                test.SetBinContent(currday, lasttestcount[country])
            test.SetBinContent(day, testcount)
            t2 = data["WRL"][1]
            t2.Fill((day-1)*86400, testcount)
            t2.SetBinError(day,0)
            lasttestcountry = country
            lasttestday[country] = day
            lasttestcount[country] = testcount
        nline += 1
#if lasttestcountry != '': print countrynames[lasttestcountry]+" "+str(lasttestcount[lasttestcountry])
print "Reading world data done"
print

# More methods
def get_date(day):
    maxdiff = 9999
    month = ""
    months = { 1:"Jan", 2:"Feb", 3:"Mar", 4:"Apr", 5:"May", 6:"Jun", 7:"Jul", 8:"Aug", 9:"Sep", 10:"Oct", 11:"Nov", 12:"Dec" }
    for key in monthtoday.keys():
        diff = day - monthtoday[key]
        if diff>0 and diff<maxdiff:
            maxdiff = diff
            month = months[key]
    return month+" "+("%02d" % maxdiff)

def get_pred_from_death(h, d, htest):
    nday_back = 12
    mortality = 0.01
    #mortality = 0.007
    ht   = h.Clone(h.GetName()+"_true")
    dt   = ht.Clone(ht.GetName()+"_gro")
    hdoc = h.Clone(h.GetName()+"_doc")
    hc   = h.Clone(h.GetName()+"_cfr")
    hef  = h.Clone(h.GetName()+"_exc")
    hed  = h.Clone(h.GetName()+"_exc")
    min_cfr = 9999
    # Moving cfr
    for binx in range(1,d.GetNbinsX()+1):
        if (h.GetBinContent(binx-1)>0 and h.GetBinContent(binx)>0 and h.GetBinContent(binx+1)>0 and
            d.GetBinContent(binx+12-1)>0 and d.GetBinContent(binx+12)>0 and d.GetBinContent(binx+12+1)>0):
            case0   = h.GetBinContent(binx-1) + h.GetBinContent(binx) + h.GetBinContent(binx+1)
            death12 = d.GetBinContent(binx+12-1) + d.GetBinContent(binx+12) + d.GetBinContent(binx+12+1)
            cfr = death12/case0
            if cfr>0 and cfr<min_cfr:
                min_cfr = cfr
            hc.SetBinContent(binx, cfr)
        else:
            hc.SetBinContent(binx,0)
        hc.SetBinError  (binx,0)
    # Excess fatality and death
    for binx in range(1,h.GetNbinsX()+1):
        death12 = d.GetBinContent(binx+12)-d.GetBinContent(binx+11)
        cfr     = hc.GetBinContent(binx)
        if cfr>0:
            hef.SetBinContent(binx, (cfr/min_cfr/2-1)*100)
            hed.SetBinContent(binx+12, (1-min_cfr/cfr/2)*death12)
        else:
            hef.SetBinContent(binx, 0)
            hed.SetBinContent(binx+12, 0)
        hef.SetBinError(binx, 0)
        hed.SetBinError(binx, 0)
    # Estimated true cases
    for binx in range(1,d.GetNbinsX()+1):
        death = d.GetBinContent(binx)
        if death>0:
            ht.SetBinContent(binx-nday_back, death/mortality)
        else:
            ht.SetBinContent(binx-nday_back, 0)
    n = 0
    for binx in range(1,h.GetNbinsX()+1):
        case  = h.GetBinContent(binx)
        true  = ht.GetBinContent(binx)
        test  = 0 if "nodata" in htest.GetName() else htest.GetBinContent(binx)
        death = d.GetBinContent(binx)
        if case>0 and true>0:
            hdoc.SetBinContent(binx,case/true*100)
            if test>0 and death>0:
                n += 1
    gtest = ROOT.TGraph(n)
    n = 0
    for binx in range(1,h.GetNbinsX()+1):
        case  = h.GetBinContent(binx)
        true  = ht.GetBinContent(binx)
        test  = 0 if "nodata" in htest.GetName() else htest.GetBinContent(binx)
        death = d.GetBinContent(binx)
        if case>0 and true>0 and test>0 and death>0:
            gtest.SetPoint(n, test/death, case/true*100)
        n += 1
    # perform linear fit to observed true ratio
    #rt = ht.Clone(ht.GetName()+"_rt")
    t_pred_nom = h.Clone(h.GetName()+"_pred_nom")
    t_pred_err = h.Clone(h.GetName()+"_err")
    vtpred = []  # histo for each vairable range fit
    vftpred = [] # fit for each range
    if ht.GetMaximum()>20:
        daily_growth(ht, dt)
        for ndays in range(-ndayspread,ndayspread+1):
            tpred,fpred = get_pred(dt,ht,ndaysback2+ndays*1)
            vtpred.append(tpred)
            vftpred.append(fpred)
        for binx in range(1,tpred.GetNbinsX()+1):
            if vtpred[0].GetBinContent(binx)==0:
                t_pred_nom.SetBinContent(binx,0)
                t_pred_nom.SetBinError  (binx,0)
                t_pred_err.SetBinContent(binx,0)
                t_pred_err.SetBinError  (binx,0)
            else:
                # error calculation
                # calculate mean
                logs = []
                vals = []
                for ifit in range(len(vtpred)):
                    logs.append(math.log(vtpred[ifit].GetBinContent(binx)))
                    vals.append(vtpred[ifit].GetBinContent(binx))
                #print binx
                #print logs
                #print vals
                #mean = math.exp(np.mean(logs))
                mean = np.mean(vals)
                #print "mean1="+str(mean)
                #print "mean2="+str(mean2)
                #print "stdev="+str(np.std(logs))
                #print "stdv2="+str(np.std(vals))
                t_pred_nom.SetBinContent(binx, mean)
                #t_pred_nom.SetBinError  (binx, np.std(vals))
                t_pred_nom.SetBinError  (binx,0)
                # calculate error up/down separately
                vup = []
                vdn = []
                for ifit in range(len(vals)):
                    if vals[ifit]>mean:
                        vup.append(vals[ifit])
                    else:
                        vdn.append(vals[ifit])
                #print binx
                #print mean
                #print (np.mean(vup)+np.mean(vdn))/2
                #print (np.mean(vup)-np.mean(vdn))/2
                #t_pred_err.SetBinContent(binx, (np.mean(vup)+np.mean(vdn))/2)
                #t_pred_err.SetBinError  (binx, (np.mean(vup)-np.mean(vdn))/2)
                err_up = ((np.mean(vup) if len(vup) else 0) - mean)/mean
                err_up = ((err_up**2 + diff**2)**0.5)*mean
                err_dn = (mean - (np.mean(vdn) if len(vdn) else 0))/mean
                err_dn = ((err_dn**2 + diff**2)**0.5)*mean
                #print err_up
                #print err_dn
                #print mean + (err_up-err_dn)/2
                #print (err_up-err_dn)/2
                t_pred_err.SetBinContent(binx, mean + (err_up-err_dn)/2)
                t_pred_err.SetBinError  (binx, (err_up+err_dn)/2)
    return [ht, dt, hdoc, hc, hef, hed, gtest, t_pred_nom, t_pred_err, vftpred]


# Do all the calculations and predictions and save results to vectors
fractions = {}
preddata = {}
for country in allcountries+["WRL"]:
    # comfirmed cases
    dddh       = data[country][0]
    ddh        = data[country][1]
    dh         = data[country][2]
    h          = data[country][3]
    diff = daily_growth_and_diff(h, dh, ddh, dddh)
    h_pred_nom, h_pred_err, vf = get_pred_and_error(h, dh, diff)
    data[country] += [h_pred_nom, h_pred_err, vf]
    # deaths
    dddd       = death[country][0]
    ddd        = death[country][1]
    dd         = death[country][2]
    d          = death[country][3]
    diff2 = daily_growth_and_diff(d, dd, ddd, dddd)
    d_pred_nom, d_pred_err, vdf = get_pred_and_error(d, dd, diff2)
    death[country] += [d_pred_nom, d_pred_err, vdf]
    # tests data
    if country in testdata:
        test = testdata[country][1]
    else:
        test = ROOT.TH1D(country+"_nodata","",1,0,1)
    # predictions from death
    preddata[country] = [test] + get_pred_from_death(h, d, test)
    # growth rates
    curr_case = h.GetBinContent(valid_day)
    curr_rate = dh.GetBinContent(valid_day)
    t_pred_nom = preddata[country][8]
    fraction = t_pred_nom.GetBinContent(valid_day+days_pred)/population[country] if population[country]>1e5 else 0
    if population[country]>min_pop and curr_case>case_threshold and country != "WRL":
        fractions[country] = fraction

# Auto select 8 worst countries based on worst projections to 14 days
if worst:
    print "Plot 8 worst pandemic threats in 14 days"
    countries = []
print "Countries with highest expected true infected fraction in "+str(days_pred)+" days:"
for ctry in reversed(sorted(fractions.items(), key=operator.itemgetter(1))[-14:]):
    #dddh, ddh, dh, h, h_pred_nom, h_pred_err, vf = data[ctry[0]]
    t_pred_nom = preddata[ctry[0]][8]
    t_pred_err = preddata[ctry[0]][9]
    fraction = t_pred_nom.GetBinContent(valid_day+days_pred)/population[ctry[0]]
    print ("%-30s - Fraction infected: %.3f%%" % (data[ctry[0]][3].GetTitle()+" ("+ctry[0]+")", fraction*100))
    for day in [7,14,21,28]:
        fraction = t_pred_nom.GetBinContent(valid_day+day)/population[ctry[0]]
        print ("%s: %d  +%d  -%d, fraction: %.3f%%" % (get_date(valid_day+day), t_pred_nom.GetBinContent(valid_day+day),
                                                     t_pred_err.GetBinContent(valid_day+day)+t_pred_err.GetBinError(valid_day+day)-t_pred_nom.GetBinContent(valid_day+day),
                                                     t_pred_nom.GetBinContent(valid_day+day)-t_pred_err.GetBinContent(valid_day+day)+t_pred_err.GetBinError(valid_day+day),
                                                     fraction*100))
    if worst and len(countries)<8: countries.append(ctry[0])
print

daily_rates = {}
for country in allcountries:
    dddh       = data[country][0]
    ddh        = data[country][1]
    dh         = data[country][2]
    h          = data[country][3]
    curr_case = h.GetBinContent(valid_day)
    curr_rate = dh.GetBinContent(valid_day)
    if population[country]>best_min_pop and curr_case>best_case_threshold:
        daily_rates[country] = curr_rate

# Auto select 8 best countries based on lowest current growth rate
if best:
    print "Plot 8 best countries currently"
    countries = []
print "Countries with lowest current growth rate:"
for ctry in sorted(daily_rates.items(), key=operator.itemgetter(1))[:14]:
    dddh, ddh, dh, h, h_pred_nom, h_pred_err, vf = data[ctry[0]]
    daily_rate = ctry[1]
    fraction = h.GetBinContent(valid_day)/population[ctry[0]]
    print ("%-30s - Current rate: %.4f  Fraction infected: %.3f%%" % (h.GetTitle()+" ("+ctry[0]+")", daily_rate, fraction*100))
    if best and len(countries)<8: countries.append(ctry[0])
print

# Save plots for selected countries to vectors
vdddh = []
vddh  = []
vdh   = []
vh    = []
vh_pred = []
vh_err  = []
vvf = []
vdddd = []
vddd  = []
vdd   = []
vd    = []
vd_pred = []
vd_err  = []
vvdf = []
vht = []
vdt = []
vhdoc = []
vhc = []
vhef = []
vhed = []
vgtest = []
vt_pred_nom = []
vt_pred_err = []
vvftpred = []
vtest = []
for i in range(len(countries)):
    dddh, ddh, dh, h, h_pred_nom, h_pred_err, vf = data[countries[i]]
    vdddh.append(dddh)
    vddh.append(ddh)
    vdh.append(dh)
    vh.append(h)
    vh_pred.append(h_pred_nom)
    vh_err .append(h_pred_err)
    vvf.append(vf)
    dddd, ddd, dd, d, d_pred_nom, d_pred_err, vdf = death[countries[i]]
    vdddd.append(dddd)
    vddd.append(ddd)
    vdd.append(dd)
    vd.append(d)
    vd_pred.append(d_pred_nom)
    vd_err .append(d_pred_err)
    vvdf.append(vdf)
    test, ht, dt, hdoc, hc, hef, hed, gtest, t_pred_nom, t_pred_err, vftpred = preddata[countries[i]]
    vtest.append(test)
    vht.append(ht)
    vdt.append(dt)
    vhdoc.append(hdoc)
    vhc.append(hc)
    vhef.append(hef)
    vhed.append(hed)
    vgtest.append(gtest)
    vt_pred_nom.append(t_pred_nom)
    vt_pred_err.append(t_pred_err)
    vvftpred.append(vftpred)

# If data is unreliable, you can validate latest comfirmed case numbers with this printout
print "Count validations for date="+get_date(valid_day)+" (day="+str(valid_day)+"):"
for i in range(len(countries)):
    print ("- %-30s: %7d (death: %6d)" %(data[countries[i]][3].GetTitle(), data[countries[i]][3].GetBinContent(valid_day),
                                         death[countries[i]][3].GetBinContent(valid_day)))
print

# Calculate some trends based on exponential slope changes
# Does not really work, seems to show the obvious
if trends:
    print "Outbrakes and trend changes:"
    for i in range(len(countries)):
        print "- "+data[countries[i]][3].GetTitle()
        dddh = data[countries[i]][0]
        outbreak = True
        slowdown = False
        for binx in range(1,dddh.GetNbinsX()+1):
            prevcont = dddh.GetBinContent(binx-1)
            cont     = dddh.GetBinContent(binx)
            nextcont = dddh.GetBinContent(binx+1)
            if cont!=0:
                if not slowdown and cont>0 and not prevcont>0:
                    slowdown = False
                    outbreak = True
                if not outbreak and cont<0 and not prevcont<0:
                    slowdown = True
                    outbreak = False
                if outbreak and cont>1:
                    if nextcont<cont:
                        print ("  + %s: Outbreak - +%.1f" % (get_date(binx), cont))
                        outbreak = False
                if slowdown and cont<-1:
                    if nextcont>cont:
                        print ("  + %s: Slowdown - %.1f" % (get_date(binx), cont))
                        slowdown = False

# Plotting of results starts here
print "Start plotting"
savedir = ""
if save:
    savedir = "results/"+(get_date(valid_day).replace(" ","_"))
    if not os.path.exists(savedir): os.makedirs(savedir)
    savedir += "/"

set_default_style_()
if trends:
    can = ROOT.TCanvas(countries[0]+"_chg", "", 1200,600)
    can.Divide(2)
    pad = can.cd(1)
    pad.SetTopMargin  (0.04)
    pad.SetLeftMargin (0.13)
    pad.SetRightMargin(0.04)
    pad.SetGrid(1,1)
    vdddh[0].GetXaxis().SetTitleSize(0.05)
    vdddh[0].GetYaxis().SetTitleSize(0.05)
    vdddh[0].GetXaxis().SetTimeDisplay(1)
    vdddh[0].GetXaxis().SetNdivisions(504)
    vdddh[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
    vdddh[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
    vdddh[0].GetYaxis().SetRangeUser(-20,20)
    leg1a = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
    leg1a.SetFillColor(0)
    leg1a.SetFillStyle(0)
    leg1a.SetBorderSize(0)
    leg1a.SetTextSize(0.04)
    for i in range(len(countries)):
        vdddh[i].SetLineColor(settings[i][2])
        vdddh[i].SetLineWidth(2)
        vdddh[i].Draw("HIST" if i==0 else "SAME HIST")
        leg1a.AddEntry(vdddh[i], vh[i].GetTitle(), "L")
    leg1a.Draw()
    
    pad = can.cd(2)
    pad.SetTopMargin  (0.04)
    pad.SetLeftMargin (0.13)
    pad.SetRightMargin(0.04)
    pad.SetGrid(1,1)
    vddh[0].GetXaxis().SetTitleSize(0.05)
    vddh[0].GetYaxis().SetTitleSize(0.05)
    vddh[0].GetXaxis().SetTimeDisplay(1)
    vddh[0].GetXaxis().SetNdivisions(504)
    vddh[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
    vddh[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
    vddh[0].GetYaxis().SetRangeUser(-20,20)
    leg1b = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
    leg1b.SetFillColor(0)
    leg1b.SetFillStyle(0)
    leg1b.SetBorderSize(0)
    leg1b.SetTextSize(0.04)
    for i in range(len(countries)):
        vddh[i].SetLineColor(settings[i][2])
        vddh[i].SetLineWidth(2)
        vddh[i].Draw("HIST" if i==0 else "SAME HIST")
        leg1b.AddEntry(vddh[i], vh[i].GetTitle(), "L")
    leg1b.Draw()
    can.SaveAs(savedir+countries[0]+"_diff.png")

can = ROOT.TCanvas(countries[0]+"comfirmed", "", 1200,600)
can.Divide(2)
pad = can.cd(1)
pad.SetTopMargin  (0.04)
pad.SetLeftMargin (0.13)
pad.SetRightMargin(0.04)
pad.SetGrid(1,1)
pad.SetLogy(logy2)
vdh[0].GetXaxis().SetTitleSize(0.05)
vdh[0].GetYaxis().SetTitleSize(0.05)
vdh[0].GetXaxis().SetTimeDisplay(1)
vdh[0].GetXaxis().SetNdivisions(504)
vdh[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
vdh[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
vdh[0].GetYaxis().SetRangeUser(yrange2[0],yrange2[1])
leg2a = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
leg2a.SetFillColor(0)
leg2a.SetFillStyle(0)
leg2a.SetBorderSize(0)
leg2a.SetTextSize(0.04)
for i in range(len(countries)):
    vdh[i].SetMarkerStyle(settings[i][0])
    vdh[i].SetMarkerColor(settings[i][2])
    vdh[i].SetMarkerSize (settings[i][3]*0.8)
    vdh[i].SetLineColor(settings[i][2])
    vdh[i].Draw("PL" if i==0 else "SAME PL")
    leg2a.AddEntry(vdh[i], vh[i].GetTitle(), "P")
    for fit in vvf[i]:
        fit.SetLineColor(settings[i][2])
        fit.SetLineStyle(7)
        fit.SetLineWidth(1)
        fit.Draw("SAME")
leg2a.Draw()

vh[0].GetXaxis().SetTimeDisplay(1)
vh[0].GetXaxis().SetNdivisions(504)
vh[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
vh[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
vh[0].GetYaxis().SetRangeUser(yrange[0],yrange[1])
leg2b = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
leg2b.SetFillColor(0)
leg2b.SetFillStyle(0)
leg2b.SetBorderSize(0)
leg2b.SetTextSize(0.04)
for i in range(len(countries)):
    leg2b.AddEntry(vh[i], vh[i].GetTitle(), "P")
    if i==0:
        pad = can.cd(2)
        pad.SetTopMargin  (0.04)
        pad.SetLeftMargin (0.13)
        pad.SetRightMargin(0.04)
        pad.SetLogy(1)
        pad.SetGrid(1,1)
    vh[i].GetXaxis().SetTitleSize(0.05)
    vh[i].GetYaxis().SetTitleSize(0.05)
    vh[i].SetMarkerStyle(settings[i][0])
    vh[i].SetMarkerColor(settings[i][2])
    vh[i].SetMarkerSize (settings[i][3])
    vh[i].SetLineColor  (settings[i][2])
    vh[i].Draw("P" if i==0 else "SAME P")
    vh_err[i].SetMarkerStyle(1)
    vh_err[i].SetMarkerColor(0)
    vh_err[i].SetMarkerSize(0)
    vh_err[i].SetFillStyle(settings[i][4])
    vh_err[i].SetFillColor(settings[i][2])
    vh_err[i].Draw("SAME E4")
    vh_pred[i].SetMarkerStyle(settings[i][1])
    vh_pred[i].SetMarkerColor(settings[i][2])
    vh_pred[i].SetMarkerSize (settings[i][3])
    vh_pred[i].SetLineColor  (settings[i][2])
    vh_pred[i].Draw("SAME P")
leg2b.Draw()
can.SaveAs(savedir+countries[0]+"_comfirmed.png")

can = ROOT.TCanvas(countries[0]+"_dead", "", 1200,600)
can.Divide(2)
pad = can.cd(1)
pad.SetTopMargin  (0.04)
pad.SetLeftMargin (0.13)
pad.SetRightMargin(0.04)
pad.SetGrid(1,1)
pad.SetLogy(logy2)
vdd[0].GetXaxis().SetTitleSize(0.05)
vdd[0].GetYaxis().SetTitleSize(0.05)
vdd[0].GetXaxis().SetTimeDisplay(1)
vdd[0].GetXaxis().SetNdivisions(504)
vdd[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
vdd[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
vdd[0].GetYaxis().SetRangeUser(yrange2[0],yrange2[1])
leg3a = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
leg3a.SetFillColor(0)
leg3a.SetFillStyle(0)
leg3a.SetBorderSize(0)
leg3a.SetTextSize(0.04)
for i in range(len(countries)):
    vdd[i].SetMarkerStyle(settings[i][0])
    vdd[i].SetMarkerColor(settings[i][2])
    vdd[i].SetMarkerSize (settings[i][3]*0.8)
    vdd[i].SetLineColor(settings[i][2])
    vdd[i].Draw("PL" if i==0 else "SAME PL")
    leg3a.AddEntry(vdd[i], vh[i].GetTitle(), "P")
    for fit in vvdf[i]:
        fit.SetLineColor(settings[i][2])
        fit.SetLineStyle(7)
        fit.SetLineWidth(1)
        fit.Draw("SAME")
leg3a.Draw()

vd[0].GetXaxis().SetTimeDisplay(1)
vd[0].GetXaxis().SetNdivisions(504)
vd[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
vd[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
vd[0].GetYaxis().SetRangeUser(yrange[0],yrange[1])
leg3b = ROOT.TLegend(0.15,y2leg-len(vd)*0.04,0.35,y2leg)
leg3b.SetFillColor(0)
leg3b.SetFillStyle(0)
leg3b.SetBorderSize(0)
leg3b.SetTextSize(0.04)
for i in range(len(countries)):
    leg3b.AddEntry(vd[i], vh[i].GetTitle(), "P")
    if i==0:
        pad = can.cd(2)
        pad.SetTopMargin  (0.04)
        pad.SetLeftMargin (0.13)
        pad.SetRightMargin(0.04)
        pad.SetLogy(1)
        pad.SetGrid(1,1)
    vd[i].GetXaxis().SetTitleSize(0.05)
    vd[i].GetYaxis().SetTitleSize(0.05)
    vd[i].SetMarkerStyle(settings[i][0])
    vd[i].SetMarkerColor(settings[i][2])
    vd[i].SetMarkerSize (settings[i][3])
    vd[i].SetLineColor  (settings[i][2])
    vd[i].Draw("P" if i==0 else "SAME P")
    vd_err[i].SetMarkerStyle(1)
    vd_err[i].SetMarkerColor(0)
    vd_err[i].SetMarkerSize(0)
    vd_err[i].SetFillStyle(settings[i][4])
    vd_err[i].SetFillColor(settings[i][2])
    vd_err[i].Draw("SAME E4")
    vd_pred[i].SetMarkerStyle(settings[i][1])
    vd_pred[i].SetMarkerColor(settings[i][2])
    vd_pred[i].SetMarkerSize (settings[i][3])
    vd_pred[i].SetLineColor  (settings[i][2])
    vd_pred[i].Draw("SAME P")
leg3b.Draw()
can.SaveAs(savedir+countries[0]+"_dead.png")

# Moving cfr 
can = ROOT.TCanvas(countries[0]+"_cfr", "", 1800,600)
can.Divide(3)
pad = can.cd(1)
pad.SetTopMargin  (0.04)
pad.SetLeftMargin (0.13)
pad.SetRightMargin(0.04)
pad.SetGrid(1,1)
vhc[0].GetXaxis().SetTitleSize(0.05)
vhc[0].GetYaxis().SetTitleSize(0.05)
vhc[0].GetXaxis().SetTimeDisplay(1)
vhc[0].GetXaxis().SetNdivisions(504)
vhc[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
vhc[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
vhc[0].GetYaxis().SetRangeUser(0,1)
vhc[0].GetYaxis().SetTitle("Case fatality ratio (moving average)")
leg4a = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
leg4a.SetFillColor(0)
leg4a.SetFillStyle(0)
leg4a.SetBorderSize(0)
leg4a.SetTextSize(0.04)
for i in range(len(countries)):
    vhc[i].SetMarkerStyle(settings[i][0])
    vhc[i].SetMarkerColor(settings[i][2])
    vhc[i].SetMarkerSize (settings[i][3]*0.8)
    vhc[i].SetLineColor(settings[i][2])
    vhc[i].Draw("PL" if i==0 else "SAME PL")
    leg4a.AddEntry(vhc[i], vh[i].GetTitle(), "P")
leg4a.Draw()

pad = can.cd(2)
pad.SetTopMargin  (0.04)
pad.SetLeftMargin (0.13)
pad.SetRightMargin(0.04)
pad.SetGrid(1,1)
vhef[0].GetXaxis().SetTitleSize(0.05)
vhef[0].GetYaxis().SetTitleSize(0.05)
vhef[0].GetXaxis().SetTimeDisplay(1)
vhef[0].GetXaxis().SetNdivisions(504)
vhef[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
vhef[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
vhef[0].GetYaxis().SetRangeUser(0,round(vhef[0].GetMaximum()*1.2/10)*10)
vhef[0].GetYaxis().SetTitle("Excess fatality (%)")
legb = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
legb.SetFillColor(0)
legb.SetFillStyle(0)
legb.SetBorderSize(0)
legb.SetTextSize(0.04)
for i in range(len(countries)):
    vhef[i].SetMarkerStyle(settings[i][0])
    vhef[i].SetMarkerColor(settings[i][2])
    vhef[i].SetMarkerSize (settings[i][3]*0.8)
    vhef[i].SetLineColor(settings[i][2])
    vhef[i].Draw("PL" if i==0 else "SAME PL")
    legb.AddEntry(vhef[i], vh[i].GetTitle(), "P")
legb.Draw()

pad = can.cd(3)
pad.SetTopMargin  (0.04)
pad.SetLeftMargin (0.13)
pad.SetRightMargin(0.04)
pad.SetGrid(1,1)
vhed[0].GetXaxis().SetTitleSize(0.05)
vhed[0].GetYaxis().SetTitleSize(0.05)
vhed[0].GetXaxis().SetTimeDisplay(1)
vhed[0].GetXaxis().SetNdivisions(504)
vhed[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
vhed[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
vhed[0].GetYaxis().SetRangeUser(0,(vhed[0].GetMaximum()*1.2))
vhed[0].GetYaxis().SetTitle("Excess death")
print "Excess death: "+str(vhed[0].Integral())
leg4c = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
leg4c.SetFillColor(0)
leg4c.SetFillStyle(0)
leg4c.SetBorderSize(0)
leg4c.SetTextSize(0.04)
for i in range(len(countries)):
    vhed[i].SetMarkerStyle(settings[i][0])
    vhed[i].SetMarkerColor(settings[i][2])
    vhed[i].SetMarkerSize (settings[i][3]*0.8)
    vhed[i].SetLineColor(settings[i][2])
    vhed[i].Draw("PL" if i==0 else "SAME PL")
    leg4c.AddEntry(vhed[i], vh[i].GetTitle(), "P")
leg4c.Draw()
can.SaveAs(savedir+countries[0]+"_cfr.png")


# True estimates
can = ROOT.TCanvas(countries[0]+"_proj", "", 1200,600)
can.Divide(2)
#can.Divide(3)
#pad = can.cd(1)
#pad.SetTopMargin  (0.04)
#pad.SetLeftMargin (0.13)
#pad.SetRightMargin(0.04)
#pad.SetGrid(1,1)
#pad.SetLogy(1)
#vh[0].GetXaxis().SetTitleSize(0.05)
#vh[0].GetYaxis().SetTitleSize(0.05)
#vh[0].GetXaxis().SetTimeDisplay(1)
#vh[0].GetXaxis().SetNdivisions(504)
#vh[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
#vh[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
#vh[0].GetYaxis().SetRangeUser(yrange[0],yrange[1])
#vh[0].GetYaxis().SetTitle("Est. true and documented cases")
#leg5a = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
#leg5a.SetFillColor(0)
#leg5a.SetFillStyle(0)
#leg5a.SetBorderSize(0)
#leg5a.SetTextSize(0.04)
#for i in range(len(countries)):
#    vh[i].SetMarkerStyle(settings[i][0])
#    vh[i].SetMarkerColor(settings[i][2])
#    vh[i].SetMarkerSize (settings[i][3]*0.8)
#    vh[i].SetLineColor(settings[i][2])
#    vh[i].Draw("P" if i==0 else "SAME P")
#    vht[i].SetLineColor(settings[i][2])
#    vht[i].SetLineWidth(2)
#    vht[i].Draw("SAME HIST")
#    leg5a.AddEntry(vh[i], vh[i].GetTitle(), "P")
#leg5a.Draw()

pad = can.cd(1)
pad.SetTopMargin  (0.04)
pad.SetLeftMargin (0.13)
pad.SetRightMargin(0.04)
pad.SetGrid(1,1)
vhdoc[0].GetXaxis().SetTitleSize(0.05)
vhdoc[0].GetYaxis().SetTitleSize(0.05)
vhdoc[0].GetXaxis().SetTimeDisplay(1)
vhdoc[0].GetXaxis().SetNdivisions(504)
vhdoc[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
vhdoc[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
vhdoc[0].GetYaxis().SetRangeUser(0,100)
vhdoc[0].GetYaxis().SetTitle("Est. documented cases (%)")
leg5b = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
leg5b.SetFillColor(0)
leg5b.SetFillStyle(0)
leg5b.SetBorderSize(0)
leg5b.SetTextSize(0.04)
for i in range(len(countries)):
    vhdoc[i].SetMarkerStyle(settings[i][0])
    vhdoc[i].SetMarkerColor(settings[i][2])
    vhdoc[i].SetMarkerSize (settings[i][3]*0.8)
    #vhdoc[i].SetLineColor(settings[i][2])
    vhdoc[i].Draw("PL" if i==0 else "SAME PL")
    leg5b.AddEntry(vhdoc[i], vh[i].GetTitle(), "PL")
leg5b.Draw()

pad = can.cd(2)
pad.SetTopMargin  (0.04)
pad.SetLeftMargin (0.13)
pad.SetRightMargin(0.04)
pad.SetGrid(1,1)
pad.SetLogx(1)
dummy = ROOT.TH1D("dummy",";Tests per death;Est. documented cases (%)",10000-1,10,100000)
dummy.GetXaxis().SetTitleSize(0.05)
dummy.GetYaxis().SetTitleSize(0.05)
dummy.GetXaxis().SetRangeUser(10,100000)
dummy.GetYaxis().SetRangeUser(0,100)
vgtest[0].GetXaxis().SetTitle("Tests per death")
vgtest[0].GetYaxis().SetTitle("Est. detected cases (%)")
dummy.Draw()
leg5c = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
leg5c.SetFillColor(0)
leg5c.SetFillStyle(0)
leg5c.SetBorderSize(0)
leg5c.SetTextSize(0.04)
first = 1
for i in range(len(countries)):
    if vgtest[i].GetN()>0:
        vgtest[i].SetMarkerStyle(settings[i][0])
        vgtest[i].SetMarkerColor(settings[i][2])
        vgtest[i].SetMarkerSize (settings[i][3])
        #vgtest[i].SetLineColor(settings[i][2])
        vgtest[i].Draw("SAME P" if first else "SAME P")
        leg5c.AddEntry(vgtest[i], vh[i].GetTitle(), "P")
        first = 0
leg5c.Draw()
can.SaveAs(savedir+countries[0]+"_detection_rate.png")

# true infected estimate end prediction
can = ROOT.TCanvas(countries[0]+"_true", "", 1200,600)
can.Divide(2)
pad = can.cd(1)
pad.SetTopMargin  (0.04)
pad.SetLeftMargin (0.13)
pad.SetRightMargin(0.04)
pad.SetGrid(1,1)
pad.SetLogy(logy2)
vdt[0].GetXaxis().SetTitleSize(0.05)
vdt[0].GetYaxis().SetTitleSize(0.05)
vdt[0].GetXaxis().SetTimeDisplay(1)
vdt[0].GetXaxis().SetNdivisions(504)
vdt[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
vdt[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
vdt[0].GetYaxis().SetRangeUser(yrange3[0],yrange3[1])
vdt[0].GetYaxis().SetTitle("Daily growth of estimated true cases (%)")
leg6a = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
leg6a.SetFillColor(0)
leg6a.SetFillStyle(0)
leg6a.SetBorderSize(0)
leg6a.SetTextSize(0.04)
for i in range(len(countries)):
    for binx in range(1,vdt[i].GetNbinsX()+1): vdt[i].SetBinError(binx,0)
    vdt[i].SetMarkerStyle(settings[i][0])
    vdt[i].SetMarkerColor(settings[i][2])
    vdt[i].SetMarkerSize (settings[i][3]*0.8)
    vdt[i].SetLineColor(settings[i][2])
    vdt[i].Draw("PL" if i==0 else "SAME PL")
    leg6a.AddEntry(vdt[i], vh[i].GetTitle(), "P")
    for fit in vvftpred[i]:
        fit.SetLineColor(settings[i][2])
        fit.SetLineStyle(7)
        fit.SetLineWidth(1)
        fit.Draw("SAME")
leg6a.Draw()

vht[0].GetXaxis().SetTimeDisplay(1)
vht[0].GetXaxis().SetNdivisions(504)
vht[0].GetXaxis().SetTimeFormat("%b %d %F2020-01-01 00:00:00")
vht[0].GetXaxis().SetRangeUser(xrange[0]*86400,xrange[1]*86400)
vht[0].GetYaxis().SetRangeUser(yrange4[0],yrange4[1])
vht[0].GetYaxis().SetTitle("Estimated total true cases")
leg6b = ROOT.TLegend(0.15,y2leg-len(vh)*0.04,0.35,y2leg)
leg6b.SetFillColor(0)
leg6b.SetFillStyle(0)
leg6b.SetBorderSize(0)
leg6b.SetTextSize(0.04)
for i in range(len(countries)):
    leg6b.AddEntry(vh[i], vh[i].GetTitle(), "P")
    if i==0:
        pad = can.cd(2)
        pad.SetTopMargin  (0.04)
        pad.SetLeftMargin (0.13)
        pad.SetRightMargin(0.04)
        pad.SetLogy(1)
        pad.SetGrid(1,1)
    vht[i].GetXaxis().SetTitleSize(0.05)
    vht[i].GetYaxis().SetTitleSize(0.05)
    vht[i].SetMarkerStyle(settings[i][0])
    vht[i].SetMarkerColor(settings[i][2])
    vht[i].SetMarkerSize (settings[i][3])
    vht[i].SetLineColor  (settings[i][2])
    vht[i].Draw("P" if i==0 else "SAME P")
    vt_pred_err[i].SetMarkerStyle(1)
    vt_pred_err[i].SetMarkerColor(0)
    vt_pred_err[i].SetMarkerSize(0)
    vt_pred_err[i].SetFillStyle(settings[i][4])
    vt_pred_err[i].SetFillColor(settings[i][2])
    vt_pred_err[i].Draw("SAME E4")
    vt_pred_nom[i].SetMarkerStyle(settings[i][1])
    vt_pred_nom[i].SetMarkerColor(settings[i][2])
    vt_pred_nom[i].SetMarkerSize (settings[i][3])
    vt_pred_nom[i].SetLineColor  (settings[i][2])
    vt_pred_nom[i].Draw("SAME P")
leg6b.Draw()
can.SaveAs(savedir+countries[0]+"_true.png")


# Predictions
print "Predictions of comfirmed cases up to 28 days from now:"
daystopred = range(1,8) if not worst else [7]
for day in [14,21,28]:
    if xrange[1]>valid_day+day:
        daystopred.append(day)

for i in range(len(countries)):
    for day in daystopred:
        print ("%s - %s: %d  +%d  -%d" % (countries[i], get_date(valid_day+day), int(vh_pred[i].GetBinContent(valid_day+day)),
                                                   int(vh_err[i].GetBinContent(valid_day+day)+vh_err[i].GetBinError(valid_day+day)-vh_pred[i].GetBinContent(valid_day+day)),
                                                   int(vh_pred[i].GetBinContent(valid_day+day)-vh_err[i].GetBinContent(valid_day+day)+vh_err[i].GetBinError(valid_day+day))))

print "Estimate of detection fraction 12 days ago:"
for i in range(len(countries)):
    print ("%s - %s: %.1f%%" % (countries[i], get_date(valid_day-12), int(vhdoc[i].GetBinContent(valid_day-12))))


# SEIR model
# alpha - incubation period
# beta  - infection probability
# gamma - removal probability

for datatype in ["sample","case", "all"]:
    can = ROOT.TCanvas("ALL_doc_"+datatype, "", 600,600)
    can.SetTopMargin  (0.04)
    can.SetLeftMargin (0.13)
    can.SetRightMargin(0.04)
    can.SetGrid(1,1)
    can.SetLogx(1)
    #can.SetLogy(1)
    dummy2 = ROOT.TH1F("dummy2_"+datatype,";Tests per death;Est. detected cases (%)",10000-1,10,100000)
    dummy2.GetXaxis().SetTitleSize(0.05)
    dummy2.GetYaxis().SetTitleSize(0.05)
    dummy2.GetXaxis().SetRangeUser(10,100000)
    #dummy2.GetYaxis().SetRangeUser(0.1,100)
    dummy2.GetYaxis().SetRangeUser(0,100)
    #dummy2.SetMinimum(0.1)
    dummy2.SetMinimum(0)
    dummy2.SetMaximum(100)
    ncountry = 0
    vxmax = []
    vymax = []
    vxlast = []
    vylast = []
    keep = []
    first = 1
    usedcountries = []
    for country in data.keys():
        if country in testdata:
            ncase  = data [country][3].GetBinContent(valid_day)
            ndeath = death[country][3].GetBinContent(valid_day)
            #print country+" "+str(ndeath)+" "+str(ncase)
            if country in preddata and (datatype == testdata[country][2] or datatype == "all") and ndeath>=10:
                ncountry += 1
                gtest = preddata[country][7]
                x = ROOT.Double(0)
                y = ROOT.Double(0)
                xmax = 0
                ymax = 0
                n = 0
                for i in range(gtest.GetN()):
                    gtest.GetPoint(i, x, y)
                    if x>10 and y>0:
                        n += 1
                        if y>ymax:
                            xmax = float(x)
                            ymax = float(y)
                if n>0:
                    usedcountries.append(country)
                    vxmax.append(xmax)
                    vymax.append(ymax)
                    vxlast.append(x)
                    vylast.append(y)
                    gtest_last  = ROOT.TGraph(1)
                    gtest_last.SetPoint(0, x, y)
                    gtest_clean = ROOT.TGraph(n)
                    n = 0
                    for i in range(gtest.GetN()):
                        gtest.GetPoint(i, x, y)
                        if x>10 and y>0:
                            gtest_clean.SetPoint(n, x, y)
                            n += 1
                    gtest_clean.SetMarkerSize(0.1)
                    gtest_clean.SetMarkerStyle(20)
                    gtest_clean.SetMarkerColor(ncountry+1)
                    gtest_clean.SetLineColor(ncountry+1)
                    if first: gtest_clean.SetHistogram(dummy2)
                    gtest_clean.Draw("AL" if first else "SAME L")
                    gtest_last.SetMarkerStyle(20)
                    gtest_last.SetMarkerColor(ncountry+1)
                    gtest_last.SetMarkerSize(1.2)
                    gtest_last.Draw("SAME P")
                    for i in range(len(countries)):
                        if country == countries[i]:
                            gtest_clean.SetLineWidth(3)
                            gtest_clean.SetLineColor(settings[i][2])
                            gtest_last.SetMarkerSize(1.5)
                            gtest_last.SetMarkerColor(settings[i][2])
                            lat = ROOT.TLatex(20, y2leg*100-5*i, vh[i].GetTitle())
                            lat.SetTextColor(settings[i][2])
                            lat.SetTextSize(0.05)
                            lat.Draw("SAME")
                            keep.append(lat)
                    keep.append(gtest_last)
                    keep.append(gtest_clean)
                    first = 0
    #gtest_all_max = ROOT.TGraph(len(vxmax))
    #for i in range(len(vxmax)):
    #    gtest_all_max.SetPoint(i, vxmax[i], vymax[i])
    #gtest_all_max.SetMarkerStyle(20)
    #gtest_all_max.SetMarkerColor(1)
    #gtest_all_max.SetMarkerSize(0.6)
    #gtest_all_max.Draw("SAME P")
    #fitmax = ROOT.TF1("allmax","[0]+exp([1]*log(x))",10,15000)
    #fitmax.SetLineColor(1)
    #gtest_all_max.Fit("allmax","RMWQ")
    gtest_all_last = ROOT.TGraph(len(vxlast))
    for i in range(len(vxlast)):
        gtest_all_last.SetPoint(i, vxlast[i], vylast[i])
    gtest_all_last.SetMarkerStyle(20)
    gtest_all_last.SetMarkerColor(1)
    gtest_all_last.SetMarkerSize(0.6)
    gtest_all_last.Draw("SAME P")
    fitlast = ROOT.TF1("alllast","[0]+exp([1]*log(x))",10,15000)
    gtest_all_last.Fit("alllast","RMWQ")
    for i in range(len(vxlast)):
        ratio = vylast[i] / fitlast.Eval(vxlast[i])
        #if datatype == "all": print ("RATIO: %s - %.2f" % (usedcountries[i], ratio))
    can.SaveAs(savedir+"ALL_detection_rate_"+datatype+".png")
