# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 12:52:27 2023

@author: Lukass Kellijs
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os 

ID = 1 

dirname = "Split" # directory with split files
dt = 5 * 60 # time between data points

# flags
SPLIT = False
MONTH_PLOTS = False
MONTH_PHASE_PLOTS = True


''' Process Files '''
if SPLIT:
    panelPosYStart = 0
    panelNegXStart = 0
    panelPosXStart = 0
    panelNegYStart = 0
    allPanelStart = 0
    
    # read all lines of the file into list
    lines = list()
    with open("data{}.csv".format(ID)) as file:
        lines = file.readlines()
        
    i = 0
    for line in lines:
        stripped = line.strip()
        if stripped == "\"Panel+Y\"":
            panelPosYStart = i
        elif stripped == "\"Panel-X\"":
            panelNegXStart = i
        elif stripped == "\"Panel+X\"":
            panelPosXStart = i 
        elif stripped == "\"Panel-Y\"":
            panelNegYStart = i
        elif stripped == "\"All Solar Panel Groups\"":
            allPanelStart = i
        i += 1
    
    
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    
    # write all individual files
    with open(dirname + "/panel+Y.csv", "w") as file:
        file.writelines(lines[panelPosYStart + 2 : panelNegXStart-1])
    
    with open(dirname + "/panel-X.csv", "w") as file:
        file.writelines(lines[panelNegXStart + 2 : panelPosXStart-1])
    
    with open(dirname + "/panel+X.csv", "w") as file:
        file.writelines(lines[panelPosXStart + 2 : panelNegYStart-1])
        
    with open(dirname + "/panel-Y.csv", "w") as file:
        file.writelines(lines[panelNegYStart + 2 : allPanelStart-1])
        
    with open(dirname + "/all.csv", "w") as file:
        file.writelines(lines[allPanelStart + 2 : -2])
    

''' Process Data ''' 
panelPosXData = pd.read_csv(dirname + "/panel+X.csv", sep=',')   
panelNegXData = pd.read_csv(dirname + "/panel-X.csv", sep=',')   
panelPosYData = pd.read_csv(dirname + "/panel+Y.csv", sep=',')
panelNegYData = pd.read_csv(dirname + "/panel-Y.csv", sep=',')   
allPanelData = pd.read_csv(dirname + "/all.csv", sep=',')

# extract power column
panelPosX = panelPosXData["Power (W)"].to_list()
panelNegX = panelNegXData["Power (W)"].to_list()
panelPosY = panelPosYData["Power (W)"].to_list()
panelNegY = panelNegYData["Power (W)"].to_list()
allPanel = allPanelData["Power (W)"].to_list()

# extract time column
time = allPanelData["Time (UTCG)"].to_list()

# extract intensity column
intensity = allPanelData["Solar Intensity"].to_list()

# create index list
i = range(len(allPanel)) # index list

# create day list for time interpretation (every element is time elapsed in days up to this index)
conversionFactor = dt / 86400 
days = [x * conversionFactor for x in i]


print("Read all files, line length is:", len(allPanel))
print("Simulation time is:", len(allPanel) * dt, "seconds =", round(len(allPanel) * dt / 86400), "days")

# find all month indexes for all panel power
months = {
    "Jan" : 0,
    "Feb" : 0,
    "Mar" : 0,
    "Apr" : 0,
    "May" : 0,
    "Jun" : 0,
    "Jul" : 0,
    "Aug" : 0,
    "Sep" : 0,
    "Oct" : 0,
    "Nov" : 0,
    "Dec" : 0    
}

# populate dictionary with time indexes for each month
for month in months:
    # find first index with new month (either at 00 or at 05 minutes)
    try:
        months[month] = time.index("1 {} 2024 00:00:00.000".format(month))
    except:
        months[month] = time.index("1 {} 2024 00:05:00.000".format(month))


''' Analyse Data ''' 
totalTime = len(allPanel) * dt

# calculate average solar power generation over 1 year for all panels
energyListAll = [x * dt for x in allPanel]
avg1Year = sum(energyListAll) / totalTime
print("Average Power over 1 year:", "%.3f" % avg1Year, "W")

# calculate average solar power generation over each month for all panels
keys = list(months.keys())
monthAverages = list()

for i in range(len(keys)):
    if not i == 11:
        thisMonth = energyListAll[months[ keys[i] ] : months[ keys[i+1] ] ]
        timeElapsed = len(thisMonth) * dt
        avg = sum(thisMonth) / timeElapsed
    else:
        thisMonth = energyListAll[months[keys[i]]:]
        timeElapsed = len(thisMonth) * dt
        avg = sum(thisMonth) / timeElapsed
    monthAverages.append(avg)
    
print("Range of averages over one month:", "%.3f" % min(monthAverages) + "-" + "%.3f" % max(monthAverages), "W")


# find: 
#   amount of sun phases in the simulation period
#   average power generated for each sun phase
#   duration of each sun phase
#   list index of each sun phase
previous = 0 # assume start out of phase
numPhases = 0
    
sunPhaseAverages = list()
sunPhaseDurations = list() # in minutes
sunPhaseIndexes = list() # list of phase starting indexes - last element of list should be disregarded
sunPhaseEnergies = list() # list of energy collected over sun phase
sunPhaseIndexes.append(0)

sumOfEnergies = 0
phaseLength = 0
phaseStart = 0 # stores starting index of phase

for i, el in enumerate(intensity):
    
    sumOfEnergies += energyListAll[i] # for each datapoint get corresponding energy (datapoints outside of sun phase are 0)
    
    if el == 0 and previous == 1: # end of previous sun phase
        # process previous phase
        sunPhaseAverages.append(sumOfEnergies / (phaseLength * dt)) # add average power over finished sun phase
        sunPhaseDurations.append((phaseLength * dt) / 60.0)
        sunPhaseEnergies.append(sumOfEnergies / 3600)
        sunPhaseIndexes.append(phaseStart)
        
        
        # reset phase parameters
        sumOfEnergies = 0
        phaseLength = 0
        
        numPhases += 1
        previous = 0
        
    if el == 1 and previous == 0: # new phase 
        phaseStart = i
        phaseLength += 1
        previous = 1
        
    elif el == 1 and previous == 1: # in phase
        phaseLength += 1
        
        
        
print("Range of sun phase power averages:", "%.3f" % min(sunPhaseAverages[1:-1]) + "-" + "%.3f" % max(sunPhaseAverages[1:-1]), "W")
print("Range of sun phase durations:", "%.3f" % min(sunPhaseDurations[1:-1]) + "-" + "%.3f" % max(sunPhaseDurations[1:-1]), "min")
print("Range of sun phase energies:", "%.3f" % min(sunPhaseEnergies[1:-1]) + "-" + "%.3f" % max(sunPhaseEnergies[1:-1]), "Wh")

averageMaxIndex = sunPhaseAverages.index(max(sunPhaseAverages[1:-1]))
averageMinIndex = sunPhaseAverages.index(min(sunPhaseAverages[1:-1]))
durationMaxIndex = sunPhaseDurations.index(max(sunPhaseDurations[1:-1]))
durationMinIndex = sunPhaseDurations.index(min(sunPhaseDurations[1:-1]))

print("Average power max index is", averageMaxIndex, "at time", time[sunPhaseIndexes[averageMaxIndex]])
print("Average power min index is", averageMinIndex, "at time", time[sunPhaseIndexes[averageMinIndex]])     
print("Duration max index:", durationMaxIndex, "at time", time[sunPhaseIndexes[durationMaxIndex]]) 
print("Duration min index:", durationMinIndex, "at time", time[sunPhaseIndexes[durationMinIndex]])


'''
# detect and remove outliers from sun phase lists - "normalize" sun phases
upThresholdDuration = 100 # threshold duration in minutes
lowThresholdDuration = 80

normSunPhaseAverages = list()
normSunPhaseDurations = list()
normSunPhaseIndexes = list()
print("While normalizing sun phase lists removed:")
for i, el in enumerate(sunPhaseDurations):
    if el > upThresholdDuration:
        print("\t [Upper Duration] Element", el, "at index", i, "at time", time[sunPhaseIndexes[i]])
        continue
    elif el < lowThresholdDuration:
        print("\t [Lower Duration] Element", el, "at index", i, "at time", time[sunPhaseIndexes[i]])
        continue
        
        
    normSunPhaseDurations.append(sunPhaseDurations[i])
    normSunPhaseAverages.append(sunPhaseAverages[i])
    normSunPhaseIndexes.append(sunPhaseIndexes[i])

print("Range of normalised sun phase power averages:", "%.3f" % min(normSunPhaseAverages[1:-1]) + "-" + "%.3f" % max(normSunPhaseAverages[1:-1]), "W")
print("Range of normalised sun phase durations:", "%.3f" % min(normSunPhaseDurations[1:-1]) + "-" + "%.3f" % max(normSunPhaseDurations[1:-1]), "min")
'''


# find average sun phase durations for each month
monthSunPhaseDurationAverages = list()

for i in range(len(keys)):
    # get list of all sun phase durations for this month
    thisMonthDurations = list()
    numSunPhases = 0
    
    if not i == 11:
        # populate the sun phase durations list for this month
        monthStart = months[keys[i]] # start index for this month
        monthEnd = months[keys[i+1]] # end index for this month
        for i, el in enumerate(sunPhaseIndexes):
            if el > monthStart and el < monthEnd: # if the sun phase is within this month
                thisMonthDurations.append(sunPhaseDurations[i])
                numSunPhases += 1
    else:
        # populate the sun phase durations list for this month
        monthStart = months[keys[i]] # start index for this month
        for i, el in enumerate(sunPhaseIndexes[:-1]):
            if el > monthStart: # if the sun phase is within this month
                thisMonthDurations.append(sunPhaseDurations[i])
                numSunPhases += 1
               
    # calculate average duration for this month
    monthSunPhaseDurationAverages.append( sum(thisMonthDurations) / numSunPhases )

print("Range of sun phase duration averages over one month:", "%.3f" % min(monthSunPhaseDurationAverages) + "-" + "%.3f" % max(monthSunPhaseDurationAverages), "min")


# find a sun phase for each month
monthSunPhaseIntervals = list() 
thisMonthPhaseDuration = 0

for i in range(len(keys)):
    # get list of all sun phase durations for this month
    thisMonthPhase = list()
    
    if not i == 11:
        # find first sun phase
        monthStart = months[keys[i]] # start index for this month
        monthEnd = months[keys[i+1]] # end index for this month
        for i, el in enumerate(sunPhaseIndexes):
            if el > monthStart and el < monthEnd: # if the sun phase is within this month
                thisMonthPhaseDuration = int(sunPhaseDurations[i] * 60 / dt) # convert from minutes to index
                monthSunPhaseIntervals.append( (el, el + thisMonthPhaseDuration) )
                
    else:
        # populate the sun phase durations list for this month
        monthStart = months[keys[i]] # start index for this month
        for i, el in enumerate(sunPhaseIndexes[:-1]):
            if el > monthStart: # if the sun phase is within this month
                thisMonthPhaseDuration = int(sunPhaseDurations[i] * 60 / dt) # convert from minutes to index
                monthSunPhaseIntervals.append( (el, el + thisMonthPhaseDuration) )


''' Plotting ''' 
if not os.path.exists("Plots"):
    os.mkdir("Plots")
    
# plotting range of averages over one month
col = ['C{}'.format(x) for x in range(12)]
fig, ax = plt.subplots(figsize=(12,8))
plt.grid()
plt.bar(keys, monthAverages, color=col)
ax.set_ylabel("Power (W) All Panels")
ax.set_xlabel("Month")
plt.title("Average Power (W) for All Panels Over Year")
plt.savefig("Plots/avgyear.png", dpi=300)

# plotting range of power averages over sun phases
fig, ax = plt.subplots(figsize=(12,8))
plt.grid()
plt.plot(range(len(sunPhaseAverages)), sunPhaseAverages, color="b")
ax.set_ylabel("Average Power of Sun Phase")
ax.set_xlabel("Sun Phase Index")
plt.title("Average Power of Sun Phase for All Sun Phases")
plt.savefig("Plots/sunphaseavgs.png", dpi=300)

# plotting range of durations of sun phases
fig, ax = plt.subplots(figsize=(12,8))
plt.grid()
plt.plot(sunPhaseIndexes[:-1], sunPhaseDurations, color="darkorchid")
ax.set_ylabel("Duration of Sun Phase")
ax.set_xlabel("Sun Phase Index")
plt.title("Duration of Sun Phase for All Sun Phases")
plt.savefig("Plots/sunphasedurations.png", dpi=300)

# plotting averages for sun phase duration over one month
col = ['C{}'.format(x) for x in range(12)]
fig, ax = plt.subplots(figsize=(12,8))
plt.grid()
plt.bar(keys, monthSunPhaseDurationAverages, color=col)
ax.set_ylabel("Average Duration of Sun Phases in Month (min)")
ax.set_xlabel("Month")
plt.title("Average Duration of Sun Phases Over Year")
plt.savefig("Plots/sunphasedurationsyear.png", dpi=300)

# histogram of energy collected over one sun phase, for all sun phases in dataset
n_bins = 100
bin_range = (0, 10)
fig, ax = plt.subplots(figsize=(12,8))
plt.grid()
N, bins, patches = plt.hist(sunPhaseEnergies, bins = n_bins, range = bin_range)
# color bins
fracs = N / N.max() 
norm = colors.Normalize(fracs.min(), fracs.max())
for thisfrac, thispatch in zip(fracs, patches):
    color = plt.cm.plasma(norm(thisfrac))
    thispatch.set_facecolor(color)
    
ax.set_ylabel("Number of Sun Phases")
ax.set_xlabel("Energy Collected in Sun Phase (Wh)")
plt.title("Histogram of Energy Collected over One Sun Phase")
#plt.savefig("Plots/sunphasehistogram.png", dpi=300)

# plotting power graphs for each month
if MONTH_PLOTS:
    print("Plotting month graphs...")
    if not os.path.exists("Plots/Months"):
        os.mkdir("Plots/Months")
        
    for i in range(len(keys)):
        if not i == 11:
            thisMonth = keys[i]
            nextMonth = keys[i+1]
            fig, ax = plt.subplots(figsize=(16,8))
            plt.grid()
            plt.plot(days[months[thisMonth] : months[nextMonth]], allPanel[months[thisMonth] : months[nextMonth]], color="m", linewidth=0.1, label = "All Panels (W)")
            plt.plot(days[months[thisMonth] : months[nextMonth]], panelPosX[months[thisMonth] : months[nextMonth]], color="g", linewidth=0.1, label = "+X Panel (W)")
            plt.plot(days[months[thisMonth] : months[nextMonth]], panelNegX[months[thisMonth] : months[nextMonth]], color="b", linewidth=0.1, label = "-X Panel (W)")
            plt.plot(days[months[thisMonth] : months[nextMonth]], panelPosY[months[thisMonth] : months[nextMonth]], color="y", linewidth=0.1, label = "+Y Panel (W)")
            plt.plot(days[months[thisMonth] : months[nextMonth]], panelNegY[months[thisMonth] : months[nextMonth]], color="darkorange", linewidth=0.1, label = "-Y Panel (W)")
            ax.set_ylabel("Power (W)")
            ax.set_xlabel("Day")
            plt.title("Power (W) All Panels vs Day in {}".format(thisMonth))
            # edit legend to make it visible
            plt.legend()
            leg = ax.legend()
            for line in leg.get_lines():
                line.set_linewidth(4.0)
            plt.savefig("Plots/Months/{}.png".format(thisMonth), dpi=300)
        else:
            thisMonth = keys[i]
            fig, ax = plt.subplots(figsize=(16,8))
            plt.grid()
            plt.plot(days[months[thisMonth] : ], allPanel[months[thisMonth] : ], color="m", linewidth=0.1, label = "All Panels (W)")
            plt.plot(days[months[thisMonth] : ], panelPosX[months[thisMonth] : ], color="g", linewidth=0.1, label = "+X Panel (W)")
            plt.plot(days[months[thisMonth] : ], panelNegX[months[thisMonth] : ], color="b", linewidth=0.1, label = "-X Panel (W)")
            plt.plot(days[months[thisMonth] : ], panelPosY[months[thisMonth] : ], color="y", linewidth=0.1, label = "+Y Panel (W)")
            plt.plot(days[months[thisMonth] : ], panelNegY[months[thisMonth] : ], color="darkorange", linewidth=0.1, label = "-Y Panel (W)")
            ax.set_ylabel("Power (W)")
            ax.set_xlabel("Day")
            plt.title("Power (W) Panels vs Day in {}".format(thisMonth))
            # edit legend to make it visible
            plt.legend()
            leg = ax.legend()
            for line in leg.get_lines():
                line.set_linewidth(4.0)
            plt.savefig("Plots/Months/{}.png".format(thisMonth), dpi=300)       
    
    print("Done")
    
# plotting power graphs for each month
if MONTH_PHASE_PLOTS:
    print("Plotting month sun phases...")
    if not os.path.exists("Plots/Month Phases"):
        os.mkdir("Plots/Month Phases")
        
    for i in range(len(keys)):
        thisMonth = keys[i]
        start, end = monthSunPhaseIntervals[i]
        
        
        
        phaseTime = [t * dt / 60 for t in range(end - start)] # set x axis to time in minutes
        
        fig, ax = plt.subplots(figsize=(12,8))
        plt.grid()
        
        plt.plot(phaseTime, allPanel[start : end], color="m", linewidth = 2, label = "All Panels (W)")
        plt.plot(phaseTime, panelPosX[start : end], color="g", linewidth = 2, label = "+X Panel (W)")
        plt.plot(phaseTime, panelNegX[start : end], color="b", linewidth = 2, label = "-X Panel (W)")
        plt.plot(phaseTime, panelPosY[start : end], color="y", linewidth = 2, label = "+Y Panel (W)")
        plt.plot(phaseTime, panelNegY[start : end], color="darkorange", linewidth = 2, label = "-Y Panel (W)")
        
        ax.set_ylabel("Power (W)")
        ax.set_xlabel("Time from Sun Phase Start (min)")
        plt.title("Power (W) All Panels vs Time for Sun Phase in {}".format(thisMonth))
        # edit legend to make it visible
        plt.legend()
        leg = ax.legend()
        for line in leg.get_lines():
            line.set_linewidth(4.0)
        plt.savefig("Plots/Month Phases/{}.png".format(thisMonth), dpi=300)
    
    print("Done")


# additional plots
'''
# plotting range of durations of normalized sun phases
fig, ax = plt.subplots(figsize=(12,8))
plt.grid()
plt.plot(normSunPhaseIndexes, normSunPhaseDurations, color="darkorchid")
ax.set_ylabel("Duration of Sun Phase")
ax.set_xlabel("Sun Phase Index")
plt.title("Duration of Sun Phase for All Normalized Sun Phases")
'''

'''
# plotting power of all panels in specific month
fig, ax = plt.subplots(figsize=(12,8))
plt.grid()
plt.plot(days[months["Jan"] : months["Feb"]], allPanel[months["Jan"] : months["Feb"]], color="m", lineWidth=0.1)
ax.set_ylabel("Power (W) All Panels")
ax.set_xlabel("Day")
plt.title("Power (W) All Panels vs Index in {}".format("Jan"))
#plt.savefig("all.png", dpi=300)
'''

'''
# plotting specific subset of indexes of all panel power
start = 0
end = 120

fig, ax = plt.subplots(figsize=(12,8))
plt.grid()
plt.plot(i[start:end], allPanel[start:end], color="m")
ax.set_ylabel("Power (W) All Panels")
ax.set_xlabel("Index")
plt.title("Power (W) All Panels vs Index")
#plt.savefig("all.png", dpi=300)
'''