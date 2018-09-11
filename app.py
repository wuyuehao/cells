import click
import pandas as pd
import numpy as np
import os
import datetime
import time
import matplotlib


@click.command()
@click.option('--ca', default=100000, prompt='Cells Amount (mL-1)', help='The cells amount in mL-1.')
@click.option('--sf', default=0.02, prompt='Sample Flow (mL/min)', help='The Sample Flow in mL/min.')
@click.option('--dt', default=0.10, prompt='Dwell Time (ms)', help='The Dwell time in milli-sec.')
@click.option('--ne', default=11, prompt='Nebulization Efficiency (%)', help='The Nebulization Efficiency in %.')
@click.option('--ds', default=2100, prompt='Detection Sensitivity (cps/µl/L)', help='The Detection Sensitivity in cps/µl/L.')
@click.option('--te_std', default='TRUE', prompt='TE_STD (T/F)', help='The boolean TE_STD value.')
@click.option('--ct', default=85000 , prompt='Cell Threshold (cps)', help='The Cells Threshold in cps.')
@click.option('--it', default=41000 , prompt='Ion threshold (cps)', help='The Ion Threshold in cps.')
@click.option('-d', default=120 , prompt='Duration (s)', help='The duration in sec.')
@click.argument('input')
@click.argument('output')

def run(ca, sf, dt, ne, ds, te_std, ct, it, d, input, output):

    output = output + datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    os.mkdir(output)

    """Calculate cell data."""

    start = time.time()




    print('loading data from '+ input +'...')

    sampleTime = pd.read_csv(input +"/SAMPLE TIME.csv", ";",  skiprows=[1], low_memory=False)
    sampleY = pd.read_csv(input +"/SAMPLE Y.csv", ";",  skiprows=[1], low_memory=False)
    blankTime = pd.read_csv(input +"/BLANK TIME.csv", ";",  skiprows=[1], low_memory=False)
    blankY = pd.read_csv(input +"/BLANK Y.csv", ";",  skiprows=[1], low_memory=False)


    data = pd.concat([sampleTime['Sample 1'], sampleY['Sample 1']],  axis=1)
    blankData = pd.concat([blankTime['Sample 1'], blankY['Sample 1']], axis=1)
    data.columns=['time', 'intensity']
    blankData.columns = ['time', 'intensity']

    count = len(data)

    print('found '+str(count) +' sample data')
    print(data.columns)

    print('processing...')
    ionBlankThreshold = blankData['intensity'].mean()

    data['Ion_Intensity'] = data.apply (lambda row: calcIonIntensity (row, ct, it, ionBlankThreshold),axis=1)

    ionIntensityAvg = data['Ion_Intensity'].mean()


    data['Blank_Substracted'] = data.apply (lambda row: calcBlankSubstracted (row, ct ,ionIntensityAvg),axis=1)

    data['Element_Mass'] = data.apply (lambda row: calcElementMass (row, dt, ds, sf, ne),axis=1)

    cells = []
    sum = 0
    for a in data['Element_Mass']:
        if not np.isnan(a):
            sum = sum + a
        else:
            if(sum > 0):
                cells.append(sum)
                sum = 0

    if sum > 0:
        cells.append(sum)




    detectedCells = len(cells)
    massMean = np.mean(cells)
    massStd = np.std(cells)
    ionicConcentration = ionIntensityAvg / ds


    result={'Ion Blank Threshold': ionBlankThreshold ,'Detected Cells': detectedCells, 'Mass Mean' : massMean, 'Mass Std' :massStd, 'Ionic Concentration':ionicConcentration}

    max = data['Element_Mass'].max()
    min = data['Element_Mass'].min()
    bins100=[]
    bins1000=[]
    for i in range (0, 101):
        bins100.append(min + i*(max-min)/100)


    for i in range (0, 1001):
        bins1000.append(min + i*(max-min)/1000)

    data['Element_Mass_bin_100']  = pd.cut(data['Element_Mass'],bins=bins100,labels=bins100[:-1],include_lowest=True)
    data['Element_Mass_bin_1000']  = pd.cut(data['Element_Mass'],bins=bins1000,labels=bins1000[:-1],include_lowest=True)

    bin100Table = pd.value_counts(data["Element_Mass_bin_100"], sort=False)
    bin1000Table = pd.value_counts(data["Element_Mass_bin_1000"], sort=False)


    print(result)

    print("writing result into directory "+output +'/  ...')
    cellsDf = pd.DataFrame(cells)
    resultDf = pd.DataFrame(list(result.items() ), columns=['key', 'value'])
    data.to_csv(output+'/sample_processed.csv', sep=',', encoding='utf-8')
    blankData.to_csv(output+'/blank.csv', sep=',', encoding = 'utf-8')
    cellsDf.to_csv(output+'/cells.csv', sep=',', encoding='utf-8')
    resultDf.to_csv(output+'/results.csv', sep=',', encoding='utf-8')


    print("generating charts ...")

    bin100Table.plot().get_figure().savefig(output+'/bin100.pdf')
    bin1000Table.plot().get_figure().savefig(output+'/bin1000.pdf')

    end = time.time()

    print("Finished in " + str(end - start) +" sec" )




def calcIonIntensity (row, ct, it, ionBlankThreshold):
    intensity =  row['intensity']
    if intensity == "" or np.isnan(intensity) or intensity > ct :
        return np.nan
    else:
        if intensity > ionBlankThreshold and intensity < it:
            return intensity
    return 0

def calcBlankSubstracted (row, ionThreashold,ionIntensityAvg ):
    intensity =  row['intensity']
    if(intensity > ionThreashold):
        return intensity - ionIntensityAvg
    return np.nan

def calcElementMass (row, dt, ds, sf, ne):

    intensity =  row['Blank_Substracted']
    if(intensity != np.nan):
        return  intensity * dt/1000/ds*sf*ne/60*10000
    return np.nan



if __name__ == '__main__':
    run()
