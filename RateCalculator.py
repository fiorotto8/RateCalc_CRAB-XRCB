import pandas as pd
import argparse
import numpy as np
import ROOT

argParser = argparse.ArgumentParser()
argParser.add_argument("-Nf", "--NISTfile", help="path to the NIST file of attenuation coeffients", default="./HeCF4_6040_NIST.txt")
def nparr(list):
    return np.array(list, dtype="d")
def graph(x,y,x_string, y_string,name=None, color=4, markerstyle=22, markersize=1,write=True):
        plot = ROOT.TGraph(len(x),  np.array(x  ,dtype="d")  ,   np.array(y  ,dtype="d"))
        if name is None: plot.SetNameTitle(y_string+" vs "+x_string,y_string+" vs "+x_string)
        else: plot.SetNameTitle(name,y_string+" vs "+x_string)
        plot.GetXaxis().SetTitle(x_string)
        plot.GetYaxis().SetTitle(y_string)
        plot.SetMarkerColor(color)#blue
        plot.SetMarkerStyle(markerstyle)
        plot.SetMarkerSize(markersize)
        if write==True: plot.Write()
        return plot
def ImportNIST(file):
    """Import NIST data from a file."""
    df = pd.read_csv(file, sep='\t') # read the csv file into
    return df
def Power_TF1_XRCB():
    func=ROOT.TF1("Power_TF1_XRCB","([0]*x*x)/((x/[1])**[2]+(x/[1])**[3])",0.1,1000)
    func.SetParameters(A,Eb,n1,n2)  # You can change the exponent value as needed
    return func
def Flux_TF1_XRCB(solidangleMult=True):
    if solidangleMult==False: func=ROOT.TF1("Flux_TF1_XRCB","([0])/((x/[1])**[2]+(x/[1])**[3])",0.1,1000)
    else: func=ROOT.TF1("Flux_TF1_XRCB","4*3.14*([0])/((x/[1])**[2]+(x/[1])**[3])",0.1,1000)
    func.FixParameter(0,A)
    func.FixParameter(1,Eb)
    func.FixParameter(2,n1)
    func.FixParameter(3,n2)
    return func
def Flux_Func_XRCB(x,solidangleMult=True):
    return 4*3.14*(A)/((x/Eb)**n1+(x/Eb)**n2)

if __name__ == "__main__":
    args = argParser.parse_args()
    main=ROOT.TFile("main.root","RECREATE")#root file creation
    main.mkdir("NIST_attCoef")
    main.mkdir("XRCB")
    main.mkdir("NuSTAR")

    #read parameters
    file = open("param.txt", "r")
    for string in file:
        exec(string)
    file.close()
    #compute gas denstiy
    gas_density=he_ratio*he_rho+cf4_ratio*cf4_rho

    #Get attenuation coeffiecnt from file (compute with RHO) &
    #compute iteraction probability for Photoeletric and compton from NIST
    NIST_data = ImportNIST(args.NISTfile)
    main.cd("NIST_attCoef")
    TG_AttPhotoeletric=graph(1E3*NIST_data[NIST_data.columns[0]],gas_density*NIST_data["PhotoelectricAbsorption[cm^2/g]"],"Photon Energy[keV]","PhotoelectricAbsorption[cm^{-1}]",name="#mu Photoelectric")
    TG_AttCompton=graph(1E3*NIST_data[NIST_data.columns[0]],gas_density*NIST_data["IncoherentScattering[cm^2/g]"],"Photon Energy[keV]","IncoherentScattering[cm^{-1}]",name="#mu Compton")
    DetSens_Ph=graph(1E3*NIST_data[NIST_data.columns[0]],1-np.exp(-1*side*gas_density*NIST_data["PhotoelectricAbsorption[cm^2/g]"]),"Photon Energy[keV]","Interaction Probability Detector Photoeletric",name="Interaction Prob Photoelectric")
    DetSens_Comp=graph(1E3*NIST_data[NIST_data.columns[0]],1-np.exp(-1*side*gas_density*NIST_data["IncoherentScattering[cm^2/g]"]),"Photon Energy[keV]","Interaction Probability Detector Compton",name="Interaction Prob Compton")

    #Read parameters of teh XRCB from file and create TF1
    main.cd("XRCB")
    file = open("FuncPar_XRCB.txt", "r")
    for string in file:
        exec(string)
    file.close()
    PowerFuncXRCB,FluxFuncXRCB=Power_TF1_XRCB(),Flux_TF1_XRCB()
    PowerFuncXRCB.Write()
    FluxFuncXRCB.Write()

    #Read file with extracted fluxes from CRAB NEBULA NuSTAR
    main.cd("NuSTAR")
    df = pd.read_csv("Extracted_NuStar_power.csv", sep=';',names=["Energy[keV]","CRABPowerFlux[keV^2 (photons/(cm2 s keV))]"]) # read the csv file into
    PowerCRAB=graph(df[df.columns[0]],df[df.columns[1]],df.columns[0],df.columns[1],name="Power Flux CRAB")
    FluxCRAB=graph(df[df.columns[0]],df[df.columns[1]]/(df[df.columns[0]]*df[df.columns[0]]),df.columns[0],"CRABFlux[photons/(cm2 s keV)]",name="Flux CRAB")

    #compute expected diff fluxes and integral fluxes
    main.cd()
    E=np.arange(en_start,en_stop,en_step)
    PhDetFluxCRAB,ComptDetFluxCRAB,PhDetFluxXRCB,ComptDetFluxXRCB=np.empty(len(E)),np.empty(len(E)),np.empty(len(E)),np.empty(len(E))
    for i,e in enumerate(E):
        PhDetFluxCRAB[i]=FluxCRAB.Eval(e,spline=0,option="S")*side*side*DetSens_Ph.Eval(e,spline=0,option="S")
        ComptDetFluxCRAB[i]=FluxCRAB.Eval(e,spline=0,option="S")*side*side*DetSens_Comp.Eval(e,spline=0,option="S")
        PhDetFluxXRCB[i]=FluxFuncXRCB.Eval(e)*side*side*DetSens_Ph.Eval(e,spline=0,option="S")
        ComptDetFluxXRCB[i]=FluxFuncXRCB.Eval(e)*side*side*DetSens_Comp.Eval(e,spline=0,option="S")
    #TGraph CRAB
    TG_PhHitTotCRAB=graph(E,PhDetFluxCRAB,"Energy[keV]","Photoelectric Detector Rate CRAB [photons s^{-1} keV^{-1})]",name="Photoelectric Hit Rate CRAB")
    TG_ComptHitTotCRAB=graph(E,ComptDetFluxCRAB,"Energy[keV]","Compton Detector Rate CRAB [photons s^{-1} keV^{-1})]",name="Compton Hit Rate CRAB")
    TG_TotHitTotCRAB=graph(E,PhDetFluxCRAB+ComptDetFluxCRAB,"Energy[keV]","Total Detector Rate CRAB [photons s^{-1} keV^{-1})]",name="Total Hit Rate CRAB")
    #TGraph XRCB
    TG_PhHitTotXRCB=graph(E,PhDetFluxXRCB,"Energy[keV]","Total Detector Rate XRCB [photons s^{-1} keV^{-1})]",name="Photoelectric Hit Rate XRCB")
    TG_ComptHitTotXRCB=graph(E,ComptDetFluxXRCB,"Energy[keV]","Total Detector Rate XRCB [photons s^{-1} keV^{-1})]",name="Compton Hit Rate XRCB")
    TG_TotHitTotXRCB=graph(E,PhDetFluxXRCB+ComptDetFluxXRCB,"Energy[keV]","Total Detector Rate XRCB [photons s^{-1} keV^{-1})]",name="Total Hit Rate XRCB")

    #Integral Hit Rate Calcualtion
    IntFluxCRAB, IntFluxXRCB=0,0
    for i,e in enumerate(E):
        IntFluxCRAB=IntFluxCRAB+(TG_TotHitTotCRAB.Eval(e,spline=0,option="S")*en_step)
        IntFluxXRCB=IntFluxXRCB+(TG_TotHitTotXRCB.Eval(e,spline=0,option="S")*en_step)
    print("Integral Flux XRCB (Hz):",IntFluxXRCB)
    print("Integral Flux CRAB (Hz):",IntFluxCRAB)






