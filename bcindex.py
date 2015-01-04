'''
This script was designed by Andrew Shatz and Dylan Broderick to
test the consistency of a resampled Species Distribution Model using the Boyce
Continuous Index. This script was originally created on 12/13/2013.

This script was updated with the Idrisi Python framework under development (idrtools/idrpy)
and was released on 1/2/2014.


'''

# Import string and Idrisi COM server
from idrtools import *
from idrtools.idrfiles import *
from idrtools.modules import modules

# Import remaining dependent libraries.
import os
import shutil
import numpy as np
from scipy import stats
from scipy.stats import spearmanr, norm #Edited 10/9/2013
from glob import glob
import csv

# Create idrtools projects and files objects
project = IdrisiExplorer(current_project[0])
idrf = IdrisiFiles(project)

# This class is the Boyce Continuous Object where all 
class bcindex():

    # Define "public" properties attributes
    # Number of classes
    bClass = 100
    # Number of width of window
    bWindow = 10
    # Size of confidence interval
    bConfInt = 90
    # Name of "bias" file
    bias = None

    #Define "private" properties
    #Total number of validation locations
    __totalValid = 0.0

    # Constructor: Set up Boyce Continuous Index project.
    def __init__(self, in_hsm, in_valid, in_mask, out_workspace):
        try:
            print("Initializing Boyce Continuous Index project...")
            #Create hsm file list
            self.hsm = self.__hsmTest(in_hsm)
            #Create valid file list
            self.valid = self.__validTest(in_valid)
            #Assign mask name
            self.mask = self.__maskTest(in_mask)
            #Assign output workspace
            self.out_workspace = self.__outputTest(out_workspace)
            #Assign dump directory path
            self.dump = os.path.join(self.out_workspace, "dump")

            print("Project initialized successfully!")
        except ValueError:
            print("Project failed to initialize due to improper arguments.")


    #Set test for hsm input
    def __hsmTest(self, in_hsm):
        if type(in_hsm) == list:
            for f in in_hsm:
                if os.path.splitext(f)[-1] != '.rst':
                    raise ValueError("Improper file type for habitat suitability model imputs")
            return in_hsm
        elif type(in_hsm) == str:
            if os.path.splitext(in_hsm)[-1] != '.rgf':
                raise ValueError("Improper group file type import.")
            else:
                tmp = (item for item in idrf.ReadRgf(in_hsm))
                return [item + ".rst" for item in tmp if os.path.splitext(item)[-1] == '']


    #Set test for valid input
    def __validTest(self, in_valid):
        if type(in_valid) == list:
            test = 1
            for f in in_valid:
                if os.path.splitext(f)[-1] != '.rst':
                    test = 0
            if test == 0:
                raise ValueError("Improper file type for validation imputs")
            else:
                return in_valid
        elif type(in_valid) == str:
            if os.path.splitext(in_valid)[-1] != '.rgf':
                raise ValueError("Improper group file type import.")
            else:
                tmp = (item for item in idrf.ReadRgf(in_valid))
                return [item + ".rst" for item in tmp if os.path.splitext(item)[-1] == '']

    #Set test for mask input
    def __maskTest(self, in_mask):
        if type(in_mask) == str:
            if os.path.splitext(in_mask)[-1] != '.rst':
                raise ValueError("Improper file type for mask imputs")
            else:
                return in_mask

    #Set test for output workspace
    def __outputTest(self, out_workspace):
        #Check to see what data type the input is
        if type(out_workspace) != str:
            raise ValueError("Output workspace parameter must be a string.")
        #Check to see if the directory exists
        if os.path.exists(out_workspace) == False:
            print("Output directory does not exist.  Creating...")
            os.mkdir(out_workspace)
        return out_workspace

    # This method find the max data of an input raster or
    # vector file using a documentation lookup method
    def __getMaxValue(self, input):
        doc = Documentation(input)
        return float(doc.MaxValue())

    # Function to produce boyce threshold image
    def BoyceImage(self, hsmImage, bclass, mask):
        # Calculate maximum HSM value from raster
        maximum = self.__getMaxValue(hsmImage)
        
        # Create Boyce image reclass file
        # Output name
        rclfile = 'Boyce_' + str(bclass) + '.rcl'
        # Create custom reclass list
        width = maximum/int(bclass)
        reclassList = []
        reclassList.append([0, -999.0, 0.0]) #Remove background values
        reclassList.append([1, 0.0, width]) #Default first row
        for val in range(2, bclass): #Rows 2 to n-1
            low = (val-1) * width
            high = low + width
            reclassList.append([val, low, high])
        reclassList.append([bclass, maximum - width, 999])
        # Write reclass file
        idrf.WriteRcl(rclfile, self.dump, reclassList)

        # Run reclass module to produce Boyce Class image draft
        dumpfile = os.path.join(self.dump, "tmp001")
        modules.reclass("rst", hsmImage, dumpfile, 1, 3, os.path.join(self.dump, rclfile))
        
        # Overaly mask to produce final Boyce Class image
        bclassImage = os.path.join(self.dump,'BClass_'+str(bclass))
        modules.overlay(dumpfile, mask, bclassImage, 3)

        # Return file path name for later consumption
        return bclassImage

    # Convert extracted predicted and expected frequencies extracted 
    # as .avl files to list objects
    def __PEreadAvl(self, bclassImage, inputImage):
        print "inputImage: " + inputImage
        output =  os.path.join(self.dump, os.path.splitext(inputImage)[0] + '.avl')
        print "output: " + output
        modules.extract(bclassImage, inputImage, 3, 1, output)
        return idrf.ReadAvl(os.path.basename(output))

    # 
    def peCalculation(self, hsmImage, validation, mask, biasfile, bclass, window, peImage):
        #Declare inital variables
        peList = [-999]*int(bclass)
        totalVal = 0.0
        totalArea = 0.0

        # Calculate predicted and expected frequencies from input imagery
        bclassImage = self.BoyceImage(hsmImage, bclass, mask)
        print "BclassImage: " + bclassImage
        hsmPF = self.__PEreadAvl(bclassImage, validation)
        if biasfile != None:
            hsmEF = self.__PEreadAvl(bclassImage, biasfile)
        else:
            hsmEF = self.__PEreadAvl(bclassImage, mask)
        for i in range(int(bclass)):
            totalVal = totalVal + hsmPF[i]
            totalArea = totalArea + hsmEF[i]
        for i in range(int(bclass)):
            hsmPF[i] = hsmPF[i]/totalVal
            hsmEF[i] = hsmEF[i]/totalArea

        # Caculate the BCI across all windowed predicted/expected frequencies
        wpeList = []
        for i in range(int(bclass)-int(window)+1):
            wpf = 0.0
            wef = 0.0
            for j in range(i, i + int(window)):
                wpf = wpf + hsmPF[j]
                wef = wef + hsmEF[j]
            wpeList.append(wpf/wef)
        for i in range(int(bclass)):
            if i < int(window)-1:
                delpe = wpeList[0:i+1]
                peList[i] = np.mean(delpe)
            if int(window)-1 <= i < int(bclass)\
               -int(window)+1:
                delpe = wpeList[i-(int(window)-1):i+1]
                peList[i] = np.mean(delpe)
            if int(bclass)-int(window)+1 <= i < int(bclass):
                delpe = wpeList[i-int(window)+1:int(bclass)-int(window)+1]
                peList[i] = np.mean(delpe)

        # Create final BCI Spatial interpretations using reclass module from
        # custom reclass parameter file
        rclfilename = 'PE_'+str(bclass)+'.rcl'
        rclList = [[item, i+1, i+2] for i, item in enumerate(peList)] # Added 12/30/2014
        idrf.WriteRcl(rclfilename, self.dump, rclList) # Added 12/30/2014
        modules.reclass('rst', bclassImage, peImage, 2, 3, rclfilename) # Added 12/30/2014

        # Add validation location count to class total
        self.__totalValid += totalVal
        return peList, peImage

    # Calculate the summary statistics for the .csv and .txt file summary outputs.
    def sumStats(self, pelistgroup, bclass, confint):
        outputMatrix = [["Boyce Class", "Mean", "Median", "Range", "Lower Bound", "Upper Bound"]]
        convertMatrix = []
        header = [["Mean", "Median", "Range", "Lower Bound", "Upper Bound"]]
        for i in range(len(pelistgroup[0])):
            tempList = []
            for list in pelistgroup:
                tempList.append(list[i])
            convertMatrix.append(tempList)
        for i in range(len(convertMatrix)):
            average = np.mean(convertMatrix[i])
            median = np.median(convertMatrix[i])
            convertMatrix[i].sort()
            ran = convertMatrix[i][-1] - convertMatrix[i][0]
            sterr = np.std(convertMatrix[i])/(len(convertMatrix[i])**0.5)
            alpha = float(confint)/100
            lbound = stats.norm.interval(alpha, average, sterr)[0]
            if lbound < 0:
                lbound = 0
            ubound = stats.norm.interval(alpha, average, sterr)[1] #Edited 10/9/2013
            outputMatrix.append([str(i+1), average, median, ran, lbound, ubound])
        return convertMatrix, outputMatrix

    # Write both summary text file and csv file for graphing.  The latter graph can be plotted using
    # either microsoft excel or R.
    def writeOutputFiles(self, prefix, bciList, totalPO, bciOverall, bclass, sumstats):
        txtOutput = os.path.join(project.workdir, prefix + ' Boyce Continuous Index Results.txt')
        writeTxt = open(txtOutput,'w')
        writeTxt.write("Boyce Continuous Index Results\n")
        writeTxt.write("\n")
        writeTxt.write("%-20s %9s\n" % ('K-Fold Partitions: ', str(len(bciList))))
        writeTxt.write("%-20s %9s\n" % ('Presence Only Sites:', str(totalPO)))
        writeTxt.write("\n")
        writeTxt.write("***************************************************\n")
        writeTxt.write("\n")
        writeTxt.write("%-20s %20s\n" % ('K-Fold Partition', 'Boyce Continuous Index'))
        writeTxt.write("\n")
        for i in range(len(bciList)):
            if bciList[i] > 0:
                writeTxt.write("%-20s %22.4s\n" % ('Partition ' + str(i+1)+':', bciList[i]))
            elif bciList[i] < 0:
                writeTxt.write("%-20s %22.5s\n" % ('Partition ' + str(i+1)+':', bciList[i]))
        writeTxt.write("\n")
        writeTxt.write("***************************************************\n")
        writeTxt.write("\n")
        if bciOverall > 0:
            writeTxt.write("%-20s %22.5s\n" % ('Overall Boyce Index:', str(bciOverall)))
        elif bciOverall < 0:
            writeTxt.write("%-20s %22.5s\n" % ('Overall Boyce Index:', str(bciOverall)))
        #Qualitative model assessment
        if bciOverall <= -0.5:
            writeTxt.write("%-20s %22s" % ('Model Quality:', 'Very Poor'))
        if -0.5 < bciOverall < -0.1:
            writeTxt.write("%-20s %22s" % ('Model Quality:', 'Poor'))
        if -0.1 <= bciOverall <= 0.1:
            writeTxt.write("%-20s %22s" % ('Model Quality:', 'Random'))
        if 0.1 < bciOverall < 0.5:
            writeTxt.write("%-20s %22s" % ('Model Quality:', 'Good'))
        if bciOverall >= 0.5:
            writeTxt.write("%-20s %22s" % ('Model Quality:', 'Very Good'))
        writeTxt.close()
        csvOutput = os.path.join(project.workdir, prefix + " Boyce Continuous Index Chart.csv")
        writeCsv = open(csvOutput, 'wb')
        writer = csv.writer(writeCsv)
        for row in sumstats:
            del(row[3])
            writer.writerow(row)
        writeCsv.close()

    # The main function of the bcindex class to perform the full geoprocessing workflow.
    # Check all class attributes to ensure that the input parameters are desired.
    def BoyceContinuousIndex(self, out_prefix):
        #Start the process, create initial empty list objects
        bciResults = []
        peListFull = []
        outputlist = []
        maxlist = []

        #Create the dump directory (remove directory if exists)
        if os.path.exists(self.dump) != True:
            os.mkdir(self.dump)

        #Add dump directory to resource folder
        project.remAllResourceDir()
        if self.dump not in project.dirlist:
            project.addResourceDir(self.dump)
        
        # Perform BCI calculation for all hsm/validation images
        # Produce BCI Spatial interpretations.
        for i, hsm_valid in enumerate(zip(self.hsm, self.valid)):
            hsmClass = os.path.join(self.dump, "HSM_B"+str(self.bClass)+'_'+str(i+1))
            peFileName = 'BCI_PE_'+str(i+1)
            outputPE = out_prefix+'_'+peFileName+'_'+str(self.bClass)+'_'+str(self.bWindow) + ".rst"
            peList = self.peCalculation(hsm_valid[0], hsm_valid[1], self.mask, 
                                        self.bias, self.bClass, self.bWindow, outputPE)[0]
            peListFull.append(peList)
            outputlist.append(outputPE)
            bciResults.append(spearmanr(peList, range(1,self.bClass+1))[0])
            tempmax = self.__getMaxValue(outputPE)
            maxlist.append(tempmax)
        maxlist.sort()
        maximum = str(maxlist[-1])
        rgfname = out_prefix+'_BCI_PE_'+str(self.bClass)+'_'+str(self.bWindow)+'.rgf'
        idrf.WriteRgf(project.workdir, outputlist, rgfname)   
        
        # Produce mean and median BCI raster images
        outname = os.path.splitext(rgfname)[0]
        modules.tstats(rgfname, outname, [True, True, False], 0, maximum, 0, len(self.hsm))
        api.DisplayFile(outname+"_Mean", 'quant')
        api.DisplayFile(outname+"_Median", 'quant')

        # Produce summary .csv and .txt files
        bciOverall = np.mean(bciResults)
        sumstatslist  = self.sumStats(peListFull, self.bClass, self.bConfInt)
        self.writeOutputFiles(out_prefix, bciResults, self.__totalValid, bciOverall, self.bClass, sumstatslist[1])
        
        # Wind things down
        project.remResourceDir(self.dump)
        self.__totalValid = 0.0

