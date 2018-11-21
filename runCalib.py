import glob
import os
import re
import yaml

# modules=["78","90","88","89",
#          "77","85","32","84"]
# modules=["69","79","67","65",
#          "83","76","36","35",
#          "70","73","86","87"]
# modules=["59","71","55","64",
#          "62","54","44","51"]
modules=["90","88","89",
         "77","85","32","84",
         "69","79","67","65",
         "83","36","35",
         "70","73","44","51",
         "86","87","54","64"
         "55","62","59","71"]
modules=["76"]
eospath='/eos/cms/store/group/dpg_hgcal/tb_hgcal/calibration/'

rootfiles=[]
yamlfiles=[]

for i in modules:
    newpath=eospath+"module"+i
    rootfiles=glob.glob(newpath+"/ana_output/*.root")
    os.system("rm calibrations/Calibration_Module"+i+".root")
    plotDir="plots/module"+i+"/"
    if os.path.isdir(plotDir)!=True:
        os.system("mkdir "+plotDir)
    for j in rootfiles:
        rfname=os.path.basename(j)
        yamlname=re.split('_pedestal',rfname)
        yamlname=newpath+"/yaml/"+yamlname[0]+".yaml"
        with open(yamlname) as fin:
            yamlobj=yaml.safe_load(fin)
            daq_options=yamlobj['daq_options']
            if daq_options['acquisitionType']!='standard':
                os.system("./bin/calibration --fileName="+j+" --module="+str(int(i))+" --channels="+str(daq_options['channelIds'][0]))
            #if daq_options['acquisitionType']!='standard' and daq_options['channelIds'][0]==36:
            #    os.system("./bin/totCalibration --fileName="+j+" --module="+str(int(i))+" --channels="+str(daq_options['channelIds'][0]))            
            #    os.system("./bin/plotter --fileName="+j+" --module="+str(int(i))+" --channels="+str(daq_options['channelIds'][0]))            
