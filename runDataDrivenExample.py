from os import system,getcwd


fitExe = 'util/myFit.exe'
configFile = 'config/dataDriven.config'

print '####################### start read input histograms ###############################'
system(fitExe + ' h ' + configFile)
print '####################### finished read input histograms ###############################'
print '####################### start prefit plots ###############################'
system(fitExe + ' d ' + configFile)
print '####################### finished prefit plots ###############################'
print '####################### start RooStats workspace ###############################'     
system(fitExe + ' w ' + configFile)
print '####################### finished RooStats workspace ###############################'
print '####################### start fit ###############################'
system(fitExe + ' f ' + configFile)
print '####################### finished fit ###############################'
print '####################### start postfit plots ###############################'     
system(fitExe + ' p ' + configFile)
print '####################### finished postfit plots ###############################'
print '####################### start limit setting ###############################'
system(fitExe + ' l ' + configFile)
print '####################### finished limit setting ###############################'
print '####################### start separation plots ###############################'
system(fitExe + ' a ' + configFile)
print '####################### finished separation plots ###############################'
