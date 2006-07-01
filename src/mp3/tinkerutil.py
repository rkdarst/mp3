# Richard Darst, 2005
#
#
#

class TinkerLogParser:
    """Class designed to aid in parsing tinker log files.
    
    """

#    MD Step      E Total   E Potential   E Kinetic       Temp       Pres
#
#         1    -4040.8772    -4812.1092    771.2320     298.19     781.16
#         2    -4041.1523    -4812.6403    771.4880     298.29     768.67

#       100    -4058.8395    -4818.9563    760.1168     293.90     -93.76
#
# Average Values for the last   100 out of      100 Dynamics Steps
#
# Simulation Time              0.1000 Picosecond
# Total Energy             -4049.8207 Kcal/mole   (+/-   4.7434)
          
    def __init__(self, log):
        if type(log) == file or ( hasattr(log, "readline") ):
            pass
        else:
            log = file(log, "r")

        # Wind ourselvel up to the "MD Step   ... E Total ... " line
        while True:
            line = log.readline()
            if line.find("    MD Step      E Total   E Potential   E Kinetic       Temp       Pres") != -1:
                break
        # ok, we have two more lines to read.
        log.readline()
        # the next line starts our data.
        self._log = log
        self.mdstep = 0

    def next(self):
        """Bring us to the next dataset in the file.
        """
        log = self._log

        while True:  
            line = log.readline()
            #print self.mdstep, line,
            line = line.split()
            try:
                mdstep = int(line[0])
                if mdstep == self.mdstep +1:
                    break
            except:
                #print line
                #print sys.exc_info()
                pass
        self.mdstep = mdstep
        self.e_total = float(line[1])
        self.e_potential = float(line[2])
        self.e_kinetic = float(line[3])
        self.temperature = float(line[4])
        self.pressure = float(line[5])


        if mdstep % 100 == 0:
            log.readline()
            line = log.readline()
            if line[0:28] != " Average Values for the last":
                print "error"
                raise "We are at an inconsintent place"
            log.readline()
            averages = {}

            line = log.readline().split()
            averages["simulation_time"] = float(line[2])
            line = log.readline().split()
            averages["total_energy"] =     float(line[2])
            averages["total_energy_err"] = float(line[5][:-1])
            line = log.readline().split()
            averages["potential_energy"] =     float(line[2])
            averages["potential_energy_err"] = float(line[5][:-1])
            line = log.readline().split()
            averages["kinetic_energy"] =     float(line[2])
            averages["kinetic_energy_err"] = float(line[5][:-1])
            line = log.readline().split()
            averages["intermolecular"] =     float(line[1])
            averages["intermolecular_err"] = float(line[4][:-1])
            line = log.readline().split()
            averages["temperature"] =     float(line[1])
            averages["temperature_err"] = float(line[4][:-1])
            line = log.readline().split()
            averages["pressure"] =     float(line[1])
            averages["pressure_err"] = float(line[4][:-1])
            line = log.readline().split()
            averages["density"] =     float(line[1])
            averages["density_err"] = float(line[4][:-1])

            #print averages
            log.readline()
            self._averages = averages

    # def mdstep(self): return self._mdstep
    # ... etc
    #
    # We do not need these anymore.  If we needed to do dynamic
    # generation, we should use `property` function.  It will allow
    # calling a function to retrieve an attribute.



if __name__ == "__main__":
    import bz2
    #log = bz2.BZ2File("/home/richard/research/dmf_tinker/automation_tests_1/run14/288/log.bz2")
    log = bz2.BZ2File("/home/richard/research/dmf_tinker/automation_tests_1/run14/288/10/0/log.bz2")
    T = TinkerLogParser(log)
    for i in range(100000):
        T.next()
        print T.mdstep()
