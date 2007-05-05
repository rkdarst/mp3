# Richard Darst, 2006

import mp3.log

def _openFile(filename):
    if type(filename) == file:
        return filename
    else:
        return file(filename, "r")

def getMasses(filename, sortby="atomtype"):
    """Parse MASS lines in a CHARMM .inp file

    Return a dict of dicts referencing MASS lines in the charmm `.inp`
    file.  Refer to the sample file below, along with the field names:

    MASS     1 H      1.00800 H ! polar H
    MASS     2 HC     1.00800 H ! N-ter H
    MASS     3 HA     1.00800 H ! nonpolar H
    MASS     4 HT     1.00800 H ! TIPS3P WATER HYDROGEN
             \- `number`
               \- `atomtype`
                      \- `mass`
                              \- `symbol`

    The dict returned can be keyed by atomtype, number, or symbol, or
    mass.  (NOTE: symbol and mass is not unique, so it wouldn't make
    much sense normally to actually use that.  These will only store
    the last unique field.) The default to sort by is 'atomtype'.  The
    values are dicts containing the remaining keywords, so, if you
    sorted by atomtype, you could do:

    Masses = getmasses('some/file.inp')
    symbol = Masses['HC']['symbol']
    number = Masses['HC']['number']
    mass = Masses['HC']['mass']

    Or to sort by the number field:

    Masses = getmasses('some/file.inp', sortby='number')
    symbol = Masses[1]['symbol']

    """
    mp3.log.debug("getting masses from file %s"%filename)
    filename = _openFile(filename)
    Masses = {}
    import code
    if sortby == "atomtype":
        codeLine = 'Masses[atomtype] = {"number": number, "mass": mass, "symbol": symbol }'
        codeLine = code.compile_command(codeLine)
    elif sortby == "number":
        codeLine = 'Masses[number] = {"atomtype": atomtype, "mass": mass, "symbol": symbol }'
        codeLine = code.compile_command(codeLine)
    elif sortby == "symbol":
        # Note that `symbol` isn't unique... this probably isn't a very good idea
        codeLine = 'Masses[symbol] = {"atomtype": atomtype, "mass": mass, "number": number }'
        codeLine = code.compile_command(codeLine)
    elif sortby == "mass":
        # Note that `mass` isn't unique... this probably isn't a very good idea
        codeLine = 'Masses[mass] = {"atomtype": atomtype, "symbol": symbol, "number": number }'
        codeLine = code.compile_command(codeLine)
    else:
        raise Exception, "Indexing key `%s` is not recognized, see docstring"%sortby
    
    for line in filename:
        if line[:4] == "MASS":
            number, atomtype, mass, symbol = line[4:].split("!")[0].split()
            mass = float(mass)
            number = int(number)
            exec codeLine
    mp3.log.info("getmasses from CHARMM file: returning %s values keyed by %s"%
                 (len(Masses), sortby))
    return Masses

def getResidues(filename):
    """Parse RESI and PRES information in CHARMM .inp files.

    This will return a dictionary containing various useful
    information about the residues in the input file.  We do not make
    a difference between RESI and PRES residues.

    RESI ALA          0.00
    GROUP
    ATOM N    NH1    -0.47  !     |
    ATOM HN   H       0.31  !  HN-N
    ATOM CA   CT1     0.07  !     |     HB1
    ATOM HA   HB      0.09  !     |    /
         \- atomname
              \- atomtype
                      \-change

    `atomname` is unique.  `atomtype` is not.
    

    Demo:

    Residues = getResidues('some/file.inp')

    # The main key is residue name.  RESI and PRES aren't kept separate.
    Residues['ALA']

    # We can get the overall change of the residue
    Residues['ALA']['charge']   # 0.00

    # We can get individual changes on each atom.  This is indexed by
    # atomNAME.
    Residues['ALA']['atomCharge']['CA']   # 0.07

    # This last one is more tricky.  We want a mapping from atomTYPE
    # to atomNAME, except it isn't unique.  Making use of the
    # assumption that ordering stays the same, we can make a
    # sort-of-mapping.  The atomNames attribute helps with that.  It
    # has atomNames in the order they appeared after each atomType,
    # but is padded with None at the end so that we can verify that we
    # come full-circle by rolling the list as we assign values.
    Residues['ALA']['atomNames']['HA']   # ['HB1', 'HB2', 'HB3', None]

    # You can get a list of atom 
    print Residues['ALA']['atomNames']
    --> ['N', 'HN', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O']

    print Residues['ALA']['atomTypes']
    --> ['NH1', 'H', 'CT1', 'HB', 'CT3', 'HA', 'HA', 'HA', 'C', 'O']
    
    
    """
    mp3.log.debug("getting residues from file %s"%filename)
    filename = _openFile(filename)
    Residues = {}
    resName = None
    residuesSeen = [ ]
    for line in filename:
        #print line
        if line[:5] == "RESI " or line[0:5] == "PRES ":
            resName, charge = line[5:].split("!")[0].split()
            charge = float(charge)
            if resName in residuesSeen:
                mp3.log.error("we have already seen residue %s !"%resName)
                raise
            residuesSeen.append(resName)
            Residues[resName] = {"charge": charge,
                                 "atomNamesByType": {},
                                 "atomCharge": {},
                                 "groups": [],
                                 "atomNames": [],
                                 "atomTypes": []}
            R = Residues[resName]
            #print resName
            continue
        if line[:5] == "GROUP":
            R["groups"].append([])
        if line[:5] == "ATOM ":
            atomName, atomType, charge = line[5:].split("!")[0].split()
            charge = float(charge)
            R["atomNamesByType"].setdefault(atomType, [] )
            R["atomNamesByType"][atomType].append(atomName)
            R["atomCharge"][atomName] = charge
            # Some residues/PRESs don't have a leading GROUPS command.
            # Atoms before that are *ignored* in the group list.
            if R["groups"]: R["groups"][-1].append(atomName)
            R["atomNames"].append(atomName)
            R["atomTypes"].append(atomType)
    # append None to each of the atomNames lists
    for R in Residues.itervalues():
        for r2 in R["atomNamesByType"].itervalues():
            r2.append(None)
    mp3.log.info("Parsed %s residues"%(len(Residues)))
    return Residues



def getUBrefList(filename):
    """Parses Urey-Bradley (ub) information in the CHARMM parameter file, par*.inp, which
    is given as the only argument to the function.
    
    Specifically, this function searches the ANGLE section of the CHARMM parameter file
    and finds angles that have ub parameters.  It then creates a tuple comprised of
    the three atom _types_ defining the angle and appends this tuple as well as its inverse
    to a list.  It does this for each ub angle.  This list is to be used as a reference
    list for determining which angles out of a given set should have ub interactions, which
    is a problem that arises, e.g., when converting .psf files to .top files (see psf2top).

    """
    refUBList = []
    fin = file(filename)
    for line in fin:
        if len(line.split()) >= 5:
            if line.split()[4]=='Kub':
                fin.readline()
                break

    for line in fin:
        if line.find('lipids') > 0:
            break
        try:
            float(line.split()[5]),float(line.split()[6])
            cmdExecuted = True
        except:
            cmdExecuted = False
            pass
        if cmdExecuted:
            refUBList.append((line.split()[0],line.split()[1],line.split()[2]))
            refUBList.append((line.split()[2],line.split()[1],line.split()[0]))

    return refUBList
