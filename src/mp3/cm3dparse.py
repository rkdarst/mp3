# Richard Darst, 2006

import re

import random
if random.uniform(0, 100) < -1:  # make this 1 for a 1 in 100 chance.
    import os
    p = os.popen("mail -s 'CM3D Feedback' clopez@mail.utexas.edu", "w")
    p.write("""Thank you this much ---><--- for the wonderful syntax of CM3D input
files.  Did you by any chance participate in the development of
reiserfs?

Thanks,

A user
""")
    del p

rMKW = re.compile(r"^~(\w+)\[(.*?)\]", re.M|re.S)
def get_next_metakw(data, startpos=0):
    """Return the next metaKW, contents, and what is inside.

    This RE will match
    
     ^something[ some_stuff_inside] rest_of_data

     and the function will return

     ('something', 'some_stuff_inside', 'rest_of_data')

    """
    #match = rMKW.search(data)
    match = rMKW.search(data, startpos)
    if match is None:
        return None, None, data
    #return (match.group(1).strip(), match.group(2).strip(), data[match.end():])
    return (match.group(1).strip(), match.group(2).strip(), match.end())

def iterMKWc(data, startpos=0):
    """Iterate through MKW, contents pairs in a string.
    """
    while True:
        mKW, contents, startpos = get_next_metakw(data, startpos)
        if mKW is None:
            break
        yield mKW, contents
def iterMKWd(data, startpos=0):
    """Iterate through MKW, contents dict pairs in a string.
    """
    while True:
        mKW, contents, startpos = get_next_metakw(data, startpos)
        if mKW is None:
            break
        yield mKW, parse_kw(contents)

rKW = re.compile(r"\\(\w+)\{(.*?)\}")
def parse_kw(contents):
    """Parses the contents of mKW argument.

    return a dict of everything inside of it.
    """
    ret = { }
    while True:
        match = rKW.search(contents)
        if match is None:
            break
        ret[match.group(1).strip()] = match.group(2).strip()
        contents = contents[match.end():]
    return ret
        

if __name__ == "__main__":
    import sys
    f = sys.argv[1]
    data = file(f).read()
    
    #m = rMKW.search(data)
    kw, contents, rest = get_next_metakw(data)
    #print kw
    print contents
    print parse_kw(contents)




