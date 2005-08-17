"""This module is designed to unify mp3's internal logging[1].  It may be
used by other packages as well if you want to unify things.

To use, you must first import the logger:
>>> from mp3.log import mp3log

Then, to log a message simply call the proper function:

>>> mp3log.debug("Hi, mom.")
>>> mp3log.info("Jacked up and good to go, sir.")
>>> mp3log.warn("You've not enough minerals")
>>> mp3log.error("The hive cluster is under attack!")
>>> mp3log.critical("He's dead, Jim.")

Look above for the different log levels.  They should be self
explanatory.

To set the verbosity level, use the function set_level("level") in
this module.


[1] So far, not all of mp3 is logged very well.

For more info, look at the source code for this module as well as
Python documentation for the"logging" module.  Really, the Python
logging module is good enough that mp3.log isn't needed, but it
provides a nice way to document what I'm doing.

"""
#
# This module is designed to unify the logging for mp3.
#
# mp3log is the logger instance.

import sys
import logging

default_log_level = "warn"   # one of "debug", "info", "warn", "error", "critical"

mp3log = logging.getLogger('mp3')
#_theloghandler = logging.StreamHandler(sys.stdout)  # leave blank for stderr
_theloghandler = logging.StreamHandler()  # leave blank for stderr
_theloghandler.setFormatter( logging.Formatter("mp3 %(filename)-15s %(levelname)-8s %(message)s"))
mp3log.addHandler(_theloghandler)
mp3log.debug('Logger Initilized')
mp3log.info('MP3 Library 0.0')
mp3log.debug('Old joke: never use version 0.x or x.0 (UT ITS UNIX Group)')
#thelog.setLevel(logging.INFO)

mp3log = mp3log

# These are for compatability with the "old" way of doing things.
#
# Example:
# >>> import mp3.log
# >>> mp3.log.debug("\strong{What?}")
#
# They used to be thin wrapper functions around mp3log.METHOD, but
# that caused problems because it looked like the log messages were
# generated by this file, not by where they really came from!  I guess
# if I need that kind of control again, I'll subclass the logger.

debug    = mp3log.debug       # don't use these.
info     = mp3log.info    
warn     = mp3log.warn    
error    = mp3log.error   
critical = mp3log.critical 


def set_level(level):
    """Set the desired logging level.

    Set the logging lever for the master mp3 logger.  Argument should
    be one of "debug", "info", "warn", "error", or "critical" to set
    the logger to that level.
    """
    level = level.lower()
    if level == "debug":
        debug("setting log level to debug.")
        mp3log.setLevel(logging.DEBUG)
    elif level == "info":
        debug("setting log level to info.")
        mp3log.setLevel(logging.WARNING)
    elif level[0:4] == "warn":
        debug("setting log level to warn.")
        mp3log.setLevel(logging.WARN)
    elif level == "error":
        debug("setting log level to error.")
        mp3log.setLevel(logging.ERROR)
    elif level == "critical":
        debug("setting log level to critical.")
        mp3log.setLevel(logging.CRITICAL)
    else:
        mp3log.warn("could not set log level to '%s' (spelling, maybe?)"%level)

set_level(default_log_level)

    