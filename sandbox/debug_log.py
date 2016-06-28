import logging
# debug_logger = logging.getLogger("debug-logger")

# import coloredlogs
# coloredlogs.install(level='DEBUG')


def autolog(tag=None, message=""):
    "Automatically log the current function details."
    import inspect, logging
    # Get the previous frame in the stack, otherwise it would
    # be this function!!!
    func = inspect.currentframe().f_back.f_code
    # Dump the message + the name of this function to the log.
    log_string = "" 
    
    log_string += "\n*** "
    if tag is not None:
        log_string += "[" + tag + "] "

    log_string += " %s() in %s:%s ***" % (func.co_name, 
            func.co_filename,
            str(func.co_firstlineno))

    log_string += "\n" + str(message)

    log_string += "\n"
    print(log_string)
