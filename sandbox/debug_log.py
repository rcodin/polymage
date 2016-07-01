from colorama import init, Fore, Style
# debug_logger = logging.getLogger("debug-logger")

# import coloredlogs
# coloredlogs.install(level='DEBUG')


def header(message):
    message = str(message)
    message = "%s" % message 
    return (Fore.BLUE + "%s\n%s\n" % (message, '-' * len(message)) + 
        Style.RESET_ALL)

def autolog(message, tag=None):
    "Automatically log the current function details."
    import inspect
    # Get the previous frame in the stack, otherwise it would
    # be this function!!!
    func = inspect.currentframe().f_back.f_code
    # Dump the message + the name of this function to the log.
    log_string = "" 
    
    log_string += Fore.GREEN
    log_string += "\n*** "
    if tag is not None:
        log_string += "[" + tag + "] "

    log_string += " %s() in %s:%s" % (func.co_name, 
            func.co_filename,
            str(func.co_firstlineno))
    log_string += " ***\n"
    log_string += Style.RESET_ALL

    log_string += str(message)

    log_string += "\n"
    print(log_string)
