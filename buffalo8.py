import time
import warnings

__all__ = ['read']

_GAMEPAD = None

try:
    import hid
    
    _GAMEPAD = hid.device()
    _GAMEPAD.open(0x0583, 0x2060)
    _GAMEPAD.set_nonblocking(True)

    _LAST_STATUS = [128, 128, 0, 0, 0, 0, 0, 0]
    _LAST_PRESS = [0 for s in _LAST_STATUS]
    
except ImportError as e:
    warnings.warn("Cannot import the 'hid' module, gamepad support disabled")
    
except OSError as e:
    _GAMEPAD = None
    
    warnings.warn("Cannot find a supported device, gamepad support disabled")
    
    
def read():
    """
    Read from a connected Buffalo Two-axis, Eight button USB game controller
    and return what button is being pressed.  If nothing is being pressed return
    an empty string.
    
    .. note::  If more than one button is pressed only the last button processed
               is returned.
    """
    
    global _LAST_STATUS
    global _LAST_PRESS
    
    output = []
    
    if _GAMEPAD is not None:
        # If we have something to read from, try reading from it
        status = _GAMEPAD.read(64)
        
        if len(status):
            ## Two-axis arrow keys - left/right
            refresh = True
            if status[0] == _LAST_STATUS[0]:
                if time.time() < _LAST_PRESS[0] + 0.2:
                    refresh = False
                    
            if refresh:
                if status[0] == 0:
                    output.append('left')
                    _LAST_PRESS[0] = time.time()
                elif status[0] == 255:
                    output.append('right')
                    _LAST_PRESS[0] = time.time()
                    
            ## Two-axis arrow keys - up/down
            refresh = True
            if status[1] == _LAST_STATUS[1]:
                if time.time() < _LAST_PRESS[1] + 0.2:
                    refresh = False
                    
            if refresh:
                if status[1] == 0:
                    output.append('up')
                    _LAST_PRESS[1] = time.time()
                elif status[1] == 255:
                    output.append('down')
                    _LAST_PRESS[1] = time.time()
                    
            ## Buttons - A/B/X/Y, L/R, and start/select
            refresh = True
            if status[2] == _LAST_STATUS[2]:
                if time.time() < _LAST_PRESS[2] + 0.2:
                    refresh = False
                    
            if refresh:
                if status[2] & 1:
                    output.append('A')
                    _LAST_PRESS[2] = time.time()
                if status[2] & 2:
                    output.append('B')
                    _LAST_PRESS[2] = time.time()
                if status[2] & 4:
                    output.append('X')
                    _LAST_PRESS[2] = time.time()
                if status[2] & 8:
                    output.append('Y')
                    _LAST_PRESS[2] = time.time()
                if status[2] & 16:
                    output.append('L')
                    _LAST_PRESS[2] = time.time()
                if status[2] & 32:
                    output.append('R')
                    _LAST_PRESS[2] = time.time()
                if status[2] & 64:
                    output.append('select')
                    _LAST_PRESS[2] = time.time()
                if status[2] & 128:
                    output.append('start')
                    _LAST_PRESS[2] = time.time()
                    
            ## "Extra" buttons - turbo/clear
            refresh = True
            if status[3] == _LAST_STATUS[3]:
                if time.time() < _LAST_PRESS[3] + 0.2:
                    refresh = False
                    
            if refresh:
                if status[3] & 16:
                    output.append('turbo')
                    _LAST_PRESS[3] = time.time()
                if status[3] & 32:
                    output.append('clear')
                    _LAST_PRESS[3] = time.time()
                    
            _LAST_STATUS = status
            
    return output
