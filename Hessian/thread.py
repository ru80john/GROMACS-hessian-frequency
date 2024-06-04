import threading
import time
import sys


def rotating_spinner(done_event):
    chars = ['-', '\\', '|', '/']
    i = 0

    while not done_event.is_set():
        sys.stdout.write('\r' + chars[i % len(chars)])
        sys.stdout.flush()
        time.sleep(0.1)
        i += 1
