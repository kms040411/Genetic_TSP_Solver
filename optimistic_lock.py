import itertools
import threading

class Atomic_int():
    def __init__(self):
        self.__num = 0
        self.__lock = threading.Lock()

    def get(self) -> int:
        result = None
        with self.__lock:
            result = self.__num
            return result
    
    def store(self, a) -> bool:
        with self.__lock:
            self.__num = a
            return True

    def inc(self) -> int:
        with self.__lock:
            self.__num += 1
            return self.__num
    
    def compare_and_inc(self, old) -> bool:
        with self.__lock:
            if self.__num != old:
                return False
            self.__num += 1
            return True

class Optimistic_lock():
    def __init__(self):
        self.current_seq = Atomic_int()

    def read_lock(self) -> int:
        while True:
            a = self.current_seq.get()
            if a % 2 != 0:
                continue
            return a

    # Do not use
    def upgrade(self, seq) -> bool:
        if seq % 2 != 0:
            print("upgrade: seq should be even")
            exit(-1)
        if self.current_seq.compare_and_inc(seq):
            return True
        return False

    def write_lock(self) -> int:
        while True:
            a = self.current_seq.get()
            if a % 2 == 0 and self.current_seq.compare_and_inc(a):
                return a

    def validate(self, seq) -> bool:
        return seq == self.current_seq.get()

    def unlock(self, seq):
        self.current_seq.store(seq + 2)

