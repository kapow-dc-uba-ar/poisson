# encoding: utf-8
'''
Created on Nov 15, 2021

@author: placiana
'''
import numpy as np
import random


class Sequence(object):
    '''
    A sequence (string) wrapper with additional info (base and name).
    
    '''
    def __init__(self, base, name, file=None, sequence=''):
        '''
        
        :param base: the sequence base
        :param name: the sequence name
        :param file: if sequence param not defined it will be read from this file
        :param sequence: a sequence as string
        '''
        self.file = file
        self.sequence = sequence
        self.base = base
        self.name = name
        if not sequence:
            if not file:
                raise AttributeError('Should pass a file or an explicit sequence in constructor')
            self.sequence = self.generate_sequence()

        self.length = len(self.sequence)

    def generate_sequence(self):
        with open(self.file) as of:
            return of.read()

    def to_file(self, filename):
        with open(filename, 'w+') as afile:
                afile.write(self.sequence)


class ThueMorseSequence(Sequence):
    def __init__(self, length):
        self.base = 2
        self.name = 'Thue-Morse'
        self.length = length
        self.sequence = self.generate_sequence()

    def generate_sequence(self):
        return ''.join([log_morse(i**2) for i in range(self.length)])


class RudinSequence(Sequence):
    def __init__(self, length):
        self.base = 2
        self.name = 'Rudin'
        self.length = length
        self.sequence = self.generate_sequence()

    def generate_sequence(self):
        return ''.join([log_rudin(i**2) for i in range(self.length)])

class FibonacciSequence(Sequence):
    def __init__(self, base, length):
        self.base = base
        self.name = f'Fibonacci b{base}'
        self.length = length
        self.sequence = self.generate_sequence()

    def generate_sequence(self):
        return fibonacci_sequence(self.length, self.base)

class RandomSequence(Sequence):
    def __init__(self, base, length):
        self.base = base
        self.name = f'Random b{base}'
        self.length = length
        self.sequence = self.generate_sequence()

    def generate_sequence(self):
        return random.choices(range(self.base), k=self.length)


# Efficient fibonacci seq in  any base
def fibonacci_sequence(length, base=10):
    have = 0
    chunks = []
    a, b = 0, 1
    while have < length:
        add = np.base_repr(a, base)
        chunks.append(add)
        have += len(add)
        a, b = b, a+b
    return ''.join(chunks)[:length]

def log_morse(pos):
    power = 1
    ans = False
    while power <= pos:
        if pos & power > 0:
            ans = not ans
        power *= 2
    return '1' if ans else '0'

def log_rudin(pos):
    power = 3
    ans = False
    while power <= pos:
        if pos & power == power:
            ans = not ans
        power *= 2
    return '0' if ans else '1'

