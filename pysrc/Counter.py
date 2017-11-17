# This is a wrapper for an int that allows the next empty index of an array to be tracked when the array is being filled by recursive functions
class Counter:
    def __init__(self, N):
        self.N = N