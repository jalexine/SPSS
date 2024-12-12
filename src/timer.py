import time

# This class allows measuring the time spent by a function
class Timer:
    # Function called right before the function execution
    def __enter__(self):
        self.t1 = time.perf_counter()  # Start time
        return self

    # Function called right after the function execution
    def __exit__(self, type, value, traceback):
        self.t2 = time.perf_counter()  # End time
        self.t = self.t2 - self.t1  # Calculate elapsed time

    # Function that prints the elapsed time in the terminal
    def print(self, template: str = "{}"):
        print(template.format(round(self.t, 2)))  # Print the time rounded to 2 decimal places
