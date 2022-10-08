# Chebyshev's bias
import math
max = 107
primes = []
primes_onemore_than_multiple_of_4 = []
primes_oneless_than_multiple_of_4 = []
def main():
    current_number = 3
    
    for i in range(max):
        aprime = True
        for prev in range(2, int(math.sqrt(current_number) + 1)):
            if current_number % prev == 0: 
                aprime = False
                break

        if aprime:
            primes.append(current_number)
            if current_number % 4 == 1:
              primes_onemore_than_multiple_of_4.append(current_number)
            else:
               primes_oneless_than_multiple_of_4.append(current_number)
            diff = (len(primes_oneless_than_multiple_of_4)) - (len(primes_onemore_than_multiple_of_4))
            if i in range(0,max +1,1000):
                print(i, diff)
        current_number += 1

main()
