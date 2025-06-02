import math

def lmd(I_y, area, length, c=1):
    r = math.sqrt(I_y/area) 
    lmd = (c*length)/r
    return round(lmd, 2)

#Running test on all functions 
if __name__ == '__main__':
    print(lmd(I_y=79820.37, area=646, length=600, c=1), "and expected: 53.97")