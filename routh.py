from pprint import pprint

## fadel 1. n check the stability of system, n count the number of poles in 
## RHS and LHS and on jw axis
## 2. handle the whole row being zero 

# -------------------- Examples --------------------

# s^5+s^4+10s^3+72s^2+152s+240 ==> in the project
# s^2+2s+1 ==> in a video on youtube
# 2s^5+3s^4+2s^3+3s^2+2s+1 ==> epsilon
# s^5+2s^4+24s^3+48s^2-25s-50 ==> whole row = 0 (unstable)
# s^5+7s^4+6s^3+42s^2+8s+56 ==> whole row = 0 (marginally stable)
# s^5+2s^4+2s^3+4s^2+11s+10
# s^4+10s^3+35s^2+50s+264 ==> 3ady

# -------------------- Global Variables --------------------
coeffs = []
powers_of_s = []
exponents = []
routh_array = []
epsilon = 1e-10

# -------------------- function definition --------------------
def parse_input(input):
    if input.find('+') != -1:
        up = input.replace('+', ' ').replace('-', ' -').split()
        print(len(up))
        for str in up:
            coeff = ''
            s_power = ''
            exponentFlag = False # not an s^ sth
            for i in str:
                if i != 's' and not exponentFlag:
                    coeff += i
                else:
                    # print("i = ", i)
                    s_power += i
                    exponentFlag = True
            if coeff == '':
                coeff = '1'
            if s_power == '':
                s_power = 's^0'
            elif s_power == 's':
                s_power = 's^1'
            coeffs.append(coeff)
            powers_of_s.append(s_power)
        for i in powers_of_s:
            i = i.replace('s^','')
            exponents.append(i)
        # print(up)
        print("coeffs = " , coeffs)
        print("powers_of_s = ", powers_of_s)
        print("exponents = ", exponents)

def initialize_routh_array(coeffs):
    first_row = []
    second_row = []
    i = 0
    for col in coeffs:
        if i%2 == 0:
            first_row.append(col)
        else:
            second_row.append(col)
        i += 1

    # #Check if the last element of first_row is zero. If not append zero
    # if first_row[-1] != 0:
    #     first_row.append(0)

    while len(second_row) < len(first_row):
        second_row.append(0)
    first_row = [float(i) for i in first_row]
    second_row = [float(i) for i in second_row]
    return first_row, second_row
        
def eval_next_row(previous, before_previous):
    # given the last 2 rows in the table, we can compute the next row
    row = []
    count = 0
    # check if a row has all elements = zero 
    for i in previous:
        if i == 0: 
            count += 1
    if count == len(previous):
        print("The whole row equals zero. So we replace it with the derivative of the previous row")
        # do something
        return [] # change this
    
    # previous row doesn't have all zeros
    # before_previous ==> c   d
    # previous ==> a    b

    a = previous[0]
    c = before_previous[0]
    
    if a == 0.0:
        print("Replacing a with epsilon")
        # replace it with very small number epsilon because we can't divide by 0
        a = epsilon
    for col in range(0, len(previous)):
        val = 0.0
        if col + 1 < len(previous):
            #
            d = before_previous[col+1]
            b = previous[col+1]
            val = (a*d-b*c)/a

        row.append(val)
    return row


# -------------------- beginning of the program --------------------   
input = input('Enter characteristic equation:\n')

parse_input(input)
first_row, second_row = initialize_routh_array(coeffs= coeffs)
routh_array.append(first_row)
routh_array.append(second_row)
for i in range(0,len(powers_of_s)-2):    
    new_row = eval_next_row(routh_array[-1], routh_array[-2])
    routh_array.append(new_row)

print("routh array = ")
pprint(routh_array)